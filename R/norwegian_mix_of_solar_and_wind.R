# --------------------------------------
# Case 2: Solar and Wind offshore farms
# --------------------------------------
# Wind data from: 
#   https://github.com/holleland/OffshoreWindPortfolios/data
# Solar data from: 
#   https://joint-research-centre.ec.europa.eu/photovoltaic-geographical-information-system-pvgis/pvgis-tools/hourly-radiation_en
# --------------------------------------
rm(list=ls())
library(tidyverse)
library(fpp3)
library(lubridate)
library(quadprog)
theme_set(theme_bw()+
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  strip.background = element_rect(fill = "transparent",  color = "transparent")))
# Solar data: 
solar_files <- list.files("data/solar/", full.names = TRUE)
PV <- tibble()
for(file in solar_files){
  PV <- 
    bind_rows(PV,
              read_csv(file, skip = 10) %>% 
                mutate(datetime = as.POSIXct(time, format = "%Y%m%d:%H%M")) %>% 
                mutate(CF = P/1000, locID = "Solar PV") %>% # CF = percentage of 1 kW
                select(datetime, locID, CF)%>% 
                filter(year(datetime) %in% 2005:2019) %>% 
                mutate(file = file))
}
#PV <- apply(solar_files, read_csv, skip = 10)
# PV <- read_csv("data/Timeseries_60.423_5.302_E5_1kWp_crystSi_14_47deg_6deg_2005_2020.csv",
#                 skip = 10) %>% 
#   mutate(datetime = as.POSIXct(time, format = "%Y%m%d:%H%M")) %>% 
#   mutate(CF = P/1000, locID = "Solar PV") %>% # CF = percentage of 1 kW
#   select(datetime, locID, CF)%>% 
#   filter(year(datetime) %in% 2005:2019)
PV %>% group_by(file) %>% summarize(mean(CF))

PV <- PV %>%  
  filter(!grepl("Trondheim|Tromso|Bodo", file)) %>% 
  group_by(datetime, locID) %>% 
  summarize(CF = mean(CF, na.rm=T))

#   
  

# NVE data and portfolio
wind <-   readRDS("data/NVE.rds") %>%
                      select(locID,value, datetime)%>% 
  filter(year(datetime) %in% 2005:2019)


port <- tibble(
  locID = 1:20,
  turbines = c(269, 146,100,208,204,41,102,194,37,0,0,63,0,25,0,0,218,32,28,333)) %>% 
  # From Table A1 of Hølleland et al (2023)
  left_join(wind, by = "locID") %>%
  group_by(datetime) %>% 
  summarize(CF = sum(turbines*value)/2000)
power <- wind %>% 
  rename(CF = value) %>% 
  mutate(locID = as.character(locID)) %>% 
  bind_rows(
    port %>% mutate(locID = "Wind offshore"),
    PV %>% filter(!is.na(datetime))
) %>% 
  as_tsibble(key = locID, index = datetime)

power <- power %>%  
  index_by(date =~as.Date(.)) %>% 
  group_by(locID) %>% 
  summarize(CF = mean(CF, na.rm =T)) %>% 
  filter(year(date) %in% 2005:2019) %>% 
  mutate(locID = factor(locID, levels = c(1:20, "Wind offshore","Solar PV")))


# Estimate matrices :
Amat <- (WP <- power %>% as_tibble() %>% 
           filter(locID %in% c("Wind offshore","Solar PV")) %>% 
           pivot_wider(names_from = locID, values_from = CF)) %>% 
  forecast::msts(seasonal.periods = 365.25) %>%
  forecast::fourier(K = 4)
Y <- WP %>% select(-date) %>% as.matrix()

# Selecting K = 4
# Use K = 5 above and run the following models: 
aic_solar <- AIC(lm(Y[,1]~1+Amat),
    lm(Y[,1]~1+Amat[,-(ncol(Amat)-1:0)]),
    lm(Y[,1]~1+Amat[,-(ncol(Amat)-3:0)]),
    lm(Y[,1]~1+Amat[,-(ncol(Amat)-5:0)]))
aic_wind <- AIC(lm(Y[,2]~1+Amat),
    lm(Y[,2]~1+Amat[,-(ncol(Amat)-1:0)]),
    lm(Y[,2]~1+Amat[,-(ncol(Amat)-3:0)]),
    lm(Y[,2]~1+Amat[,-(ncol(Amat)-5:0)]))
rownames(aic_solar)<- paste0("K=",5:2)
names(aic_solar)<- c("df", "Solar")
aic_solar$Wind = aic_wind$AIC
aic_solar$Solar <- aic_solar$Solar-min(aic_solar$Solar)
aic_solar$Wind <- aic_solar$Wind-min(aic_solar$Wind)
t(round(aic_solar[4:1,-1],1)) %>% xtable::xtable()

ffit <- lm(Y~ 1 + Amat)

A <- t(coef(ffit)[-1, ])

A <- A[,c(seq(2,ncol(A),2),seq(1,ncol(A)-1,2))]

SigS <- nrow(WP)/(nrow(WP)-1) * A%*%t(A)/2
SigZ <- cov(residuals(ffit))
Sig  <- cov(Y)

cov2cor(SigS)
cov2cor(SigZ)
cov2cor(Sig)
# store fitted values for purpose of plotting:
print(Sig*1e2) %>% round(2)
print(SigS*1e2) %>% round(2)
print(SigZ*1e2)%>% round(2)


.fitted <- fitted(ffit)
fitted.df <- left_join(reshape2::melt(Y) %>% rename("CF"=value), 
                       reshape2::melt(.fitted) %>% 
                         rename(fitted = value)) %>% 
  mutate(date = WP$date[Var1]) %>% 
  rename(locID= Var2) %>% 
  select(-Var1) %>% 
  as_tibble()


# 
# 
# fig1 <- fitted.df %>% 
#   mutate(locID2 = ifelse(locID == "Wind offshore", "Offshore Wind", "Solar PV")) %>% 
#   dplyr::select(-locID) %>% 
#   rename(locID=locID2) %>% 
#   ggplot(aes(x =date, y = fitted, color = locID))+
#   geom_line(lwd = .7)+
#   geom_point(aes(x= date, y = CF, color = locID), size = .2, show.legend = FALSE)+
#   scale_y_continuous(name = "Capacity factor", limits =c(-0.01,1.01), expand = c(0,0))+
#   scale_x_date(expand = c(0,0), limits = c(as.Date("2004-12-16"),as.Date("2019-12-31")),
#                date_breaks = "1 year",
#                date_labels = "%Y",
#                name = "") +
#   scale_color_manual(values = c( "blue","red"))+
#   guides(color = guide_legend(override.aes = list(size = 2, lwd =5)))+
#   theme(legend.title = element_blank(), 
#         legend.position = "bottom",
#         legend.margin=margin(-20, 0, 0, 0),
#         legend.background = element_rect(fill = "transparent", color = "transparent"))
# fig1
# ggsave(fig1, file = "figures/case2_wind_solar_time_plot_season.pdf", width = 6, height = 3)



wppPV <- power %>% 
  filter(!(locID %in% 1:20))
empcov <- wppPV %>% 
  as_tibble() %>% 
  pivot_wider(names_from = locID, values_from = CF) %>% select(-1) %>%
  cov()
empcor <- Y %>% 
  cor()

diag(1/sqrt(diag(SigS)))%*%SigS%*%diag(1/sqrt(diag(SigS)))
diag(1/sqrt(diag(Sig)))%*%Sig%*%diag(1/sqrt(diag(Sig)))
diag(1/sqrt(diag(Sig-SigS)))%*%(Sig-SigS)%*%diag(1/sqrt(diag(Sig-SigS)))
diag(1/sqrt(diag(SigZ)))%*%SigZ%*%diag(1/sqrt(diag(SigZ)))


mu <- coef(ffit)[1,]

MPT <- function(mu, cov, target, return_value = TRUE){
  Amat <- cbind(1, mu, diag(length(mu)))
  bvec <- c(1, target, rep(0,length(mu)))
  sol <- solve.QP(Dmat = cov, 
            dvec = rep(0, length(mu)),
            Amat= Amat,
            bvec = bvec, 
            meq = 2)
  if(return_value) return(sol$value)
  return(c("mean"=target, "value" = sol$value, "weights" = sol$solution))
}
targets <- seq(min(mu)+0.01, max(mu)-0.01, by = 0.01)
emp <- sapply(targets,
              MPT, cov = Sig, mu = mu)
season <- sapply(targets,
              MPT, cov = SigS, mu = mu)
seasonadjusted <- sapply(targets,
                 MPT, cov = SigZ, mu = mu)
ports_by_cov <- tibble(
  targets = targets, 
  "Balanced" = emp, 
  "Seasonal" = season,
  "Season-adj." = seasonadjusted
) %>% 
  pivot_longer(cols = 2:4) %>% 
  mutate(name = factor(name, levels = c("Season-adj.", "Balanced", "Seasonal")))
  
fig2 <- ports_by_cov %>% 
  ggplot(aes(x = value, y = targets, color = name, group = name)) + 
  geom_path() + 
  scale_x_continuous("Risk")+
  scale_y_continuous(latex2exp::TeX("Portfolio capacity factor ($\\mu^*$)"))+
  geom_segment(data = ports_by_cov %>% group_by(name) %>% 
               filter(value == min(value)),
               aes(x = -Inf, xend = value, y = targets, yend= targets, color = name),
               arrow = arrow(length = unit(0.02, "npc")))+
  geom_vline(xintercept = 0, lty = 2, color = "grey")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.title = element_blank(), 
        legend.position.inside = TRUE,
        legend.position = c(.65,.15),
        legend.text  = element_text(size = 12),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
ggsave(fig2, file = "figures/case2_efficient_frontiers.png", width = 8, height = 4)

targets[c(which.min(emp), which.min(season), which.min(seasonadjusted))]
weights <- (MPTsummary <- rbind(
  "Balanced" = MPT(mu,Sig, target = targets[which.min(emp)], return_value = FALSE),
  "Seasonal" = MPT(mu,SigS, target = targets[which.min(season)], return_value = FALSE),
  "Season-adjusted" = MPT(mu,SigZ, target = targets[which.min(seasonadjusted)], return_value = FALSE)) %>% 
  
  as.data.frame() %>% 
  rownames_to_column("covariance") %>% 
  as_tibble() %>% 
  mutate(covariance = factor(covariance, levels = c("Season-adjusted", "Balanced", "Seasonal")))) %>% 
  rename("Solar PV" = "weights1", 
         "Wind offshore" = "weights2") %>% 
  select(-c(2:3)) %>% 
  pivot_longer(cols = 2:3, names_to = "locID", values_to = "weights") %>% 
  mutate(locID = factor(locID, levels = c("Solar PV" ,"Wind offshore")))

portfolios_ts <- wppPV %>% as_tibble() %>% left_join(weights, by = "locID", relationship = "many-to-many") %>% 
  group_by(date, covariance) %>% 
  summarize(totalCF = sum(weights*CF)) %>% 
  as_tsibble(key = covariance, index = date) %>% 
  
  mutate(covariance = factor(covariance, levels = c("Season-adjusted", "Balanced", "Seasonal")))



SD <- function(w, sigma){
  w <- c(w,1-w)
  return((t(w)%*%sigma%*%w))
}
SRS <- function(w, season, sigmaZ){
  w <- c(w,1-w)
  return((t(w)%*%season%*%w / t(w)%*%(season+sigmaZ)%*%w))
}

(timeplot <- portfolios_ts %>% 
  filter(year(date) %in% 2014:2019) %>% 
  ggplot(aes(x=date, y = totalCF, color = covariance)) + 
  geom_line()+
  facet_wrap(~covariance, ncol = 1) +
  geom_text(data= weights %>% filter(locID == "Solar PV"),
            aes(x = as.Date("2017-01-01"), y = Inf, 
                label = paste0("Solar weight: ", round(100*weights,1),"%    ",
                "SRS: ",round(sapply(weights,SRS, season= SigS,sigmaZ= SigZ),3),
                "    Var: ", round(sapply(weights, SD, sigma = Sig),3))),
            vjust = 1.5, hjust=.5)+
  scale_x_date(expand = c(0,0), date_breaks = "1 year", date_labels = "%Y")+
  scale_y_continuous(name = "Portfolio capacity factor", limits = c(0.09,.58), breaks = seq(.1,.6,.1))+
  theme(axis.title.x = element_blank(),
        legend.position = "none")+
  geom_hline(data=MPTsummary, aes(yintercept = mean), lty = 2, col = "black"))

# acfplot <- portfolios_ts %>% 
#   ACF(totalCF, lag_max = 400) %>%
#   autoplot(color = covariance)  +
#   scale_x_cf_lag(breaks = seq(0, 400, 30), limits = c(0,401), expand = c(0,0))+
#   scale_y_continuous(name = "Autocorrelation")

acf_data <- portfolios_ts %>% 
  ACF(totalCF, lag_max = 400) %>% 
  as_tibble()

conf_level <- 0.05
conf_bound <- qnorm(1 - conf_level / 2) / sqrt(nrow(portfolios_ts))
# Plot using ggplot2
acfplot <- ggplot(acf_data, aes(x = lag, y = acf, fill = covariance)) +
  geom_col(width = 1) +
  #scale_fill_manual(values = c("positive" = "blue", "negative" = "red")) +
  scale_x_continuous(breaks = seq(0, 400, 30), limits = c(0, 401), expand = c(0, 0)) +
  scale_y_continuous(name = "Autocorrelation") +
  facet_wrap(~covariance, ncol = 1, strip.position = "right")+
  geom_hline(yintercept = c(-conf_bound, conf_bound),
             linetype = "dashed", color = "black")+
  theme(legend.position = "none")

fig3 <- ggpubr::ggarrange(timeplot+ theme(strip.text =element_blank()),
                  acfplot+ theme(axis.title.x = element_blank()),
                  labels = c("C","D"))
ggsave(fig3, file = "figures/case2_Time_ACF_combined.png", width = 10, height = 8)




weights$covariance <- as.character(weights$covariance)
weights$covariance[5:6] <- "Season-adj."
weights$covariance <- factor(weights$covariance, 
                                levels = c("Season-adj.","Balanced", "Seasonal"
                                           ))
SRSfig <- tibble(
  w = seq(0,1,0.01),
  SRS = sapply(w, SRS, season = SigS, sigmaZ = SigZ)
) %>% 
  ggplot(aes(x=w, y = SRS)) + geom_path()+
  geom_vline(data = weights %>% filter(locID == "Solar PV"),
             aes(xintercept=weights, color = covariance), lty = 2)+
  ggtext::geom_richtext(data = weights %>% filter(locID == "Solar PV"),
             aes(x=weights,y=Inf, color = covariance, label = covariance),angle=90, 
             show.legend = FALSE, hjust=1)+#-0.01)+
  geom_segment(data = weights %>% filter(locID == "Solar PV"),
               aes(x=-Inf, xend = weights, 
                   y = sapply(weights,SRS, season= SigS,sigmaZ= SigZ), 
                   yend = sapply(weights,SRS, season= SigS,sigmaZ= SigZ),
                   color = covariance), lty = 2)+
  scale_x_continuous(expand = c(0,0), name = "Solar PV weight",
                     breaks= seq(0,1,.25),
                     labels = c(seq(0,.75,.25),"1.0"))+
  scale_y_continuous(expand = c(0,0), name = "Seasonal risk score",
                     limits = c(0,.85))+
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent"))
ggsave(SRSfig, file = "figures/case2_SRS.png", width = 8, height = 6)




library(ggpubr)
case2<-ggarrange(ggarrange(SRSfig,fig2,ncol = 2, labels = c("A","B")),
          fig3, ncol = 1, heights = c(1,2))
ggsave(case2, file = "figures/case2_combined.pdf", width = 10, height = 10)


par(mfrow=c(1,2))
qqnorm(ffit$residuals[,1], main = "Solar PV")
qqline(ffit$residuals[,1], col = 2)
qqnorm(ffit$residuals[,2], main = "Wind portfolio")
qqline(ffit$residuals[,2], col = 2)

resid <- ffit$residuals
sfig1 <- resid %>% reshape2::melt() %>% 
  mutate(Var3 = ifelse(Var2 == "Wind offshore", "Offshore Wind", "Solar PV")) %>% 
  dplyr::select(-Var2) %>% 
  rename(Var2=Var3) %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ggplot(aes(sample=standardized)) +
  stat_qq() + 
  stat_qq_line()+
  facet_wrap(~Var2, ncol = 1, strip.position = "left")+
  xlab("Theoretical quantiles")+
  ylab("Residual quantiles")+
  theme(strip.placement = "outside")
ggsave("figures/case2_appendix_qqplots_per_energysource.png", width = 8, height = 4)
#resid <- ffit2$residuals
# resid %>% reshape2::melt() %>% 
#   group_by(Var2) %>% 
#   mutate(standardized = (value-mean(value))/sd(value)) %>% 
#   ggplot(aes(sample=standardized)) +
#   stat_qq() + 
#   stat_qq_line()+
#   facet_wrap(~Var2)+
#   xlab("Theoretical quantiles")+
#   ylab("Residual quantiles")
# ggsave("figures/case2_appendix_qqplots_per_energysource_transformed.png", width = 8, height = 4)

sfig2 <- resid %>% reshape2::melt() %>% 
  mutate(Var3 = ifelse(Var2 == "Wind offshore", "Offshore Wind", "Solar PV")) %>% 
  dplyr::select(-Var2) %>% 
  rename(Var2=Var3) %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ggplot(aes(x=standardized)) +
  stat_density(fill = "skyblue")+
  geom_function(fun = dnorm, col = 2, lwd = 1.2)+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.48))+
  scale_x_continuous(expand = c(0,0), limits = c(-3.3,2.8))+
  facet_wrap(~Var2, ncol = 1, strip.position = "left")+
  xlab("Standardized residuals")+
  ylab("Density")+
  theme(strip.placement = "outside",
        strip.text = element_blank())
  
sfig3 <- resid %>% reshape2::melt() %>% 
  mutate(Var3 = ifelse(Var2 == "Wind offshore", "Offshore Wind", "Solar PV")) %>% 
  dplyr::select(-Var2) %>% 
  rename(Var2=Var3) %>% 
  as_tsibble(key = Var2, index = Var1) %>% 
  ACF(value, lag_max = 400) %>% 
  autoplot() +
  scale_x_cf_lag(breaks = seq(0, 400, 100), limits = c(0,400), expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")+
  theme(axis.title.x = element_blank())+
  facet_wrap(~Var2, ncol = 1)+
  xlab("lag")+
  theme(strip.text = element_blank())
ggpubr::ggarrange(sfig1,sfig2,sfig3, ncol = 3,
                  labels = c("A","B","C"))
ggsave("figures/case2_appendix_residuals.png", width = 10, height = 10*1/2)
ggsave("figures/case2_appendix_residuals.pdf", width = 12, height = 8)


# --------------
# Consquences
# --------------
weights_wide <- weights %>% 
  pivot_wider(names_from = locID, values_from = weights) %>% 
  mutate(windpower = 30,
         solarpower = windpower * mu[2]/mu[1]*`Solar PV`/`Wind offshore`,
         area = solarpower*1e3 / 200,
         if_30p_CF_on_solar = windpower * mu[2]/.3*`Solar PV`/`Wind offshore`,)

weights_wide


test <- wppPV %>% as_tibble() %>% left_join(weights, by = "locID", relationship = "many-to-many") %>% filter(covariance == "Seasonal") %>% 
  mutate(weights = weights[1]) %>% 
  pivot_wider(names_from = "locID", values_from = "CF") %>% 
  mutate(Production = `Wind offshore`*30+30*mu[2]/mu[1]*weights/(1-weights))
test %>% 
  ggplot(aes(x=date, y = Production)) + geom_line()


test <- wppPV %>% as_tibble() %>% left_join(weights, by = "locID", relationship = "many-to-many") %>% filter(covariance == "Balanced") %>% 
  mutate(weights = weights[1]) %>% 
  pivot_wider(names_from = "locID", values_from = "CF") %>% 
  mutate(Production = `Wind offshore`*30+30*mu[2]/mu[1]*weights/(1-weights))
test %>% 
  ggplot(aes(x=date, y = Production)) + geom_line()

test <- wppPV %>% as_tibble() %>% left_join(weights, by = "locID", relationship = "many-to-many") %>% filter(covariance == "Season-adj.") %>% 
  mutate(weights = weights[1]) %>% 
  pivot_wider(names_from = "locID", values_from = "CF") %>% 
  mutate(Production = `Wind offshore`*30+30*mu[2]/mu[1]*weights/(1-weights))
test %>% 
  ggplot(aes(x=date, y = Production)) + geom_line()


summary(ffit)
weights %>%  filter(locID == "Solar PV") %>% 
  mutate(install = 30*mu[2]/mu[1]*weights/(1-weights))

spectrum(Y[,2], method = "ar")
fft_result <- fft(Y[,1])
power <- Mod(fft_result)^2
n <- length(Y[,1])
freq <- (0:(n - 1)) / n

# Plot only up to Nyquist frequency
half_n <- floor(n / 500)
plot(freq[1:half_n], power[1:half_n], type = "l",
     xlab = "Frequency", ylab = "Power", main = "Power Spectrum")

x  <- Y/1e3
#x <- x
n <- nrow(x)
fft_result <- fft(x)
power <- Mod(fft_result)^2 / n  # Normalize
# Frequencies (in cycles per sample)
freq <- (0:(n - 1)) / n
j = 0:21; 
q = 15
k = j*q
# Plot only up to Nyquist frequency
half_n <- floor(n / 2)
tibble(j = j,
       "Solar PV" = power[k+1,1],
       "Wind offshore" = power[k+1,2]) %>%
  filter(j>0) %>% 
  pivot_longer(cols = c(`Solar PV`, `Wind offshore`)) %>% 
  ggplot(aes(x=j, y = (value), color = name, fill = name)) + 
  geom_bar(stat="identity", position = "dodge", width = .9)+
  #scale_y_continuous(name = "Spectral Power of Fourier frequencies", expand = c(0,0), limits = c(0, 58))+
  scale_x_continuous(name = "k", breaks = 1:8,expand = c(0,0))+
  theme(legend.title = element_blank(),
        legend.position.inside = TRUE,
        legend.position = c(.9,.9))
ggsave("figures/power_spectrum_of_K.pdf",
       width = 8, height = 4)
plot(j, log(power[k+1]), type = "l",
     xlab = "Frequency", ylab = "Power", main = "Power Spectrum")

# To get dominant frequencies:
dominant_freqs <- n*(freq)[order(power, decreasing = TRUE)][1:50]

k <- 1:half_n
k %% 365
tibble(
  k = 1:half_n,
  freq = freq[k],
  power = power[1:half_n],
  k_trun = (k %% 365)+1
) %>% 
  group_by(k_trun) %>% 
  summarize(power = mean(power)) %>% 
  ggplot(aes(x=k_trun, y = power)) + geom_point()
which.max(power[1:half_n])



targets <- seq(min(mu)+0.01, max(mu)-0.01, by = 0.01)
emp <- sapply(targets,
              MPT, cov = Sig, mu = mu, return_value = FALSE) %>% t() %>% 
  as_tibble() %>% 
  mutate(strategy = "Balanced")
season <- sapply(targets,
                 MPT, cov = SigS, mu = mu, return_value = FALSE) %>% t() %>% 
  as_tibble() %>% 
  mutate(strategy = "Seasonal")
seasonadjusted <- sapply(targets,
                         MPT, cov = SigZ, mu = mu, return_value = FALSE) %>% t() %>% 
  as_tibble() %>% 
  mutate(strategy = "Season-adjusted")

weights_full <- bind_rows(
  emp, 
  season,
  seasonadjusted
)
weights_full %>% 
  ggplot(aes(x=weights1, y = value, col = strategy)) + geom_line()+
  scale_y_continuous("Risk")+
  scale_x_continuous("Solar weight")+
  geom_vline(aes(xintercept = weights1, col = strategy),
             lty = 2,
             data = weights_full %>% group_by(strategy) %>% 
               filter(value == min(value)))


weights_full %>% 
  ggplot(aes(x=weights1, y = mean, col = strategy)) + geom_line()+
  scale_y_continuous("Risk")+
  scale_x_continuous("Solar weight")


# Power spectrum stuff:
X <- Y[,1:2]
pow<- matrix(NA_real_, ncol = 2, nrow=10)
for(k in 1:(nrow(pow))){
  pow[k,] <- Mod(colMeans(X*exp(-complex(1,0,1)*2*pi/365.25*k*(1:length(X)))))^2
}
colnames(pow) <- colnames(Y)
colnames(pow)[2] <- "Wind"
tibble(k=1:nrow(pow)) %>% 
  bind_cols(pow) %>% 
  pivot_longer(cols = 2:3) %>% 
  group_by(name) %>% 
  mutate(cumsum = cumsum(value)/sum(value)) %>% 
  ggplot(aes(x=k, y=cumsum, color = name)) +
  #geom_line()+
  geom_point() + geom_line(show.legend=FALSE)+
  geom_hline(yintercept=0.99, color = "steelblue", linewidth = .7,lty=1)+
  geom_vline(xintercept = 4, lty = 2, color = "steelblue")+
  scale_x_continuous(breaks =seq(1,10,1),
                     limits = c(.5,10.5),
                     expand = c(0,0))+
  scale_y_continuous(labels = scales::percent, 
                     breaks = c(seq(.9,1,.05),.99),
                     name = expression(P[k](omega)),
                     limits =c(.9, 1),
                     expand =c(0,0))+
  
  scale_color_manual(values = c("red","blue"))+
  theme(legend.title = element_blank(),
        legend.position.inside = T,
        legend.position = c(.85,.1),
        legend.background = element_rect(fill = "transparent",
                                         color = "transparent"))
ggsave("figures/power_spectrum_of_K_energymix.pdf", width = 5, height = 4)  
