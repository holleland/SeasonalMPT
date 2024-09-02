# Case study
# Holleland et al (2023)
# Data available at: "https://github.com/holleland/OffshoreWindPortfolios/raw/main/data/
rm(list=ls())
library(tidyverse)
library(quadprog)
library(fpp3)
library(lubridate)
theme_set(theme_bw()+
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  strip.background = element_rect(fill = "transparent",  color = "transparent")))
maxWeights <- readxl::read_excel("data/NVE_areal.xlsx")$maxWeight
windpower <- readRDS(file = "data/NVE.rds") %>% 
  select(datetime, value, locID, lon,lat) %>% 
  rename(CF = value) %>%
  as_tsibble(key = locID, 
             index = datetime)
windpower <- windpower %>% 
  group_by(locID) %>% 
  index_by(date = ~yearweek(.)) %>% 
  summarize(CF = mean(CF))
windpower %>% ungroup() %>% 
  filter(year(date) == 2014) %>% 
  autoplot(CF)

# Estimate matrices: 

Amat <- (WP <- windpower %>% as_tibble() %>% 
  pivot_wider(names_from = locID, values_from = CF)) %>% 
  forecast::msts(seasonal.periods = 52) %>%
  forecast::fourier(K = 10)
Y <- WP %>% select(-date) %>% as.matrix()
ffit <- lm(Y~ 1 + Amat)
# -- Could apply transformation: 
# Y2 <- Y
# Y2[which(Y2==0, arr.ind=T)] <- 0.001
# Y2[which(Y2==1, arr.ind=T)] <- 0.999
# Ytrans <- log(Y2/(1-Y2))
# ffit <- lm(Ytrans~ 1 + Amat)



A <- t(coef(ffit)[-1, ])

A <- A[,c(seq(2,ncol(A),2),seq(1,ncol(A)-1,2))]

SigS <- nrow(WP)/(nrow(WP)-1) * A%*%t(A)/2
SigZ <- cov(residuals(ffit))
Sig  <- .5*SigS+.5*SigZ

# store fitted values for purpose of plotting:
.fitted <- fitted(ffit)
fitted.df <- left_join(reshape2::melt(Y) %>% rename("CF"=value), 
          reshape2::melt(.fitted) %>% 
            rename(fitted = value)) %>% 
  mutate(date = WP$date[Var1]) %>% 
  rename(locID= Var2) %>% 
  select(-Var1) %>% 
  as_tibble()


# ----------------
(matlib::gaussianElimination(A, diag(nrow(A))))[1:20,1:20]

fig1 <- fitted.df %>% 
  filter(year(date) %in% c(2013:2016)) %>% 
  #select(-.model) %>% 
 # group_by(locID) %>% 
  ggplot(aes(x=date, y = fitted, color = factor(locID))) +geom_line()+
  geom_point(aes(x = date, y = CF, color = factor(locID)), alpha = .5, size = .1)+
  scale_y_continuous(name = "Capacity factor", 
                     limits = c(-.01,1.01), 
                     expand =c(0,0))+
   scale_x_yearweek(limits = c(as.Date("2012-12-25"), as.Date("2016-12-31")), 
                    expand = c(0,0),
                    date_breaks = "year")+
  scale_color_discrete(name = "Wind farm\nlocations")+
  theme(axis.title.x = element_blank())
ggsave(fig1, file = "figures/case1_time_plots_weekly.png", width = 15, height = 7)


empcor <- Y %>% 
  cor()
season_cor <- diag(1/sqrt(diag(SigS))) %*% SigS %*% diag(1/sqrt(diag(SigS)))
rownames(season_cor)<-colnames(season_cor)<- 1:20
covariances <- left_join(
Sig %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
  pivot_longer(cols =  2:21, names_to = "locID2", values_to = "Empirical covariance"),
empcor %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
  pivot_longer(cols =  2:21, names_to = "locID2", values_to = "Empirical correlation")) %>% 
  left_join(
    SigS %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
      pivot_longer(cols =  2:21, names_to = "locID2", values_to = "Seasonal covariance")
  ) %>% 
  left_join(
    season_cor %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
      pivot_longer(cols =  2:21, names_to = "locID2", values_to = "Seasonal correlation")
  ) %>% 
  pivot_longer(cols = 3:6, names_to = "source", values_to = "covr") %>% 
  mutate(type = ifelse(grepl("covariance", source, fixed = TRUE), "Covariance", "Correlation")) %>% 
  mutate(locID1 = as.numeric(locID1),
         locID2 = as.numeric(locID2))


#sqrt(emp.turbines %*% season_cov %*%emp.turbines /
#emp.turbines %*% empcov %*%emp.turbines )

.type <- "Correlation"
fig2 <- covariances %>% 
  filter(type == .type) %>% 
  ggplot(aes(x=locID1, y = locID2, fill = covr, color = covr)) + geom_tile() +
  scale_y_reverse(name = "Location ID", expand = c(0,0), breaks =seq(1,20,1))+
  scale_x_continuous(name = "Location ID", expand =c(0,0), breaks =seq(1,20,1))+
  scale_fill_viridis_c(name = .type)+
  scale_color_viridis_c(name = .type)+
  facet_wrap(~source, ncol = 2)+
  guides(fill = guide_colorbar(barheight = 12))
ggsave(fig2, file = "figures/case1_correlation_matrices_weekly.png", width = 8, height = 4)

# covariances %>% 
#   filter(type == .type) %>% 
#   filter(locID1 %in% c(1,5,8,12,14,20),
#          locID2 %in% c(1,5,8,12,14,20)) %>% 
#   ggplot(aes(x=locID1, y = locID2, fill = covr, color = covr)) + geom_tile() +
#   scale_y_reverse(expand = c(0,0), breaks =seq(1,20,1))+
#   scale_x_continuous(expand =c(0,0), breaks =seq(1,20,1))+
#   scale_fill_viridis_c(name = .type)+
#   scale_color_viridis_c(name = .type)+
#   facet_wrap(~source, ncol = 2)+
#   guides(fill = guide_colorbar(barheight = 15))
# 
# det(season_cov[c(1,5,8,12,14,20),c(1,5,8,12,14,20)])


# Mean value: 
mustar <- .6
mu <- coef(ffit)[1,]
#mu <- season_parameters %>%  arrange(locID) %>% pull(intercept)
Amat <- cbind(1, mu, diag(length(mu)),-diag(length(mu)))
bvec <- c(1, mustar, rep(0,length(mu)),-maxWeights)
(emp.turbines <- solve.QP(Dmat=Sig, 
         dvec = rep(0, length(mu)),
         Amat= Amat,
         bvec = bvec, 
         meq = 2)$solution * 2000) %>% round(digits = 0)

include <- 1:20
(emp.turbines <- solve.QP(Dmat=Sig[include,include], 
                          dvec = rep(0, length(mu[include])),
                          Amat= Amat[include,c(1:2,2+include, 22+include)],
                          bvec = c(1, mustar, rep(0,length(include)),-maxWeights[include]), 
                          meq = 2)$solution * 2000) %>% round(digits = 0)
#(include <- sample(1:20, replace = FALSE, size = 7)
#)
emp.turbines %>% round(0)

(season.turbines <- solve.QP(Dmat=SigS[include,include], 
          dvec = rep(0, length(include)),
          Amat= Amat[include,c(1:2,2+include, 22+include)],
          bvec = c(1, mustar, rep(0,length(include)),-maxWeights[include]), 
          meq = 2)$solution * 2000) %>% round(digits = 0)
(season.turbines.unrestricted <- solve.QP(Dmat=SigS[include,include], 
                             dvec = rep(0, length(include)),
                             Amat= Amat[include,c(1:2,2+include)],
                             bvec = c(1, mustar, rep(0,length(include))), 
                             meq = 2)$solution * 2000) %>% round(digits = 0)

  (seasonadj.turbines <- solve.QP(Dmat=SigZ[include,include], 
                             dvec = rep(0, length(include)),
                             Amat= Amat[include,c(1:2,2+include, 22+include)],
                             bvec = c(1, mustar, rep(0,length(include)),-maxWeights[include]), 
                             meq = 2)$solution * 2000) %>% round(digits = 0)


portfolios <- (weights <- tibble(locID = include,
       "Seasonal" = round(season.turbines,0),
       "Seasonal (unrestricted)" = round(season.turbines.unrestricted,0),
       "Naive" = round(emp.turbines,0),
       "Season-adjusted"=round(seasonadj.turbines,0))) %>% 
  right_join(windpower, by = "locID") %>% 
  # If transformation: 
 # mutate(CF = ifelse(CF == 0, log(0.001/.999), ifelse(CF == 1, -log(0.001/.999), log(CF/(1-CF))))) %>% 
  pivot_longer(cols = 2:5, names_to = "source", values_to = "turbines") %>% 
  group_by(date, source) %>% 
  summarize(CFport= sum(turbines*CF,na.rm=T)/2000) %>% 
  mutate(source = factor(source, 
                         levels = c("Naive", "Seasonal","Season-adjusted", "Seasonal (unrestricted)")))
SRS <- function(w, season, sigmaZ){
  #w <- c(w,1-w)
  return(as.numeric(sqrt(t(w)%*%season%*%w / t(w)%*%(season+sigmaZ)%*%w)))
}
srs <- tibble(
  "Seasonal"= SRS(weights$Seasonal, 
                  season = SigS[include,include],
                  sigmaZ = SigZ[include,include]),
  "Naive"= SRS(weights$Naive, 
                  season = SigS[include,include],
                  sigmaZ = SigZ[include,include]),
  "Season-adjusted"= SRS(weights$`Season-adjusted`, 
                  season = SigS[include,include],
                  sigma = Sig[include,include]),
  "Seasonal (unrestricted)"= SRS(weights$`Seasonal (unrestricted)`, 
                  season = SigS[include,include],
                  sigmaZ = SigZ[include,include])) %>% 
  pivot_longer(cols = 1:4,
               names_to = "source", values_to = "srs") %>% 
  mutate(source = factor(source, 
                         levels = c("Naive","Seasonal",
                                    "Season-adjusted", "Seasonal (unrestricted)"))) %>% 
  right_join(
    portfolios %>% group_by(source) %>% 
      summarize(sd = sd(CFport))
  ) %>% 
  arrange(source)
srs$SSAC400=portfolios %>% as_tsibble(index = date, key = source) %>% 
  ACF(CFport, lag_max = 52) %>% as_tibble() %>% 
  group_by(source) %>% 
  summarize("SSAC400"=sum(acf^2)) %>% 
  arrange(source) %>% 
  pull(SSAC400)
  
dates <- seq(as.Date("2012-01-01"),as.Date("2018-01-01"), by = "1 year")
porttimeplots <- portfolios %>% 
  filter(year(date)%in% 2013:2017) %>% 
  ggplot(aes(x=date, y = CFport)) + geom_line() +
  scale_y_continuous(name = "Portfolio capacity factor")+
  scale_x_yearweek(expand = c(0,0), date_breaks = "1 year")+
  facet_wrap(~source, ncol = 1, strip.position = "right")+
  theme(axis.title.x = element_blank(),
        strip.placement = "none")+
  geom_hline(yintercept = .6, lty = 2, col = 2)
(portACFplots <- portfolios %>% 
  as_tsibble(key = source, index = date) %>% 
  ACF(CFport, lag_max = 52) %>% 
  autoplot() +
  scale_x_cf_lag(breaks = seq(0, 52, 10), limits = c(0,53), expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")+
  theme(axis.title.x = element_blank())+
  geom_text(data = srs, aes(x = 26, y = Inf, label = paste("SRS: ", round(srs,3),"\tSD: ",round(sd,3),
                                                            "\n",
                            "SS52: ", round(SSAC400,3))),
            hjust = 0.5, vjust = 1.5))
figz<-ggpubr::ggarrange(porttimeplots+theme(strip.text =element_blank()),portACFplots, ncol = 2, labels = c("C","D"))
ggsave(figz, file = "figures/case1_Time_ACF_combined_weekly.png", width = 10, height = 8)


res <- weights %>% 
  mutate("Max turbines" = round(maxWeights[include]*2000,0)) %>% 
  relocate(locID, `Max turbines`, Naive, Seasonal, `Season-adjusted`)
res <- as_tibble(cbind(nms = names(res), t(res))) %>% 
  mutate("srs" = c(NA_real_,NA_real_,srs$srs)) 
#  mutate("sd"  = c(NA_real_,NA_real_, srs$sd))
res$sd<-  c(NA_real_,NA_real_, srs$sd)

library(xtable)
res[,c(1, 22:23, 2:10)] %>% xtable() %>% print(include.rownames = F, include.colnames = FALSE)
res[,c(1,11:21)] %>% xtable() %>% print(include.rownames = F, include.colnames = FALSE)

# Combined figure:

library(ggpubr)
combo <- ggarrange(fig1+theme(legend.position = "none"),
                   fig2+theme(legend.background = element_rect(fill = "transparent", color = "transparent")),
                   figz, ncol = 1, heights = c(1,1.9,3.5), labels = c("A","B",NA))
ggsave(combo, file = "figures/case1_all_combined_weekly.pdf", width = 8, height = 11)


portfolios
library(TSA)
par(mfrow=c(2,2))
periodogram(portfolios %>% filter(source=="Naive") %>% pull(CFport),
            main = "Naive", ylim = c(0,38))
periodogram(portfolios %>% filter(source=="Seasonal") %>% pull(CFport), 
            main = "Seasonal", ylim = c(0,38))
periodogram(portfolios %>% filter(source=="Season-adjusted") %>% pull(CFport), 
            main = "Seasonal adjust", ylim = c(0,38))
periodogram(portfolios %>% filter(source=="Seasonal (unrestricted)") %>% pull(CFport), 
            main = "Seasonal (unrestricted)", ylim = c(0,38))

weights %>% 
  pivot_longer(cols = 2:5) %>% 
  group_by(locID) %>% 
  filter(any(value>0)) %>% 
  ungroup() %>% 
  #)       name != "Seasonal (unrestricted)") %>% 
  mutate(name = factor(name, 
                       levels = c("Naive", "Seasonal", "Season-adjusted", "Seasonal (unrestricted)"))) %>% 
  ggplot(aes(x=factor(locID), y = value, fill = name, group=name))+
  geom_bar(stat="identity", position = "dodge")+
  scale_y_continuous(name = "Number of turbines", limits = c(0,1500), breaks = seq(0,1500,100), expand=c(0,0))+
  xlab("Non-zero-turbine locations")+
  theme(legend.position = c(.9,.85),
        legend.background = element_rect(fill = "transparent",color = "transparent"),
        legend.title =element_blank(),
        legend.direction = "vertical")
ggsave("figures/case1_turbine_locations_weekly.png", width = 9, height = 5)

library(sf)

proj <-  "+proj=lcc +lat_0=66.3 +lon_0=-42 +lat_1=66.3 +lat_2=66.3 +no_defs +R=6.371e+06"
  locations <-  readRDS(file = "data/NVE.rds") %>% 
    select(locID, lon,lat) %>% distinct()
  
land <- rnaturalearth::ne_countries(country = "Norway", 
                                      returnclass = "sf", scale = 10) %>%
    st_transform(crs=4326) %>% 
    select(geometry)
locsSF <- st_as_sf(locations, coords = c("lon","lat"), crs= 4326)  %>% st_transform(crs=4326)
locsSF$lon <- st_coordinates(locsSF)[,"X"]
locsSF$lat <- st_coordinates(locsSF)[,"Y"]
land %>%
  ggplot() + 
  geom_sf(pch = 15,fill = "transparent")+#fill = "wheat") +
  coord_sf(expand = FALSE, xlim = range(locsSF$lon)+c(-1,1), 
           ylim = range(locsSF$lat)++c(-1,1)) + 
  
  geom_text(data = locsSF, aes(x=lon,y=lat,label=locID), color = "black", size = 6)+
  theme(axis.title =element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.border=element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.background = element_rect(fill = "transparent", color = "transparent"))
ggsave("figures/norway_locations_transparentMap_week.png", bg = "transparent")



resid <- ffit$residuals
resid %>% reshape2::melt() %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ggplot(aes(sample=standardized)) +
  stat_qq() + 
  stat_qq_line()+
  facet_wrap(~Var2)+
  xlab("Theoretical quantiles")+
  ylab("Residual quantiles")
ggsave("figures/case1_appendix_qqplots_per_region_weekly.png", width = 10*.8,
       height = 8*.8)

resid %>% reshape2::melt() %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ggplot(aes(x=standardized)) +
  stat_density(fill = "skyblue")+
  geom_function(fun = dnorm, col = 2, lwd = 1.2)+
  facet_wrap(~Var2)
ggsave("figures/case1_appendix_density_plots_residuals.png", width = 8, height =6.4)
resid %>% reshape2::melt() %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ACF(x=standardized) +
  stat_density(fill = "skyblue")+
  geom_function(fun = dnorm, col = 2, lwd = 1.2)+
  facet_wrap(~Var2)
resid %>% reshape2::melt() %>% 
  as_tsibble(key = Var2, index = Var1) %>% 
  ACF(value, lag_max = 52) %>% 
  autoplot() +
  scale_x_cf_lag(breaks = seq(0, 52, 10), limits = c(0,53), expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")+
  theme(axis.title.x = element_blank())+
  facet_wrap(~Var2, ncol = 5)
ggsave("figures/case1_appendix_ACF_plots_residuals.png", width = 8, height =6.4)



appfig1 <- resid %>% reshape2::melt() %>% 
  filter(Var2 %in% seq(1, 20, by = 5)) %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ggplot(aes(sample=standardized)) +
  stat_qq() + 
  stat_qq_line()+
  facet_wrap(~Var2, ncol = 5)+
  xlab("Theoretical quantiles")+
  ylab("Residual quantiles")

appfig2 <- resid %>% reshape2::melt() %>% 
  filter(Var2 %in% seq(1, 20, by = 5)) %>% 
  group_by(Var2) %>% 
  mutate(standardized = (value-mean(value))/sd(value)) %>% 
  ggplot(aes(x=standardized)) +
  xlab("Standardized residuals")+
  ylab("Density")+
  stat_density(aes(fill = "skyblue" ),color = "black")+
  geom_function(fun = dnorm, col = 2, lwd = 1.2)+
  facet_wrap(~Var2, ncol = 5)+
  theme(strip.text = element_blank(),
        legend.background = element_rect(fill = "grey70", color = "grey70"))


appfig3 <- resid %>% reshape2::melt() %>% 
  filter(Var2 %in% seq(1, 20, by = 5)) %>% 
  as_tsibble(key = Var2, index = Var1) %>% 
  ACF(value, lag_max = 52) %>% 
  autoplot() +
  scale_x_cf_lag(name = "Lag (weeks)", 
                 breaks = seq(0, 52, 10),
                 limits = c(0,53), 
                 expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")+
  theme(strip.text = element_blank())+
  facet_wrap(~Var2, ncol = 5)
ggpubr::ggarrange(appfig1,appfig2,appfig3, ncol = 1,
                  labels = c("A","B","C"))
ggsave("figures/case1_appendix_residual_plots.png", width = 8, height = 6)

resid %>% reshape2::melt() %>% 
  filter(Var2 %in% seq(1, 20, by = 5)) %>% 
  as_tsibble(key = Var2, index = Var1) %>% 
  PACF(value) %>% 
  autoplot() +facet_wrap(~Var2)+
  theme(legend.background = element_rect(fill = "grey", color = "grey"))
