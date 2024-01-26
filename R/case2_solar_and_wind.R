# --------------------------------------
# Case 2: Solar and Offshore Wind farms
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
PV <- read_csv("data/SolarPV/nhh_Timeseries_60.423_5.302_E5_1kWp_crystSi_14_47deg_6deg_2005_2020.csv",
                skip = 10) %>% 
  mutate(datetime = as.POSIXct(time, format = "%Y%m%d:%H%M")) %>% 
  mutate(CF = P/1000, locID = "Solar PV") %>% # CF = percentage of 1 kW
  select(datetime, locID, CF)%>% 
  filter(year(datetime) %in% 2005:2019)
  


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
    port %>% mutate(locID = "WP portfolio"),
    PV %>% filter(!is.na(datetime))
) %>% 
  as_tsibble(key = locID, index = datetime)

power <- power %>%  
  index_by(date =~as.Date(.)) %>% 
  group_by(locID) %>% 
  summarize(CF = mean(CF, na.rm =T)) %>% 
  filter(year(date) %in% 2005:2019) %>% 
  mutate(locID = factor(locID, levels = c(1:20, "WP portfolio","Solar PV")))

fitted_seasons <- power %>% 
  model(TSLM(CF ~ fourier(K = 3, period = "year"))) %>% 
  fitted() %>% 
  filter(year(date) %in% 2014:2016) %>% 
  select(-.model)  
fig1 <- ggplot(data = fitted_seasons %>% filter(locID %in% 1:20), aes(x =date, y = .fitted, group = locID)) + 
  geom_line(color = "grey", alpha = .5) +
  geom_line(data=fitted_seasons%>% filter(!(locID %in% 1:20)), aes(x =date, y = .fitted, color = locID),
            lwd = .9)+
  geom_point(data=power %>% filter(!(locID %in% 1:20),
                                   year(date) %in% 2014:2016),
             aes(x= date, y = CF, color = locID), size = .2)+
  scale_y_continuous(name = "Capacity factor", limits =c(-0.01,1.01), expand = c(0,0))+
  scale_x_date(expand = c(0,0), limits = c(as.Date("2013-12-16"),as.Date("2016-12-31")), 
               name = "") +
  scale_color_manual(values = c("blue", "red"))+
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        legend.margin=margin(-20, 0, 0, 0),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
ggsave(fig1, file = "figures/case2_wind_solar_time_plot_season.png", width = 8, height = 4)

wppPV <- power %>% 
  filter(!(locID %in% 1:20))
empcov <- wppPV %>% 
  as_tibble() %>% 
  pivot_wider(names_from = locID, values_from = CF) %>% select(-1) %>%
  cov()
empcor <- wppPV %>% 
  as_tibble() %>% 
  pivot_wider(names_from = locID, values_from = CF) %>% select(-1) %>%
  cor()
season_parameters <- wppPV %>%   
  model(TSLM(CF ~ fourier(K = 3, period = "year"))) %>% 
  coef() %>%select(locID, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)
names(season_parameters)[-1] <- c("intercept", paste0(c("a","b"), rep(1:((ncol(season_parameters)-2)/2), each = 2)))
season_cov <- season_parameters[,-(1:2)] %>% as.matrix() %*%t(season_parameters[,-(1:2)] %>% as.matrix()) /2
colnames(season_cov)<-rownames(season_cov)<-season_parameters$locID
season_cov <- season_cov[2:1,2:1]

diag(1/sqrt(diag(season_cov)))%*%season_cov%*%diag(1/sqrt(diag(season_cov)))
diag(1/sqrt(diag(empcov)))%*%empcov%*%diag(1/sqrt(diag(empcov)))
diag(1/sqrt(diag(empcov-season_cov)))%*%(empcov-season_cov)%*%diag(1/sqrt(diag(empcov-season_cov)))



mu <- season_parameters %>%  arrange(desc(locID)) %>% pull(intercept)

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
              MPT, cov = empcov, mu = mu)
season <- sapply(targets,
              MPT, cov = season_cov, mu = mu)
seasonadjusted <- sapply(targets,
                 MPT, cov = empcov - season_cov, mu = mu)
ports_by_cov <- tibble(
  targets = targets, 
  "Empirical" = emp, 
  "Seasonal" = season,
  "Season-adj." = seasonadjusted
) %>% 
  pivot_longer(cols = 2:4) %>% 
  mutate(name = factor(name, levels = c("Empirical", "Seasonal", "Season-adj.")))
  
fig2 <- ports_by_cov %>% 
  ggplot(aes(x = value, y = targets, color = name, group = name)) + 
  geom_path() + 
  scale_x_continuous("Variance measure")+
  scale_y_continuous("Portfolio capacity factor")+
  geom_segment(data = ports_by_cov %>% group_by(name) %>% 
               filter(value == min(value)),
               aes(x = -Inf, xend = value, y = targets, yend= targets, color = name),
               arrow = arrow(length = unit(0.02, "npc")))+
  geom_vline(xintercept = 0, lty = 2, color = "grey")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        legend.title = element_blank(), 
        legend.position = c(.7,.15),
        legend.text  = element_text(size = 12),
        legend.background = element_rect(fill = "transparent", color = "transparent"))
ggsave(fig2, file = "figures/case2_efficient_frontiers.png", width = 8, height = 4)

targets[c(which.min(emp), which.min(season), which.min(seasonadjusted))]
weights <- (MPTsummary <- rbind(
  "Empirical" = MPT(mu,empcov, target = targets[which.min(emp)], return_value = FALSE),
  "Seasonal" = MPT(mu,season_cov, target = targets[which.min(season)], return_value = FALSE),
  "Season-adjusted" = MPT(mu,empcov-season_cov, target = targets[which.min(seasonadjusted)], return_value = FALSE)) %>% 
  
  as.data.frame() %>% 
  rownames_to_column("covariance") %>% 
  as_tibble() %>% 
  mutate(covariance = factor(covariance, levels = c("Empirical", "Seasonal", "Season-adjusted")))) %>% 
  rename("Solar PV" = "weights1", 
         "WP portfolio" = "weights2") %>% 
  select(-c(2:3)) %>% 
  pivot_longer(cols = 2:3, names_to = "locID", values_to = "weights") %>% 
  mutate(locID = factor(locID, levels = c("Solar PV" ,"WP portfolio")))

portfolios_ts <- wppPV %>% as_tibble() %>% left_join(weights, by = "locID", relationship = "many-to-many") %>% 
  group_by(date, covariance) %>% 
  summarize(totalCF = sum(weights*CF)) %>% 
  as_tsibble(key = covariance, index = date) %>% 
  
  mutate(covariance = factor(covariance, levels = c("Empirical", "Seasonal", "Season-adjusted")))



SD <- function(w, sigma){
  w <- c(w,1-w)
  return(sqrt(t(w)%*%sigma%*%w))
}
SRS <- function(w, season, sigma){
  w <- c(w,1-w)
  return(sqrt(t(w)%*%season%*%w / t(w)%*%sigma%*%w))
}

(timeplot <- portfolios_ts %>% 
  filter(year(date) %in% 2014:2019) %>% 
  ggplot(aes(x=date, y = totalCF)) + 
  geom_line()+
  facet_wrap(~covariance, ncol = 1) +
  geom_text(data= weights %>% filter(locID == "Solar PV"),
            aes(x = as.Date("2017-01-01"), y = Inf, 
                label = paste0("Solar weight: ", round(100*weights,1),"%    ",
                "SRS: ",round(sapply(weights,SRS, season= season_cov,sigma= empcov),3),
                "    SD: ", round(sapply(weights, SD, sigma = empcov),3))),
            vjust = 1.5, hjust=.5)+
  scale_x_date(expand = c(0,0), date_breaks = "1 year", date_labels = "%Y")+
  scale_y_continuous(name = "Portfolio capacity factor", limits = c(0.09,.58), breaks = seq(.1,.6,.1))+
  theme(axis.title.x = element_blank())+
  geom_hline(data=MPTsummary, aes(yintercept = mean), lty = 2, col = 2))

acfplot <- portfolios_ts %>% 
  ACF(totalCF, lag_max = 400) %>%
  autoplot()  +
  scale_x_cf_lag(breaks = seq(0, 400, 30), limits = c(0,401), expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")
fig3 <- ggpubr::ggarrange(timeplot+ theme(strip.text =element_blank()),
                  acfplot+ theme(axis.title.x = element_blank()),
                  labels = c("D","E"))
ggsave(fig3, file = "figures/case2_Time_ACF_combined.png", width = 10, height = 8)


SRSfig <- tibble(
  w = seq(0,1,0.01),
  SRS = sapply(w, SRS, season = season_cov, sigma = empcov)
) %>% 
  ggplot(aes(x=w, y = SRS)) + geom_path()+
  geom_vline(data = weights %>% filter(locID == "Solar PV"),
             aes(xintercept=weights, color = covariance), lty = 2)+
  ggtext::geom_richtext(data = weights %>% filter(locID == "Solar PV"),
             aes(x=weights,y=0, color = covariance, label = covariance),angle=90,
             show.legend = FALSE, hjust=-0.01)+
  geom_segment(data = weights %>% filter(locID == "Solar PV"),
               aes(x=-Inf, xend = weights, 
                   y = sapply(weights,SRS, season= season_cov,sigma= empcov), 
                   yend = sapply(weights,SRS, season= season_cov,sigma= empcov),
                   color = covariance), lty = 2)+
  scale_x_continuous(expand = c(0,0), name = "Solar PV weight",
                     breaks= seq(0,1,.25),
                     labels = c(seq(0,.75,.25),"1.0"))+
  scale_y_continuous(expand = c(0,0), name = "Seasonal risk score",
                     limits = c(0,.75))+
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent"))
ggsave("figures/case2_SRS.png", width = 8, height = 6)


library(ggpubr)
case2<-ggarrange(ggarrange(fig1,SRSfig,fig2,ncol = 3, labels = c("A","B","C")),
          fig3, ncol = 1, heights = c(1,2))
ggsave(case2, file = "figures/case2_combined.pdf", width = 8, height = 10)
