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
  index_by(date = ~as.Date(.)) %>% 
  summarize(CF = mean(CF))
windpower %>% ungroup() %>% 
  filter(year(date) == 2014) %>% 
  autoplot(CF)

season_models <- windpower %>% 
  model(TSLM(CF ~fourier(K=3, period = "1 year"))) 
season_parameters <- season_models %>% 
  coef() %>% 
  select(locID, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)
names(season_parameters) <- c("locID", "intercept", "a","b")
names(season_parameters)[-1] <- c("intercept", paste0(c("a","b"), rep(1:((ncol(season_parameters)-2)/2), each = 2)))
A <- season_parameters[,-(1:2)] %>% as.matrix()
season_cov <- A%*% t(A)/2
rownames(season_cov)<-colnames(season_cov)<- 1:20
fig1<-season_models %>% 
  fitted() %>%  
  filter(year(date) %in% 2014:2016) %>% 
  select(-.model) %>% 
  group_by(locID) %>% 
  ggplot(aes(x=date, y = .fitted, color = factor(locID))) +geom_line()+
  geom_point(data = windpower %>% filter(year(date) %in% 2014:2016),
             aes(x = date, y = CF, color = factor(locID)), alpha = .5, size = .1)+
  scale_y_continuous(name = "Capacity factor", 
                     limits = c(-.01,1.01), 
                     expand =c(0,0))+
  scale_x_date(limits = c(as.Date("2013-12-31"), as.Date("2017-01-01")), expand = c(0,0),
               date_breaks = "6 months", date_labels = "%b-%Y")+
  scale_color_discrete(name = "Wind farm\nlocations")+
  theme(axis.title.x = element_blank())
ggsave(fig1, file = "figures/case1_time_plots.png", width = 15, height = 7)


empcov <- windpower %>% 
  as_tibble() %>% 
  pivot_wider(names_from = locID, values_from = CF) %>% 
  select(-date) %>% 
  cov()
empcor <- windpower %>% 
  as_tibble() %>% 
  pivot_wider(names_from = locID, values_from = CF) %>% 
  select(-date) %>% 
  cor()
season_cor <- diag(1/sqrt(diag(season_cov))) %*% season_cov %*% diag(1/sqrt(diag(season_cov)))
rownames(season_cor)<-colnames(season_cor)<- 1:20
covariances <- left_join(
empcov %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
  pivot_longer(cols =  2:21, names_to = "locID2", values_to = "Empirical covariance"),
empcor %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
  pivot_longer(cols =  2:21, names_to = "locID2", values_to = "Empirical correlation")) %>% 
  left_join(
    season_cov %>%  as_tibble() %>% rownames_to_column("locID1") %>% 
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


sqrt(emp.turbines %*% season_cov %*%emp.turbines /
emp.turbines %*% empcov %*%emp.turbines )

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
ggsave("figures/case1_correlation_matrices.png", width = 8, height = 4)

covariances %>% 
  filter(type == .type) %>% 
  filter(locID1 %in% c(1,5,8,12,14,20),
         locID2 %in% c(1,5,8,12,14,20)) %>% 
  ggplot(aes(x=locID1, y = locID2, fill = covr, color = covr)) + geom_tile() +
  scale_y_reverse(expand = c(0,0), breaks =seq(1,20,1))+
  scale_x_continuous(expand =c(0,0), breaks =seq(1,20,1))+
  scale_fill_viridis_c(name = .type)+
  scale_color_viridis_c(name = .type)+
  facet_wrap(~source, ncol = 2)+
  guides(fill = guide_colorbar(barheight = 15))

det(season_cov[c(1,5,8,12,14,20),c(1,5,8,12,14,20)])


# Mean value: 
mu <- season_parameters %>%  arrange(locID) %>% pull(intercept)
Amat <- cbind(1, mu, diag(length(mu)),-diag(length(mu)))
bvec <- c(1, .6, rep(0,length(mu)),-maxWeights)
(emp.turbines <- solve.QP(Dmat=empcov, 
         dvec = rep(0, length(mu)),
         Amat= Amat,
         bvec = bvec, 
         meq = 2)$solution * 2000) %>% round(digits = 0)

include <- c(1, 4, 5, 8, 20)
include <- c(1,5,8,12,14, 20)
(emp.turbines <- solve.QP(Dmat=empcov[include,include], 
                          dvec = rep(0, length(mu[include])),
                          Amat= Amat[include,c(1:2,2+include, 22+include)],
                          bvec = c(1, .6, rep(0,length(include)),-maxWeights[include]), 
                          meq = 2)$solution * 2000) %>% round(digits = 0)
#(include <- sample(1:20, replace = FALSE, size = 7)
#)
emp.turbines %>% round(0)

(season.turbines <- solve.QP(Dmat=season_cov[include,include], 
          dvec = rep(0, length(include)),
          Amat= Amat[include,c(1:2,2+include, 22+include)],
          bvec = c(1, .6, rep(0,length(include)),-maxWeights[include]), 
          meq = 2)$solution * 2000) %>% round(digits = 0)
(season.turbines.unrestricted <- solve.QP(Dmat=season_cov[include,include], 
                             dvec = rep(0, length(include)),
                             Amat= Amat[include,c(1:2,2+include)],
                             bvec = c(1, .6, rep(0,length(include))), 
                             meq = 2)$solution * 2000) %>% round(digits = 0)

  (seasonadj.turbines <- solve.QP(Dmat=empcov[include,include]-season_cov[include,include], 
                             dvec = rep(0, length(include)),
                             Amat= Amat[include,c(1:2,2+include, 22+include)],
                             bvec = c(1, .6, rep(0,length(include)),-maxWeights[include]), 
                             meq = 2)$solution * 2000) %>% round(digits = 0)

(solve.QP(Dmat=season_cov, 
          dvec = rep(0, length(mu)),
          Amat= Amat,
          bvec = c(1, .6, rep(0,length(mu))), 
          meq = 2)$solution * 2000) %>% round(digits = 0)

portfolios <- (weights <- tibble(locID = include,
       "Seasonal" = round(season.turbines,0),
       "Seasonal (unrestricted)" = round(season.turbines.unrestricted,0),
       "Empirical" = round(emp.turbines,0),
       "Season-adjusted"=round(seasonadj.turbines,0))) %>% 
  right_join(windpower, by = "locID") %>% 
  pivot_longer(cols = 2:5, names_to = "source", values_to = "turbines") %>% 
  group_by(date, source) %>% 
  summarize(CFport= sum(turbines*CF,na.rm=T)/2000) %>% 
  mutate(source = factor(source, 
                         levels = c("Empirical", "Seasonal","Season-adjusted", "Seasonal (unrestricted)")))
SRS <- function(w, season, sigma){
  #w <- c(w,1-w)
  return(as.numeric(sqrt(t(w)%*%season%*%w / t(w)%*%sigma%*%w)))
}
srs <- tibble(
  "Seasonal"= SRS(weights$Seasonal, 
                  season = season_cov[include,include],
                  sigma = empcov[include,include]),
  "Empirical"= SRS(weights$Empirical, 
                  season = season_cov[include,include],
                  sigma = empcov[include,include]),
  "Season-adjusted"= SRS(weights$`Season-adjusted`, 
                  season = season_cov[include,include],
                  sigma = empcov[include,include]),
  "Seasonal (unrestricted)"= SRS(weights$`Seasonal (unrestricted)`, 
                  season = season_cov[include,include],
                  sigma = empcov[include,include])) %>% 
  pivot_longer(cols = 1:4,
               names_to = "source", values_to = "srs") %>% 
  mutate(source = factor(source, 
                         levels = c("Empirical","Seasonal",
                                    "Season-adjusted", "Seasonal (unrestricted)"))) %>% 
  right_join(
    portfolios %>% group_by(source) %>% 
      summarize(sd = sd(CFport))
  )


porttimeplots <- portfolios %>% 
  filter(year(date)%in% 2014:2017) %>% 
  ggplot(aes(x=date, y = CFport)) + geom_line() +
  scale_y_continuous(name = "Portfolio capacity factor")+
  facet_wrap(~source, ncol = 1, strip.position = "right")+
  theme(axis.title.x = element_blank(),
        strip.placement = "none")+
  geom_hline(yintercept = .6, lty = 2, col = 2)
(portACFplots <- portfolios %>% 
  as_tsibble(key = source, index = date) %>% 
  ACF(CFport, lag_max = 400) %>% 
  autoplot() +
  scale_x_cf_lag(breaks = seq(0, 400, 50), limits = c(0,401), expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")+
  theme(axis.title.x = element_blank())+
  geom_text(data = srs, aes(x = Inf, y = Inf, label = paste("SRS: ", round(srs,3),"\nSD: ",round(sd,3))),
            hjust = 1, vjust = 1.2))
figz<-ggpubr::ggarrange(porttimeplots+theme(strip.text =element_blank()),portACFplots, ncol = 2, labels = c("C","D"))
ggsave(file = "figures/case1_Time_ACF_combined.png", width = 10, height = 8)


res <- weights %>% 
  mutate("Max turbines" = round(maxWeights[include]*2000,0)) %>% 
  relocate(locID, `Max turbines`, Empirical, Seasonal, `Season-adjusted`)
res <- as_tibble(cbind(nms = names(res), t(res))) %>% 
  mutate("srs" = c(NA_real_,NA_real_,srs$srs[c(2,1,3,4)])) 
#  mutate("sd"  = c(NA_real_,NA_real_, srs$sd))
res$sd<-  c(NA_real_,NA_real_, srs$sd[c(2,1,3,4)])

library(xtable)
res %>% xtable() %>% print(include.rownames = F)


# Combined figure:

library(ggpubr)
combo <- ggarrange(fig1+theme(legend.position = "none"),
                   fig2+theme(legend.background = element_rect(fill = "transparent", color = "transparent")),
                   figz, ncol = 1, heights = c(1,1.9,3.5), labels = c("A","B",NA))
ggsave(combo, file = "figures/case1_all_combined.pdf", width = 8, height = 11)
