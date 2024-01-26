# Case 3
rm(list=ls())
library(tidyverse)
library(readxl)
library(fpp3)
library(quadprog)
theme_set(theme_bw()+
            theme(strip.background = element_rect(fill = "white",
                                                  color = "transparent"),
                  panel.grid.minor = element_blank()))
Total <- read_excel("data/08403_20231127-125058.xlsx", range = "B4:W457") %>% 
  mutate(purpose = "total")
Conferance <- read_excel("data/08403_20231127-125058.xlsx", range = "X4:AR457") %>% 
  mutate(purpose = "conferance",
         yearmonth =unlist(Total[,1]))
Occupation <- read_excel("data/08403_20231127-125058.xlsx", range = "AS4:BM457") %>% 
  mutate(purpose = "occupation",
         yearmonth = unlist(Total[,1]))
Holiday <- read_excel("data/08403_20231127-125058.xlsx", range = "BN4:CH457") %>% 
  mutate(purpose = "holiday",
         yearmonth = unlist(Total[,1]))
names(Total)[1] <- "yearmonth"
population <- read_excel("data/01222_20231127-140229.xlsx", range = "A4:V108")
tourism <- bind_rows(Total, Conferance, Occupation, Holiday)
names(tourism) <- 
  c("yearmonth", 
    "01 Ostfold", 
    "02 Akershus",
    "03 Oslo",
    "04 Hedmark",
    "05 Oppland",
    "06 Buskerud",
    "07 Vestfold",
    "08 Telemark",
    "09 Aust-Agder",
    "10 Vest-Agder",
    "11 Rogaland",
    "12 Hordaland",
    "14 Sogn og Fjordane",
    "15 More og Romsdal",
    "50 Trondelag",
    "16 Sor-Trondelag",
    "17 Nord-Trondelag",
    "18 Nordland", "19 Troms",
    "20 Finnmark",
    "21 Svalbard",
    "purpose")
names(tourism)

unique(tourism$purpose)

names(population) <- 
  c("yearquarter", 
    "01 Ostfold", 
    "02 Akershus",
    "03 Oslo",
    "04 Hedmark",
    "05 Oppland",
    "06 Buskerud",
    "07 Vestfold",
    "08 Telemark",
    "09 Aust-Agder",
    "10 Vest-Agder",
    "11 Rogaland",
    "12 Hordaland",
    "14 Sogn og Fjordane",
    "15 More og Romsdal",
    "50 Trondelag",
    "16 Sor-Trondelag",
    "17 Nord-Trondelag",
    "18 Nordland", "19 Troms",
    "20 Finnmark",
    "21 Svalbard")
names(population)




tourism <- tourism %>% #
  mutate(
    yearmonth =
      tsibble::make_yearmonth(year = as.numeric(substr(yearmonth, 0,4)),
                              month = as.numeric(substr(yearmonth, 6,7)))) %>% 
# filter(!is.na(yearmonth)) %>% 
  pivot_longer(names_to = "county", cols = 2:22, values_to = "guestnights") %>% 
  filter(!(county %in% c("50 Trondelag", "21 Svalbard")))
tourism %>% 
  filter(year(yearmonth)<=2017) %>% 
  as_tsibble(key = c(purpose, county), index =yearmonth) %>% 
  autoplot() +
  theme(legend.position = "none")

population <- population %>% 
  filter(!is.na(yearquarter)) %>% 
  mutate(yearquarter = tsibble::make_yearquarter(
    year = as.numeric(substr(yearquarter, 1,4)),
    quarter = as.numeric(substr(yearquarter, 6,6)))) %>% 
  pivot_longer(names_to = "county", cols = 2:22, values_to = "population") 
data <- left_join(tourism %>% mutate(yearquarter = yearquarter(yearmonth)), 
                  population, 
                  by= c("county", "yearquarter")) %>% 
  select(-yearquarter) %>% 
  filter(!is.na(population),!is.na(guestnights)) 

data <- data %>% as_tibble() %>% 
  mutate(countycode = as.numeric(substr(county, 1,2)),
         countrypart = ifelse(countycode <= 8, "East",
                              ifelse(countycode <= 10, "South",
                                     ifelse(countycode <= 15, "West",
                                            ifelse(countycode <= 17, "Mid",
                                                   "North"))))) #%>% 
  # group_by(purpose, countrypart, yearmonth) %>% 
  # summarize(guestnights = sum(guestnights),
  #           population = sum(population))

data <- data %>%  
  filter(yearmonth <= yearmonth("2017 Sep")) %>% 
  as_tsibble(index = yearmonth,
             key = c(county, purpose))
data %>% 
  filter(year(yearmonth)%in% 2008,
         purpose == "total"
         ) %>% 
  
  autoplot(guestnights/population*1000) 
unique(data$county)
data %>% 
  filter(#year(yearmonth) <= 2015,
         county %in% c("03 Oslo"),
         #"12 Hordaland",
         purpose != "total"
  ) %>% 
  
  ggplot(aes(x = yearmonth, 
             y = (guestnights/population),
             col = purpose, 
             linetype = county)) + geom_line()
oslo <- data %>% 
  filter(
    county %in% c("03 Oslo"),
    !(purpose %in% c("total"))#, "holiday"))
  ) %>% 
  mutate(holiday = ifelse(purpose == "holiday", "pleasure", "buisness")) %>% 
  as_tibble() %>% 
  group_by(yearmonth, holiday, county) %>% 
  summarize(population = population[1], 
            guestnights = sum(guestnights)) %>% 
  select(yearmonth, guestnights, population, holiday, county) %>% 
  as_tsibble(index = yearmonth, key = c(holiday, county)) 
oslo %>% autoplot(guestnights)
oslo %>% gg_season(guestnights/1000)



data %>% 
  filter(#year(yearmonth)%in% 2008:2017,
    county %in%c("03 Oslo"),
    purpose != "total"
  ) %>% 
  mutate(gnpp = guestnights/population) %>% 
  as_tibble() %>% 
  select(yearmonth, purpose, gnpp) %>% 
  pivot_wider(names_from = purpose, values_from = gnpp)
  

oslo %>%
  filter(county == "03 Oslo")  %>% 
  gg_season(guestnights/population)
data %>%
  filter(purpose == "holiday")  %>% 
  gg_subseries(log(guestnights/population))
data %>% as_tibble() %>% 
  filter(purpose == "holiday",
         countycode %in% c(3,5,6,19)) %>% 
  mutate(mth =month(yearmonth)) %>% 
  group_by(countycode, mth) %>% 
  summarize(mnights = mean(1000*guestnights/population)) %>% 
 # filter(countycode == 12) %>% 
  ggplot(aes(x= mth, y = mnights, col = factor(countycode))) + geom_line()+
  geom_text(aes(label = countycode))


# INCLUDE 
# 03 OSLO
# 05 Oppland
# 06 Buskerud
# 08 Telemark
# 12 Hordaland
# 14 Sogn og Fjordane
# 16 Sør-Trøndelag
# 19 Troms
# 20 Finnmark

include <-c(3,5,6,19)# c(3,5,10,14,16,19)

holiday <- data %>% 
  filter(purpose == "holiday",
         countycode %in% include)  %>% 
  mutate(gpc = guestnights*1e-3)
holiday %>% autoplot(gpc)
holiday %>% gg_season(gpc)

holiday <- oslo %>% 
  mutate(gpc = 1000*guestnights/population)

models <- holiday %>% 
#  filter(county == "03 Oslo") %>% 
  model(#m1 = TSLM(guestnights/population ~ trend()+fourier(K = 1, period = "year")),
        #m2 = TSLM(guestnights/population ~ trend()+fourier(K = 2, period = "year")),
        #m3 = TSLM(guestnights/population ~ trend()+fourier(K = 3, period = "year")),
        #m4 = TSLM(guestnights/population ~ trend()+fourier(K = 4, period = "year")),
        #m5 = TSLM(guestnights/population ~ trend()+fourier(K = 5, period = "year")),
        #m6 = TSLM(guestnights/population ~ fourier(K = 6, period = "year")),
        #m1ut = TSLM(guestnights/population ~ fourier(K = 1, period = "year")),
        #m2ut = TSLM(guestnights/population ~ fourier(K = 2, period = "year")),
        #m3ut = TSLM(guestnights/population ~ fourier(K = 3, period = "year")),
        #m4ut = TSLM(guestnights/population ~ fourier(K = 4, period = "year")),
        #m5ut = TSLM(guestnights/population ~ fourier(K = 5, period = "year")),
        m6ut = TSLM(gpc ~fourier(K = 3, period = "year"))) 
models %>% glance() %>%
  group_by(county) %>%
  filter(AICc == min(AICc)) 

models %>% 
  forecast(holiday) %>% 
  filter(year(yearmonth) %in% 2010:2015) %>% 
  ggplot(aes(x=yearmonth, y = .mean, linetype = holiday))+
  geom_line() + geom_line(aes(y=1000*guestnights/population), col = "red")
  autoplot()+
  theme(legend.position ="none")+
  geom_line(aes(y= guestnights *1e-3))
coeffs <- models %>% coef()

season_coeffs <- coeffs %>% 
  select(holiday, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) 

names(season_coeffs)<- c(
  "county", "mu", #"trend",
  paste0(c("a","b"), rep(1:((ncol(season_coeffs)-2)/2), each = 2)
  ))
  
plot(season_coeffs[,-1])


A <- season_coeffs[,-c(1:2)] %>% as.matrix()
Sigma_season <- A %*% t(A) /2 
colnames(Sigma_season) <-rownames(Sigma_season) <- season_coeffs$county
det(Sigma_season)

empcov <- holiday %>%
  as_tibble() %>% 
  select(yearmonth,gpc, holiday ) %>% 
  pivot_wider(names_from = holiday, values_from = gpc) %>% 
  select(-yearmonth) %>% 
  cov()
empcor <- holiday %>% as_tibble() %>% 
  select(yearmonth,gpc, holiday ) %>% 
  pivot_wider(names_from = holiday, values_from = gpc) %>% 
  select(-yearmonth) %>% 
  cor()

Cor_season <- diag(1/sqrt(diag(Sigma_season)))%*%Sigma_season%*%diag(1/sqrt(diag(Sigma_season)))
colnames(Cor_season) <-rownames(Cor_season) <- season_coeffs$county
Cor_season %>% 
  as_tibble() %>% 
  mutate(county = rownames(Sigma_season)) %>% 
  pivot_longer(cols = 1:nrow(Sigma_season), 
               names_to = "county2", 
               values_to = "cor") %>% 
  mutate(source = "Season") %>% 
  bind_rows(
    empcor %>% 
      as_tibble() %>% 
      mutate(county = rownames(Sigma_season)) %>% 
      pivot_longer(cols =1:nrow(Sigma_season), 
                   names_to = "county2", 
                   values_to = "cor") %>% 
      mutate(source = "Empirical")

  ) %>% 
  mutate(county = factor(county),
         county2 = factor(county2, levels = rev(levels(county)))) %>% 
  ggplot(aes(x=county, y = county2, fill = cor)) + geom_tile() +
  scale_fill_viridis_c(direction = -1, name = "Correlation") + 
  scale_x_discrete(expand =c(0,0)) +
  scale_y_discrete(expand =c(0,0)) +
  guides(fill = guide_colorbar(barheight = 20))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = -60, hjust = 0))+
  facet_wrap(~source, ncol = 2)
Sigma_season %>% 
  as_tibble() %>% 
  mutate(county = rownames(Sigma_season)) %>% 
  pivot_longer(cols = 1:nrow(Sigma_season), 
               names_to = "county2", 
               values_to = "cor") %>% 
  mutate(source = "Season") %>% 
  bind_rows(
    empcov %>% 
      as_tibble() %>% 
      mutate(county = rownames(Sigma_season)) %>% 
      pivot_longer(cols = 1:nrow(Sigma_season), 
                   names_to = "county2", 
                   values_to = "cor") %>% 
      mutate(source = "Empirical")
    
  ) %>% 
  mutate(county = factor(county),
         county2 = factor(county2, levels = rev(levels(county)))) %>% 
  ggplot(aes(x=county, y = county2, fill = cor)) + geom_tile() +
  scale_fill_viridis_c(direction = -1, name = "Covariance") + 
  scale_x_discrete(expand =c(0,0)) +
  scale_y_discrete(expand =c(0,0)) +
  guides(fill = guide_colorbar(barheight = 20))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = -60, hjust = 0))+
  facet_wrap(~source, ncol = 2)



mu <- season_coeffs %>% pull(mu)
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
range(mu)
#targets <- seq(min(mu)+0.001, max(mu)-0.001, by =.01)
targets <- seq(250, max(mu)-0.001, by =.01)
#targets <- seq(26, 110,.1)
#Sigma_season <- matrix(covariances$cov, ncol = 2)
emp <- sapply(targets,
              MPT, cov = empcov, mu = mu)
season <- sapply(targets,
                 MPT, cov = Sigma_season, mu = mu)
seasonadjusted <- sapply(targets,
                         MPT, cov = empcov - Sigma_season, mu = mu)
ports_by_cov <- tibble(
  targets = targets, 
  "Empirical" = emp, 
  "Seasonal" = season,
  "Season-adjusted" = seasonadjusted
) %>% 
  pivot_longer(cols = 2:4) %>% 
  mutate(name = factor(name, levels = c("Empirical", "Seasonal", "Season-adjusted")))
tst <- ports_by_cov %>% 
  group_by(name) %>% 
  filter(value == min(value))

weights <- (MPTsummary <- rbind(
  "Empirical" = MPT(mu,empcov, target = targets[which.min(emp)], return_value = FALSE),
  "Seasonal" = MPT(mu,Sigma_season, target = targets[which.min(season)], return_value = FALSE),
  "Season-adjusted" = MPT(mu,empcov-Sigma_season, target = targets[which.min(seasonadjusted)], return_value = FALSE)) %>% 
    
    as.data.frame() %>% 
    rownames_to_column("covariance") %>% 
    as_tibble() %>% 
    mutate(covariance = factor(covariance, levels = c("Empirical", "Seasonal", "Season-adjusted")))) %>% 
  rename("buisness" = "weights1", 
         "pleasure" = "weights2") %>% 
  select(-c(2:3)) %>% 
  pivot_longer(cols = 2:3, names_to = "purpose", values_to = "weights") 

portfolios_ts <- holiday %>% as_tibble() %>% 
  rename("purpose" = "holiday")%>% 
  left_join(weights, by = "purpose", relationship = "many-to-many") %>% 
  group_by(yearmonth, covariance) %>% 
  summarize(portgpc = sum(weights*gpc)) %>%
  as_tsibble(key = covariance, index = yearmonth) %>% 
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
    #filter(year(yearmonth) %in% 2014:2019) %>% 
    ggplot(aes(x=yearmonth, y = portgpc)) + 
    geom_line()+
    facet_wrap(~covariance, ncol = 1) +
    geom_text(data= weights %>% filter(purpose == "buisness"),
              aes(x = as.Date("1998-01-01"), y = Inf, 
                  label = paste0("buisness weight: ", round(100*weights,1),"%    ",
                                 "SRS: ",round(sapply(weights,SRS, season= Sigma_season,sigma= empcov),3),
                                 "    SD: ", round(sapply(weights, SD, sigma = empcov),3))),
              vjust = 1.5, hjust=0)+
    #scale_x_date(expand = c(0,0), date_breaks = "1 year", date_labels = "%Y")+
    scale_x_yearmonth(date_breaks = "2 year", date_labels = "%Y", date_minor_breaks = "1 year")+
    scale_y_continuous(name = "Portfolio guestnight per capita")+
    theme(axis.title.x = element_blank(),
          panel.grid.major.x = element_blank())+
    geom_hline(data=MPTsummary, aes(yintercept = mean), lty = 2, col = 2))

acfplot <- portfolios_ts %>% 
  ACF(portgpc, lag_max = 36) %>%
  autoplot()  +
  scale_x_cf_lag(breaks = seq(0, 400, 3), limits = c(0,37), expand = c(0,0))+
  scale_y_continuous(name = "Autocorrelation")+
  theme(axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
ggpubr::ggarrange(timeplot+ theme(strip.text =element_blank()),
                  acfplot+ theme(axis.title.x = element_blank()))
ggsave(file = "figures/case3_Time_ACF_combined.png", width = 12, height = 7)

ports_by_cov %>% 
  ggplot(aes(x = (value), y = targets, color = name, group = name)) + 
  geom_path() + 
  scale_x_continuous("Variance measure", limits = c(0, 1210), expand = c(0,0))+
  scale_y_continuous("Hotel guestnights per thousand capita", breaks = seq(0,600,5))+
  geom_segment(data = ports_by_cov %>% group_by(name) %>% 
                 filter(targets == 125 & name == "Seasonal"|targets == 80 & name == "Empirical"|
                          targets== 85 & name == "Season-adjusted"),
               aes(x = -Inf, xend = value, y = targets, yend= targets, color = name),
               arrow = arrow(length = unit(0.02, "npc")))+
  geom_hline(data = tst, aes(yintercept = targets, color = name, x= NULL), lty = 2)+
  theme(legend.title = element_blank(), 
        legend.position = c(.12,.88),
        legend.text  = element_text(size = 12),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ggtitle("Oslo hotel guestnights")
ggsave("figures/case3_oslo_buisness_pleasure_efficient_frontiers.png", width = 8, height = 4)
  
season <- MPT(mu = mu, cov = Sigma_season, target = 270, return_value = FALSE)
weights <- (MPTsummary <- rbind(
  "Empirical" = MPT(mu,empcov, target = 270, return_value = FALSE),
  "Seasonal" = MPT(mu,Sigma_season, target = 270, return_value = FALSE),
  "Season-adjusted" = MPT(mu,empcov-Sigma_season, target = 270, return_value = FALSE)) %>% 
    
    as.data.frame() %>% 
    rownames_to_column("covariance") %>% 
    as_tibble() %>% 
    mutate(covariance = factor(covariance, levels = c("Empirical", "Seasonal", "Season-adjusted")))) %>% 
  select(-c(2:3)) %>% 
  pivot_longer(cols = 2:5, names_to = "locID", values_to = "weights")# %>% 
 # mutate(locID = factor(locID, levels = c("Solar PV" ,"WP portfolio")))





rm(tourism)
data("tourism", package = "fpp3")
holidays <- tourism |>
  filter(Purpose == "Holiday") |>
  group_by(State) |>
  summarise(Trips = sum(Trips))

holidays <- holidays %>% 
  group_by(State) %>% 
  mutate(Trips = Trips/mean(Trips))
holidays %>% autoplot((Trips))
empcov <- holidays %>% filter(State %in% c("ACT","Northern Territory", "Queensland")) %>% 
  #group_by(State) %>% 
 # mutate(Trips = Trips-mean(Trips)) %>%  
  pivot_wider(names_from = State, values_from = Trips) %>% as_tibble() %>% select(-Quarter) %>% cov()
empcorlong <- empcov %>%
  as_tibble() %>%
  mutate(loc1 = colnames(empcov)) 

mu <- holidays %>% filter(State %in% c("ACT","Northern Territory", "Queensland")) %>% 
  as_tibble() %>% 
  group_by(State) %>% summarize(m = mean(Trips)) %>% pull(m)
  

#names(empcorlong)[-9] <- colnames(empcov)
empcorlong <- empcorlong %>%   pivot_longer(cols = 1:3, names_to = "loc2")
empcorlong %>%
  filter((loc1 %in%c("ACT","Northern Territory", "Queensland")),
         (loc2 %in%c("ACT","Northern Territory", "Queensland"))) %>% 
  ggplot(aes(x=loc1,y=loc2, fill = value))+geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
holidays %>% 
  model(m2=TSLM((Trips) ~fourier(K = 2, period = "year"))) %>% 
  fitted() %>% 
  filter(State %in% c("ACT","Northern Territory", "Queensland")) %>% 
  #group_by(State) %>%
  #mutate(m = mean(.fitted)) %>% 
  filter(year(Quarter) %in% 2005:2006) %>% 
  autoplot(.fitted)
A <- 
  (coeffs <- (model <-holidays %>%
  filter(State %in% c("ACT","Northern Territory", "Queensland")) %>% 
 #   group_by(State) %>% 
  #  mutate(Trips = Trips-mean(Trips)) %>% 
  model(m2=TSLM((Trips) ~0+fourier(K = 2, period = "year")))) %>% 
  coef() %>% 
  select(State, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate)) %>% 
  column_to_rownames("State") %>%
  as.matrix()
colnames(A) <- c("a1","b1", "a2")
names(coeffs) <- c("State","a1","b1", "a2")
#coeffs %>% mutate(r = sqrt(a1^2 +b1^2),
#                  phi = ifelse(a1>0, atan(b1/a1),ifelse(b1>=0, atan(b1/a1)+pi, atan(b1/a1)-pi)))
det(Sigma_season <- A%*%t(A)/2)
det(empcov)
MPT(mu = rep(0,3), cov = Sigma_season, target = 0)
Sigma_season
empcov
(Sigma_season_cor<-(diag(1/sqrt(diag(Sigma_season)))%*%Sigma_season%*%diag(1/sqrt(diag(Sigma_season)))))


MPT(mu = mu, cov = empcov, target = 1, return_value = FALSE)
MPT(mu = mu, cov = Sigma_season, target = 500, return_value = FALSE)



lines(sss <- sapply(targets, MPT, mu = A[,1], cov = Sigma_season), targets)
targets[which.min(sss)]
targets[which.min(ssemp)]



aus_arrivals %>% 
  mutate(diff = c(NA,diff(Arrivals, 1))) %>% 
  #select(-.model) %>%
  filter(between(Quarter, yearquarter("1982 Q1"), yearquarter("2011 Q4"))) %>% 
  #gg_season(diff)
  model(TSLM(diff~fourier(K = 2, period = "year"))) %>% 
  coef()
  residuals() %>% 
  #gg_season(.resid)
  autoplot(.resid) 

  