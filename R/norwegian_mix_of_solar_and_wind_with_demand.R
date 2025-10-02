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
  SD <- function(w, sigma){
    w <- c(w,1-w)
    return((t(w)%*%sigma%*%w))
  }
  SRS <- function(w, season, sigmaZ){
    w <- c(w,1-w)
    return((t(w)%*%season%*%w / t(w)%*%(season+sigmaZ)%*%w))
  }
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
    # From Table A1 of HC8lleland et al (2023)
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
    filter(locID %in% c( "WP portfolio","Solar PV"))
  
  
  # Demand from statnett: https://driftsdata.statnett.no/Web/Download/
  consumption <- lapply(list.files(path = "data/statnett", full.names =TRUE),
                 read_csv2) %>% 
    bind_rows() %>% 
    mutate(datetime = as.POSIXct(`Time(Local)`,format = "%d.%m.%Y %H:%M:%S")) %>% 
    select(datetime, Consumption)
  
  
  daily_consumption <- consumption %>% 
    mutate(date=date(datetime)) %>% 
    group_by(date) %>% 
    summarize(Consumption=sum(Consumption,na.rm=T))
  
  size_of_solar_park = 10000
  
  
  power <- power %>% 
    pivot_wider(names_from = "locID", values_from = "CF") %>% 
    mutate(`Solar PV` = `Solar PV`*2000*24*size_of_solar_park/1e6,
           # We assume 7500 solar panels per unit solar and 2000W per panel
           `WP portfolio`=`WP portfolio`*15e6*24/1e6) %>% 
    # Units are MWh (per day)
    left_join(daily_consumption) %>% 
    filter(year(date)>2005) %>% as_tibble() %>% 
    filter(Consumption > 2e5) %>% 
    na.omit()
  
  
  power %>%  ggplot(aes(x=date, y = Consumption)) + geom_line()+
    geom_smooth(method = "lm")
  
  # Detreding consumption
  detrend <- power %>% 
    mutate(t = 1:n())
  con.lm <- lm(Consumption ~ t, data = detrend)
  power$Consumption <- power$Consumption + coef(con.lm)[2]*nrow(power)- coef(con.lm)[2]*(1:nrow(power))
  
  # Fixed at 2020 level.
  
  # Estimate matrices :
  Amat <- power %>% as_tibble() %>% 
    forecast::msts(seasonal.periods = 365.25) %>%
    forecast::fourier(K = 4)
  Y <- power %>% select(-date) %>% as.matrix()
  
  # Selecting K = 4
  # Use K = 5 above and run the following models: 
  AIC(lm(Y[,1]~1+Amat),
      lm(Y[,1]~1+Amat[,-(ncol(Amat)-1:0)]),
      lm(Y[,1]~1+Amat[,-(ncol(Amat)-3:0)]),
      lm(Y[,1]~1+Amat[,-(ncol(Amat)-5:0)]))
  AIC(lm(Y[,2]~1+Amat),
      lm(Y[,2]~1+Amat[,-(ncol(Amat)-1:0)]),
      lm(Y[,2]~1+Amat[,-(ncol(Amat)-3:0)]),
      lm(Y[,2]~1+Amat[,-(ncol(Amat)-5:0)]))
  AIC(lm(Y[,3]~1+Amat),
      lm(Y[,3]~1+Amat[,-(ncol(Amat)-1:0)]),
      lm(Y[,3]~1+Amat[,-(ncol(Amat)-3:0)]),
      lm(Y[,3]~1+Amat[,-(ncol(Amat)-5:0)]))
  
  
  ffit <- lm(Y~ 1 + Amat)
  A <- t(coef(ffit)[-1, ])
  
  A <- A[,c(seq(2,ncol(A),2),seq(1,ncol(A)-1,2))]
  
  SigS <- nrow(power)/(nrow(power)-1) * A%*%t(A)/2
  SigZ <- cov(residuals(ffit))
  Sig  <- cov(Y)
  
  cov2cor(SigS)
  cov2cor(SigZ)
  cov2cor(Sig)
  # store fitted values for purpose of plotting:
  
  
  plot(A[3,])
  .fitted <- fitted(ffit)
  fitted.df <- left_join(reshape2::melt(Y) %>% rename("CF"=value), 
                         reshape2::melt(.fitted) %>% 
                           rename(fitted = value)) %>% 
    mutate(date = power$date[Var1]) %>% 
    rename(locID= Var2) %>% 
    select(-Var1) %>% 
    as_tibble()
  
  
  
  
  fig1 <- fitted.df %>% 
    mutate(isConsumption = ifelse(locID=="Consumption",
                                  "Electricity consumption (GWh)",
                                  "Capacity factor")) %>% 
    mutate(locID = ifelse(locID == "WP portfolio", "Offshore Wind", as.character(locID))) %>% 
    mutate(CF = ifelse(locID == "Solar PV",
                       CF *1e6/(2000*24*size_of_solar_park), 
                       ifelse(locID == "Offshore Wind",
                              CF /(15*24), 
                              CF/1e3))
                       ,
           fitted  = ifelse(locID == "Solar PV",
                            fitted  *1e6/(2000*24*size_of_solar_park), 
                       ifelse(locID == "Offshore Wind",
                              fitted  /(15*24), 
                              fitted/1e3 ))
    ) %>% 
    mutate(locID = factor(locID,
                          levels = c(
                            "Offshore Wind",
                            "Solar PV",
                            "Consumption"
                          ))) %>% 
    
   # dplyr::select(-locID) %>% 
    #rename(locID=locID2) %>% 
    ggplot(aes(x =date, y = fitted, color = locID))+
    geom_line(lwd = .7)+
    geom_point(aes(x= date, y = CF, color = locID), size = .2, show.legend = FALSE)+
      facet_wrap(~isConsumption,
                 scales= "free_y", 
                 ncol = 1,
                 strip.position = "left") +
    scale_y_continuous(expand = c(0,0))+
    scale_x_date(expand = c(0,0), limits = c(as.Date("2005-12-16"),as.Date("2019-12-31")),
                 date_breaks = "1 year",
                 date_labels = "%Y",
                 name = "") +
    #scale_color_manual(values = c( "blue","red", "green"))+
    guides(color = guide_legend(override.aes = list(size = 2, lwd =5)))+
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          legend.margin=margin(-20, 0, 0, 0),
          legend.background = element_rect(fill = "transparent", color = "transparent"),
          strip.placement = "outside",
          axis.title.y = element_blank())
  fig1
  ggsave(fig1, file = "figures/case2_wind_solar_time_plot_season.pdf", width = 8, height = 5)
  
  
  
  
  
  
  mu <- coef(ffit)[1,]
  
  MPT <- function(mu, cov, p_target, return_value = TRUE, return_weight=FALSE){
    
    Amat <- cbind(c(0,0,1), mu, c(1,0,0),c(0,1,0))
    bvec <- c(-1, (p_target-1)*mu[3], rep(0,length(mu)-1))
    sol <- solve.QP(Dmat = cov, 
              dvec = rep(0, length(mu)),
              Amat= Amat,
              bvec = bvec, 
              meq = 2)
    if(return_value) return(sol$value)
    if(return_weight) return(sol$solution)
    return(c("p_target"=p_target, "value" = sol$value, "weights" = sol$solution))
  }
  p_targets <- seq(0.01, 1, by = 0.01)
  emp <- sapply(p_targets,
                MPT, cov = Sig, mu = mu)
  season <- sapply(p_targets,
                MPT, cov = SigS, mu = mu)
  seasonadjusted <- sapply(p_targets,
                   MPT, cov = SigZ, mu = mu)
  ports_by_cov <- tibble(
   p_targets = p_targets, 
    "Balanced" = emp, 
    "Seasonal" = season,
    "Season-adj." = seasonadjusted
  ) %>% 
    pivot_longer(cols = 2:4) %>% 
    mutate(name = factor(name, levels = c( "Season-adj.", "Balanced", "Seasonal")))
    
  # --------------
  p_targets <- seq(0.01, 10, by = 0.01)
  emp <- sapply(p_targets,
                MPT, cov = Sig, mu = mu, return_weight=T, return_value = F) %>% t() %>% 
    bind_cols(p_targets) %>% 
    mutate(stragegy = "Balanced") 
  season <- sapply(p_targets,
                   MPT, cov = SigS, mu = mu, return_weight=T, return_value = F) %>% t()%>% 
    bind_cols(p_targets) %>% 
    mutate(stragegy = "Seasonal")
  seasonadjusted <- sapply(p_targets,
                           MPT, cov = SigZ, mu = mu, return_weight=T, return_value = F) %>% t()%>% 
    bind_cols(p_targets) %>% 
    mutate(stragegy = "Season-adj.")
  ports_by_cov <-  bind_rows(emp, season, seasonadjusted)
  names(ports_by_cov) <- c("Solar PV", "WP portfolio", "Consumption", "p_target", "strategy")
  ports_by_cov <- ports_by_cov %>% 
    mutate(mix = `WP portfolio`/(`WP portfolio`+`Solar PV`),
           strategy = factor(strategy, levels = c("Season-adj.", "Balanced", "Seasonal"))) 
  thresholds <- ports_by_cov %>% group_by(strategy) %>% 
    filter(round(mix,5) <1) %>% 
    filter(p_target == min(p_target))
  
  more_than_30gw <- ports_by_cov %>% 
    group_by(strategy) %>% 
    mutate(facet=factor("Solar PV proportion",
                        levels = c("Solar PV proportion","SRS"))) %>% 
    filter(`WP portfolio`>=2000) %>% 
    filter(p_target == min(p_target))
  
  
  W <- ports_by_cov[,1:3] %>% as.matrix()
  SRS <- function(w, SigS, SigZ){
    return((t(w)%*%SigS%*%w / t(w)%*%(SigS+SigZ)%*%w))
  }
  
  ports_by_cov$SRS = apply(W, 1, SRS, 
                           SigS=SigS, SigZ=SigZ)
  
  
  
  ports_by_cov %>% filter(p_target <=3.06) %>% 
    mutate("Solar PV proportion" = 1-mix) %>% 
    pivot_longer(cols = c(`Solar PV proportion`, "SRS"), names_to = "facet", values_to = "value") %>% 
    mutate(facet = factor(facet, 
                          levels = c("Solar PV proportion", "SRS"))) %>% 
    filter(facet == "Solar PV proportion") %>% 
    ggplot(aes(x=p_target, y = value, color = strategy)) + geom_line()+
    scale_x_continuous(expand = c(0,0),
                       breaks = c(seq(0,5,.5), more_than_30gw$p_target[1:3]),
                       labels = scales::percent,
                       name = expression(p[target]))+
    scale_y_continuous(name = "Solar PV proportion",
                       labels = scales::percent,
                       breaks = c(seq(0,1,.25), 1-more_than_30gw$mix[c(1,3)]),
                       limits = c(-0.01,1.0), expand = c(0,0))+
    #facet_wrap( ~facet, ncol = 1, strip.position = "left")+
    geom_segment(aes(x = p_target,xend = p_target, y = -Inf, yend = 1-mix, color = strategy) , data = more_than_30gw, lty = 2, show.legend = FALSE)+
    geom_segment(aes(x = -Inf,xend = p_target, y = 1-mix, yend = 1-mix, color = strategy) , data = more_than_30gw, lty = 2, show.legend = FALSE)+
    # geom_vline(data = more_than_30gw %>% 
    #              mutate(facet = factor("SRS", levels = c("Solar PV proportion", "SRS"))),
    #            aes(xintercept = p_target, color = strategy), lty = 2, show.legend = FALSE)+
    # geom_point(aes(x = p_target, y = mix), 
    #            data = more_than_30gw, show.legend = FALSE, color = "black",
    #            shape = "x", size = 3)+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = "transparent", 
                                           color = "transparent"),
          plot.margin = margin(5, 5, 5, 5),
          legend.spacing.y = unit(0, "cm"),
          legend.box.spacing = unit(0, "cm"),
          legend.key = element_rect(fill = "transparent"),
          strip.placement = "outside")
  ggsave("figures/demand_energymix_per_target.pdf", width = 8, height = 3)
  
  
  
  power_long <- power %>% 
      pivot_longer(cols = 2:4, names_to ="asset", values_to = "power") 
  
  time_plots <- ports_by_cov %>% 
    filter(p_target == 1) %>% 
    pivot_longer(cols = 1:3, names_to ="asset", values_to = "weight") %>% 
    right_join (power_long, relationship = "many-to-many") %>% 
    mutate(Production = power*weight) %>% 
    filter(asset != "Consumption") %>% 
    group_by(date, strategy) %>% 
    summarize(Production = sum(Production)) %>% 
    left_join(daily_consumption) %>% 
    mutate("Net-balance" = Production-Consumption) %>% 
    pivot_longer(cols = 3:5) %>% 
    #filter(year(date) == 2019) %>% 
    mutate(yday = date-years(year(date))) %>% 
    group_by(yday, strategy, name) %>% 
    summarize(value = mean(value, na.rm=T)) %>% 
    mutate(name = paste("Average",name)) %>% 
    mutate(name = factor(name, 
                         levels = paste("Average",
                                        c("Production",
                                          "Consumption",
                                          "Net-balance"))))
  
  time_plots %>% 
    ggplot(aes(x=yday, y = value/1e3, color = name)) +
    geom_line()+ 
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey", 
               data = tibble(name = "Average Net-balance"))+
    facet_grid(name == "Average Net-balance"~strategy, scales = "free_y")+
    scale_x_date(date_labels = "%b", expand = c(0,0))+
    scale_y_continuous("Power (GWh)")+
    #scale_color_manual(values = c("red","blue","darkgreen"))
    theme(strip.text.y = element_blank(),
          legend.title = element_blank(),
          legend.position =  "top",
          axis.title.x = element_blank())
  ggsave(file = "figures/demand_portfolios_and_netbalance_averaged.pdf", width = 8, height = 5)
  
  ports_by_cov %>% 
    filter(p_target == 1) %>% 
    mutate(mix = 1-mix)
  ports_by_cov %>% 
    filter(p_target == 1) %>% 
    pivot_longer(cols = 1:3, names_to ="asset", values_to = "weight") %>% 
    right_join (power_long, relationship = "many-to-many") %>% 
    mutate(Production = power*weight) %>% 
    filter(asset != "Consumption") %>% 
    group_by(date, strategy) %>% 
    summarize(Production = sum(Production)) %>% 
    left_join(daily_consumption) %>% 
    mutate("Net-balance" = Production-Consumption) %>% 
    pivot_longer(cols = 3:5) %>% 
    filter(name == "Net-balance") %>% 
    group_by(strategy) %>% 
    summarize(sqrt(mean(value^2)))
  
  W <- ports_by_cov[,1:3] %>% as.matrix()
  SRS <- function(w, SigS, SigZ){
    return((t(w)%*%SigS%*%w / t(w)%*%(SigS+SigZ)%*%w))
  }
  
  ports_by_cov$SRS = apply(W, 1, SRS, 
        SigS=SigS, SigZ=SigZ)
  
  
  ports_by_cov %>% 
    #filter(p_target <3.01) %>% 
    ggplot(aes(x=p_target, y = SRS, color = strategy)) + 
    geom_line()+
    scale_x_continuous(expand = c(0,0),
                       breaks = c(seq(0,5,.5)),
                       labels = scales::percent,
                       name = expression(p[target]))
    
  # Plot only up to Nyquist frequency
  half_n <- floor(n / 500)
  plot(freq[1:half_n], power[1:half_n], type = "l",
       xlab = "Frequency", ylab = "Power", main = "Power Spectrum")
  
  
  
  # --------- POWER SPECTRUM STUFF -----------
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
  
  X <- Y[,1:3]
  pow<- matrix(NA_real_, ncol = 3, nrow=10)
  for(k in 1:(nrow(pow))){
    pow[k,] <- Mod(colMeans(X*exp(-complex(1,0,1)*2*pi/365.25*k*(1:length(X)))))^2
  }
  colnames(pow) <- colnames(Y)
  colnames(pow)[2] <- "Offshore wind"
  tibble(k=1:nrow(pow)) %>% 
    bind_cols(pow) %>% 
    pivot_longer(cols = 2:4) %>% 
    group_by(name) %>% 
    mutate(cumsum = cumsum(value)/sum(value)) %>% 
    mutate(name = factor(name, levels = c( "Offshore wind","Solar PV","Consumption"))) %>% 
    ggplot(aes(x=k, y=cumsum, color = name)) +
    #geom_line()+
    geom_point() + geom_line(show.legend=FALSE)+
    #geom_hline(yintercept=0.99, color = "steelblue", linewidth = .7,lty=1)+
    #geom_vline(xintercept = 6, lty = 2, color = "steelblue")+
    scale_x_continuous(breaks =seq(1,10,1),
                       limits = c(.5,10.5),
                       expand = c(0,0))+
    scale_y_continuous(labels = scales::percent, 
                       breaks = c(seq(.9,1,.02)),
                       name = expression(P[k](omega)),
                       limits =c(.9, 1),
                       expand =c(0,0))+
    #scale_color_manual(values = c("red","blue","darkgreen"))+
    theme(legend.title = element_blank(),
          legend.position.inside = T,
          legend.position = c(.85,.25),
          legend.background = element_rect(fill = "transparent",
                                           color = "transparent"))
ggsave("figures/power_spectrum_of_K_energymix_demand.pdf", width = 6, height = 2.5)  
