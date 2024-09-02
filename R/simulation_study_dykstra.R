# Simulation study
library(tidyverse)
library(quadprog)
library(fpp3)
library(ggpubr)
theme_set(theme_bw()+
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  strip.background = element_rect(fill = "transparent",  color = "transparent")))
restab <- tibble()

SRS <- function(w, season, sigma){
  w <- c(w,1-w)
  return(as.numeric(sqrt(round(t(w)%*%season%*%w / t(w)%*%sigma%*%w, digits = 8))))
}


phis <- rep(c(0.05, 0.5, 1),2)
cors <- c(rep(0,3), rep(-.8,3))
for(i in 1:length(phis)){
  # What we adjust: 
  phi1 <- phis[i]*pi 
  #phi1 <- .01*pi 
  #phi1 <- .5*pi 
  r1 <- 3
  rescor <- cors[i] # residual correlation
  #ressd <- 2 # residual sd
  ressd <- 3 # residual sd
  
  # -------------------
  # Fixed parameters: 
  # -------------------
  r2 <- 1
  sample_size <- 365*15                                       
  sample_meanvector <- c(0,0)
  mu <- c(100, 100)
  sample_covariance_matrix <- matrix(c(ressd^2,
                                       ressd^2*rescor,
                                       ressd^2*rescor, 
                                       ressd^2),
                                     ncol = 2)
  Amat <- cbind(1, mu, diag(length(mu)))
  dat <- tibble(
    t = 1:sample_size,
    date = as.Date("2000-12-31")+t)
  set.seed(2134)
  sample_distribution <- MASS::mvrnorm(n = sample_size,
                                       mu = sample_meanvector, 
                                       Sigma = sample_covariance_matrix)
  dat2 <- dat %>% mutate(
    x = mu[1]+r1*cos(2*pi*t/365-phi1)+sample_distribution[,1],
    y = mu[2]+r2*cos(2*pi*t/365)+sample_distribution[,2]
  ) 
  timeplotExample <-  dat2 %>% 
    rename("X1"="x","X2"="y") %>% 
    pivot_longer(cols = c("X1","X2")) %>% 
    filter(year(date)%in% 2013:2015) %>% 
    ggplot(aes(x=date, y = value, col = name)) + geom_line()+
    #labs(title = "Time plot")+
    scale_x_date(date_breaks = "1 years",
                 date_labels = "%Y",
                 expand = c(0,0), limits = c(as.Date("2012-11-30"),as.Date("2015-12-31")))+
    scale_y_continuous(name = "", breaks = seq(0,150,2))+
    geom_hline(yintercept = mu[1], lty = 2)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          legend.position = c(0.95,0.93),
          legend.background = element_rect(fill = "transparent"))+
    geom_segment(aes(x=as.Date("2012-12-15"),xend = as.Date("2012-12-15"),
                     y = 95, yend = 105), color = "gray", lwd = 5)
  if(i == 3)
    ggsave(timeplotExample,file = "figures/simulation_example.png", width = 8, height = 4)
  dat2 %>% 
    ggplot(aes(x=x, y = y)) + geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Scatterplot")
  
  sol <- Dykstra::dykstra(Dmat = cov(dat2[,c("x","y")]), 
                  dvec = rep(0, 2),
                  Amat = Amat[,-2],
                  bvec = c(1,#mu[1],
                           rep(0,length(mu))),#rangemu[which.min(js)]), 
                  meq  = 1)
  sol.clim <- Dykstra::dykstra(Dmat = matrix(
    c(r1^2, r1*r2*cos(phi1),
      r1*r2*cos(phi1),r2^2)/2,
    ncol=2, nrow=2), 
    dvec = rep(0, 2),
    Amat = Amat[,-2],
    bvec = c(1,#mu[1],
             rep(0,length(mu))),#rangemu[which.min(js)]), 
    meq  = 1)
  sol.noclim <- Dykstra::dykstra(Dmat = cov(dat2[,c("x","y")])-matrix(
    c(r1^2, r1*r2*cos(phi1),
      r1*r2*cos(phi1),r2^2)/2,
    ncol=2, nrow=2), 
    dvec = rep(0, 2),
    Amat = Amat[,-2],
    bvec = c(1,#mu[1],
             rep(0,length(mu))),#rangemu[which.min(js)]), 
    meq  = 1)
  weights <- rbind("Naive" = sol$solution,
                   "Seasonal" = sol.clim$solution,
                   "Seasonal-adjusted"= sol.noclim$solution)
  colnames(weights) <- c("X", "Y")
  weights
  
  weights <-as_tibble(weights) %>% mutate(SRS=sapply(X, SRS, season = matrix(
    c(r1^2, r1*r2*cos(phi1),
      r1*r2*cos(phi1),r2^2)/2,
    ncol=2, nrow=2), sigma = cov(dat2[,c("x","y")])))
  dat3 <- dat2 %>% 
    mutate("Naive" = weights$X[1]*x + weights$Y[1]*y,
           "Seasonal" = weights$X[2]*x + weights$Y[2]*y,
           "Seasonal-adjusted" = weights$X[3]*x + weights$Y[3]*y) %>% 
    pivot_longer(cols = 5:7, values_to = "power", names_to = "source") %>% 
    mutate(source = factor(source, 
                           levels = c("Naive", "Seasonal", "Seasonal-adjusted"))) 
  
  if(i == 3){
    fitted <- dat3 %>% 
      as_tsibble(key = source, index = date) %>% 
      model(LM = TSLM(power~fourier(K = 1, period = "year"))) %>% 
      fitted()
    
    
    (  timeplot <- dat3 %>% 
        filter(lubridate::year(date)%in%2013:2015) %>% 
        ggplot(aes(x=date, y= power,group = source)) + 
        geom_line(col = "blue")+
        scale_y_continuous("Portfolio output")+
        scale_x_date(date_breaks = "1 years",
                     date_labels = "%Y",
                     expand = c(0,0), limits = c(as.Date("2012-11-01"),as.Date("2015-12-31")))+
        theme(axis.title.x = element_blank(),
              legend.title = element_blank())+
        facet_wrap( ~ source, nrow = 3, strip.position = "right")+
        geom_line(data= fitted %>% filter(year(date)%in%2013:2015),
                  aes(x = date, y = .fitted), lty = 2, lwd = .9, color = "orange")+
        geom_segment(aes(x=as.Date("2012-12-01"),xend = as.Date("2012-12-01"),
                         y = 95, yend = 105), color = "gray", lwd = 5))
    ggsave(timeplot, file = "figures/simulation_timeplots.png", width = 8, height = 5)
    
    (  acfplot <- dat3 %>% 
        as_tsibble(index = date, key = source) %>% 
        ACF(power, lag_max = 400) %>% autoplot() +
        scale_x_cf_lag( breaks = seq(0,400,30), limits = c(0,401), expand = c(0,0))+
        ylab("Autocorrelation"))
    
    ggsave(acfplot, file = "figures/simulation_ACF.png", width = 8, height = 5)
    
    ggarrange(timeplot+ theme(strip.text =element_blank()), acfplot+ theme(axis.title.x = element_blank()))
    ggsave(file = "figures/simulation_Time_ACF_combined.png", width = 10, height = 8)
    ggarrange(timeplotExample, 
              ggarrange(timeplot+ theme(strip.text =element_blank()), acfplot+ theme(axis.title.x = element_blank()),
                        labels = c("B","C")),
              labels = c("A",NA),
              nrow = 2)
    ggsave(file = "figures/simulation_figure_combined.pdf", width = 10, height = 10)
  }
  restab <- restab %>% 
    bind_rows(
      dat3 %>% 
        as_tsibble(index = date, key = source) %>% 
        ACF(power, lag_max = 400) %>% as_tibble() %>% 
        mutate(lag = as.numeric(lag)) %>% 
        group_by(source) %>% 
        summarize("lag365" = acf[lag==365],
                  "SSlag1-400" = sum(acf^2),
                  "SSlag355-75" = sum(acf[lag %in% 355:375])) %>% 
        mutate(phi1 = phi1/pi, r1 = r1, rescor = rescor, ressd = ressd,
               weightX = weights$X,
               SRS= weights$SRS) %>% 
        right_join(
          dat3 %>% group_by(source) %>% 
            summarize(sd = sd(power))
        ))
}
library(xtable)
xtable(
  restab %>% relocate(source, phi1, rescor, weightX, sd, SRS) %>% arrange(rescor, source, phi1) %>%
    select(-r1,-ressd) %>% arrange(source, desc(rescor)) %>%
    filter(!(source == "Seasonal-adjusted" & phi1 %in% c(0.05,.5)))
  #filter(source != "Climate adjusted")
) %>% 
  print(include.rownames = FALSE)#, file = "output/simulationtable.tex")


