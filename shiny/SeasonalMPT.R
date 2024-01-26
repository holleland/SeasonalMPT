# Shiny application for Seasonal MPT simulation
# Author: Sondre Hølleland

library(shiny)
library(tidyverse)
library(quadprog)
theme_set(theme_bw()+
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  strip.background = element_rect(fill = "white", color = "transparent")))
sample_size <- 365*15                                       
sample_meanvector <- c(0,0)
set.seed(2134)
sample_covariance_matrix <- matrix(c(5,2.5,2.5, 5),
                                   ncol = 2)

# create bivariate normal distribution
sample_distribution <- MASS::mvrnorm(n = sample_size,
                                     mu = sample_meanvector, 
                                     Sigma = sample_covariance_matrix)

MPT <- function(mu, cov, target, return_value = TRUE){
  Amat <- cbind(1, mu, diag(length(mu)))
  bvec <- c(1, target, rep(0,length(mu)))
  sol <- solve.QP(Dmat = cov, 
                  dvec = rep(0, length(mu)),
                  Amat= Amat,
                  bvec = bvec, 
                  meq = 2)
  if(return_value) return(sol$value)
  return(c("mean" = target, "value" = sol$value, "weights" = sol$solution))
}
SD <- function(w, sigma){
  w <- c(w,1-w)
  return(sqrt(t(w)%*%sigma%*%w))
}
SRS <- function(w, season, sigma){
  w <- c(w,1-w)
  return(sqrt(t(w)%*%season%*%w / t(w)%*%sigma%*%w))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Seasonal MPT simulation"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("r1",
                  "Amplitude X:",
                  min = 0.1, 
                  max = 6, 
                  value = 3,
                  step = 0.1),
      sliderInput("phaseshift",
                  "Phase shift (as fraction of pi):",
                  min = 0,
                  max = .99,
                  step = 0.01,
                  value = .99),
      sliderInput("rescor",
                  "Residual correlation",
                  min = -.99,
                  max = .99,
                  step = 0.01,
                  value = 0),
      sliderInput("ressd",
                  "Residual standard deviation:",
                  min = 0.1,
                  max = 3,
                  step = 0.1,
                  value = 3),
      textOutput("dailyphase"),
      tableOutput("empiricalcorrelation"),
      tableOutput("portfoliosSummary")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("timePlot"), 
      plotOutput("portfolioPlot")
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  dat <- tibble(
    t = 1:sample_size,
    date = as.Date("2009-12-31")+t)
  
  output$empiricalcorrelation <- renderTable({
    set.seed(2134)
    sample_covariance_matrix <- matrix(c(input$ressd^2,
                                         input$ressd^2*input$rescor,
                                         input$ressd^2*input$rescor, 
                                         input$ressd^2),
                                       ncol = 2)
    sample_distribution <- MASS::mvrnorm(n = sample_size,
                                         mu = sample_meanvector, 
                                         Sigma = sample_covariance_matrix)
    dat2 <- dat %>% mutate(y = cos(2*pi*(t/365-input$phaseshift/2))+sample_distribution[,2],
                           x = input$r1*cos(2*pi*t/365)+sample_distribution[,1]) 
    tab <- cov(dat2[,c("x","y")])
    colnames(tab)<-c("Ex","Ey")
    season <- matrix(c(input$r1^2, input$r1*cos(pi*input$phaseshift),input$r1*cos(pi*input$phaseshift),1^2),
                  ncol=2, nrow=2)/2
    rownames(season)<-colnames(season)<-c("Sx","Sy")
    tab <- tab %>% bind_cols(season) %>% as.matrix()
    rownames(tab) <- c("x","y")
    tab
  },rownames = TRUE, caption = "Covariance matrix (E = Emprical, S = Seasonal)")
  output$dailyphase <- renderText({
    paste0("Phase shift in days: \t", round(input$phaseshift*365/2,0))
  })
  output$portfoliosSummary <- renderTable({
    set.seed(2134)
    sample_covariance_matrix <- matrix(c(input$ressd^2,
                                         input$ressd^2*input$rescor,
                                         input$ressd^2*input$rescor, 
                                         input$ressd^2),
                                       ncol = 2)
    seasonal <- matrix(c(input$r1^2, input$r1*cos(pi*input$phaseshift),input$r1*cos(pi*input$phaseshift),1^2),
                  ncol=2, nrow=2)/2
    sample_distribution <- MASS::mvrnorm(n = sample_size,
                                         mu = sample_meanvector, 
                                         Sigma = sample_covariance_matrix)
    dat2 <- dat %>% mutate(y = cos(pi*(t/365-input$phaseshift/2))+sample_distribution[,2],
                           x = input$r1*cos(pi*t/365)+sample_distribution[,1]) 
    COV <- cov(dat2[,c("x","y")])
    res <- rbind(
      "Empirical"= MPT(mu=c(1,1), cov = COV, target = 1, return_value = FALSE),
    "Seasonal" = MPT(c(1,1), cov = seasonal, target = 1, return_value = FALSE),
    "Season-adj"=MPT(c(1,1), cov = sample_covariance_matrix, target = 1, return_value = FALSE))[,3]
    res <- cbind(
      "Weight X"= res,
      "SRS" =  as.numeric(sapply(res,SRS, season = seasonal, sigma = COV)),
      "SD" = c(SD(res[1], COV),
               SD(res[2], COV),
               SD(res[3], COV))
    )
    res
  },rownames = TRUE, caption = "Portfolio summary")
  output$timePlot <- renderPlot({
    set.seed(2134)
    sample_covariance_matrix <- matrix(c(input$ressd^2,
                                         input$ressd^2*input$rescor,
                                         input$ressd^2*input$rescor, 
                                         input$ressd^2),
                                       ncol = 2)
    sample_distribution <- MASS::mvrnorm(n = sample_size,
                                         mu = sample_meanvector, 
                                         Sigma = sample_covariance_matrix)
    dat2 <- dat %>% mutate(y = cos(2*pi*(t/365-input$phaseshift/2))+sample_distribution[,2],
                           x = input$r1*cos(2*pi*t/365)+sample_distribution[,1]) 
    fig1 <- dat2 %>% 
      pivot_longer(cols = c("x","y")) %>% 
      filter(year(date)%in% 2010:2014) %>% 
      ggplot(aes(x=date, y = value, col = name)) + geom_line()+
      labs(title = "Time plot")+
      scale_x_date(date_breaks = "1 year",
                   date_labels = "%Y") + 
      theme(legend.position = c(.9,.9),
            axis.title = element_blank(),
            legend.title = element_blank())
    fig2 <- dat2 %>% 
      ggplot(aes(x=x, y = y)) + geom_point() + geom_smooth(method = "lm")+
      labs(title = "Scatterplot")
    ggpubr::ggarrange(fig1,fig2,ncol = 2)
  })
  
  output$portfolioPlot <- renderPlot({
    set.seed(2134)
    sample_covariance_matrix <- matrix(c(input$ressd^2,
                                         input$ressd^2*input$rescor,
                                         input$ressd^2*input$rescor, 
                                         input$ressd^2),
                                       ncol = 2)
    seasonal <- matrix(c(input$r1^2, input$r1*cos(pi*input$phaseshift),input$r1*cos(pi*input$phaseshift),1^2),
                       ncol=2, nrow=2)/2
    sample_distribution <- MASS::mvrnorm(n = sample_size,
                                         mu = sample_meanvector, 
                                         Sigma = sample_covariance_matrix)
    dat2 <- dat %>% mutate(y = cos(2*pi*(t/365-input$phaseshift/2))+sample_distribution[,2],
                           x = input$r1*cos(2*pi*t/365)+sample_distribution[,1]) 
    COV <- cov(dat2[,c("x","y")])
    #COV <- seasonal + sample_covariance_matrix
    
    weights <- rbind(
      "Empirical"= MPT(mu=c(1,1), cov = COV, target = 1, return_value = FALSE),
      "Seasonal" = MPT(c(1,1), cov = seasonal, target = 1, return_value = FALSE),
      "Season-adj"=MPT(c(1,1), cov = sample_covariance_matrix, target = 1, return_value = FALSE))[,3]
    sample_distribution <- MASS::mvrnorm(n = sample_size,
                                         mu = sample_meanvector, 
                                         Sigma = sample_covariance_matrix)
    dat2 <- dat %>% mutate(y = cos(2*pi*(t/365-input$phaseshift/2))+sample_distribution[,2],
                           x = input$r1*cos(2*pi*t/365)+sample_distribution[,1]) %>% 
      mutate("Empirical" = weights[1]*x+(1-weights[1])*y,
             "Seasonal" = weights[2]*x+(1-weights[2])*y,
             "Season-adj" = weights[3]*x+(1-weights[3])*y) %>% 
      select(-x,-y, -t) %>% 
      pivot_longer(cols = 2:4) %>% 
      mutate(name = factor(name, levels = c("Empirical",
                                            "Seasonal",
                                            "Season-adj"))) %>% 
      as_tsibble(index = date, key = name)
    fig1 <- dat2 %>% 
      filter(year(date)%in% 2010:2014) %>% 
     ggplot(aes(x=date, y = 100+value))+geom_line()+facet_wrap(~name, ncol = 1)+
      theme(legend.position = "none",
            strip.text = element_blank())+
      scale_x_date(date_breaks = "1 year", name = "Time")+
      ylab("Portfolio value")
    fig2 <- dat2 %>% ACF(value, lag_max =400) %>% autoplot()+
      scale_x_cf_lag(breaks = seq(0,400,50))
    ggpubr::ggarrange(fig1,fig2, ncol = 2)  
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
