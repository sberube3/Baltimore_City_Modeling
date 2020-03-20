rm(list=ls())

setwd('~/Dropbox/COVID-BaltimoreCity//')

require(xlsx)
require(magrittr)
require(stringr)
require(reshape2)
require(dplyr)
require(ggplot2)

theme_set(theme_classic(base_size=12))

## load state demog data and clean it
#data <- read.xlsx('MO_demographic_2017.xlsx',sheetIndex = 1)
data <- read.csv('MO_demographic_2017.csv')

data$Estimate %>% as.character() %>% str_remove( ",") %>%
  str_remove(' ') %>% str_remove(',') %>% as.numeric() -> data$Estimate

data$Margin_of_Error %>% as.character() %>%
  str_remove("\\+") %>% str_remove('\\/') %>%
  str_remove('\\-') %>% str_remove(',') %>% str_remove(' ') %>% as.numeric() -> data$Margin_of_Error

data$Percent %>% as.character() %>% as.numeric() -> data$Percent

data$Percent_Margin_of_Error %>% as.character() %>%
  str_remove("\\+") %>% str_remove('\\/') %>% str_remove('\\-') %>%
  str_remove(',') %>% str_remove(' ') %>% as.numeric() -> data$Percent_Margin_of_Error


data[5:17,] %>%
  mutate(Variable = factor(Variable,levels = data[5:17,]$Variable)) %>%
  ggplot(aes(Variable,Estimate))+geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin=Estimate - Margin_of_Error,ymax=Estimate+Margin_of_Error)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab(NULL)

## Pull out population data
MO.pop <- data[5:17,]$Estimate

## aggregate the data such that oldest class is 75 and up
## this is to match with the socialmixr matrix
MO.pop[12] <- MO.pop[12] + MO.pop[13]
MO.pop <- MO.pop[1:12]
Ncomp <- length(MO.pop)

## setup polymod matrix
## use age classes in the data
## for now, hard code it into the function, but can change later
## use the UK mixing pattern as I don't believe US is available
library(socialmixr)
data(polymod)
W <- contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0,5,10,15,20,25,35,45,55,60,65,75,85,90))$matrix

## test cases to check behavior of R0
## sanity checks in the case of equal mixing and equal populations
# MO.pop <- rep(mean(MO.pop),Ncomp)
# W <- waifw <- matrix(.4,Ncomp,Ncomp)

## we want the max eigen value of our input matrix to be 1 such that we can say R0 is just beta / gamma
## the following formula gives the R0 calculation for age strucutred matrix models
## ref : http://www.sherrytowers.com/towers_feng_2012.pdf page 242 right side (matrix is called C_ij)

A <- matrix(0,Ncomp,Ncomp)

for (i in 1:Ncomp){
  for (j in 1:Ncomp){
    A[i,j] <- W[i,j]*MO.pop[i]/MO.pop[j]
  }
}

## compute spectral radius / scaling parameter
r <- eigen(A)
lam <- r$values
alpha <- max(Re(lam))

## now check to make sure this worked with a dummy matrix
A1 <- matrix(0,Ncomp,Ncomp)

for (i in 1:Ncomp){
  for (j in 1:Ncomp){
    A1[i,j] <- W[i,j]*MO.pop[i]/MO.pop[j]*(1/alpha)
  }
}

r.1 <- eigen(A1)
lam.1 <- r.1$values
specrad.1 <- max(Re(lam.1))
stopifnot(all.equal(specrad.1,1))

## now the matrix is rescaled have R0  = 1, so beta0 can be scaled to be real transmission
W <- W / alpha

## set initial conditions
## I'm not entirely sure what to put for these currently...
ICs <- c(
  S = MO.pop * 1,
  A = rep(100,length(MO.pop)),
  I = rep(1,length(MO.pop)),
  R = MO.pop,
  incid = rep(0,Ncomp)
)

## set the R compartment
## crudely just set anything negative to be zero but we can fine tune this more depending on what we assume S0 does
ICs[(3*Ncomp+1):(4*Ncomp)] <- MO.pop -  ICs[1:Ncomp] - ICs[(Ncomp+1):(2*Ncomp)] - ICs[(2*Ncomp+1):(3*Ncomp)]
ICs[(3*Ncomp+1):(4*Ncomp)] [ICs[(3*Ncomp+1):(4*Ncomp)]  < 0 ] <- 0

## prop_serious : proportion of cases requiring serious attention / hospitalization
## should this be age dependent (probably), but how?
## this is a good assumption to plot in our diagnostics.
prop_serious <- 0.19

## set prop_symtomatic
## right now just set it to be 5% but
prop_symptomatic <- 0.40

## declare model parameters
## all model parameters should be in units of days!

## population sizes by demographic data
N <- MO.pop

## units in days
## infectious period
gamma <- 1/6.5

## R0
R0 <- 2.5

## set beta based on that value
beta0 <- R0*gamma

## seasonal forcing should be modest here
beta1 <- 0.1

## when should seasonal forcing peak?
phase <- 0

## set births to be zero currently
mu <- 0

## set natural death rate to be zero currently
v <- 0

## now check to make sure the R0 we get is the R0 we put in
## using the same formula as above
R0.mat <- matrix(0,Ncomp,Ncomp)

for (i in 1:Ncomp){
  for (j in 1:Ncomp){
    R0.mat[i,j] <- W[i,j]*MO.pop[i]/MO.pop[j]*beta0 / gamma
  }
}

eigen(R0.mat)$values[1]

## how long to run the model for?
## currently set to be 50 days integrated at by hour
delta.t <- 1/24
time <- seq(1,20,by = delta.t)

seir_step <- function(stoch = F) {

  ## set up a matrix to store values in by variable and time
  ## each X[it,] is the variable at one hour
  x <- matrix(NA,length(time),Ncomp * 5)

  x[1,] <- round(ICs)

  ## susceptible individuals
  S <- x[,1:Ncomp];
  ## asypm individuals
  A <- x[,(Ncomp+1):(2*Ncomp)];
  ## symp individuals
  I <- x[,(2*Ncomp+1):(3*Ncomp)];
  ## recovered individuals
  R <- x[,(3*Ncomp+1):(4*Ncomp)]
  ## incidence
  incid <- x[,(4*Ncomp+1):(5*Ncomp)];

  ## seasonal transmission
  seas <- beta0 * (1 + beta1 * cos(2 * pi * time/365 - phase))

  for(it in 1:(length(time) - 1)){

    WI <- W%*%(A[it,] + I[it,])
    WI[!is.finite(WI)] <- 0

    births <-rep(0,Ncomp)
    births[1] <- mu
    deaths <- rep(v,Ncomp)

    ## add stochasticity to FOI
    if(stoch == T){
      dw <- rnorm(Ncomp,mean = 1, sd = 0.05)
    }else{
      dw <- 1
    }

    ## declare transitions in model
    foi_prob <- 1 - exp( - seas[it] * WI/N * dw * delta.t)
    inf_prob <- 1 - exp( - gamma * delta.t)
    death_prob <- 1 - exp( - deaths * delta.t)

    ## stochastic formulation of the model
    if(stoch == T){

      new_inf <- rbinom(n = Ncomp, size = round(S[it,]), prob = foi_prob)
      new_rec_A <- rbinom(n = Ncomp, size = round(A[it,]), prob = inf_prob)
      new_rec_I <- rbinom(n = Ncomp, size = round(I[it,]), prob = inf_prob)

      S[it + 1, ] <- S[it,] +  births*delta.t - new_inf - rbinom(n = Ncomp, size = round(S[it]), prob = death_prob)
      A[it + 1, ] <- A[it,] +  (1 - prop_symptomatic) * new_inf - new_rec_A - rbinom(n = Ncomp, size = round(A[it]), prob = death_prob)
      I[it + 1, ] <- I[it,] +  prop_symptomatic * new_inf - new_rec_I - rbinom(n = Ncomp, size = round(I[it]), prob = death_prob)
      R[it + 1, ] <- R[it,] +  new_rec_I + new_rec_A - rbinom(n = Ncomp, size = round(R[it]), prob = death_prob)

      incid[it, ] <-  new_rec_I + new_rec_A

    }
    ## deterministic equations to check
    if(stoch == F){

      S[it + 1, ] <- S[it,] + delta.t * (births - seas[it] * WI * S[it,] * dw / N  - deaths*S[it,])
      A[it + 1, ] <- A[it,] + delta.t * ( (1 - prop_symptomatic) *  seas[it] * WI * S[it,] * dw / N - A[it,]*(gamma - deaths))
      I[it + 1, ] <- I[it,] + delta.t * (  prop_symptomatic *  seas[it] * WI * S[it,] * dw / N - I[it,]*(gamma - deaths))
      R[it + 1, ] <- R[it,] + delta.t * (A[it,]*gamma+ I[it,]*gamma - R[it,]* deaths)

      incid[it,] <-  delta.t*(A[it,]*gamma+ I[it,]*gamma )

    }

  }

  out <- data.frame(cbind(time,S,A,I,R,incid))
  names(out) <- c('time',names(ICs))
  return(out)
}

## how many simulations to run?
nsim <- 10
for(n in 1:nsim){

  ## run the simulation one time
  single.sim <- seir_step(stoch = T)

  single.sim %>%
    ggplot(aes(time,I5))+geom_line()

  single.sim %>%
    dplyr::select(paste0('I',1:Ncomp)) %>%
    rowSums() -> totalI

  single.sim %>%
    dplyr::select(paste0('incid',1:Ncomp)) %>%
    rowSums() -> totalincid

  daily_incid <- unname(tapply(totalincid, (seq_along(totalincid)-1) %/% (1/delta.t), sum))

  cases_requiring_attention <- totalI * prop_serious

  daily_cases_requiring_attention <- unname(tapply(cases_requiring_attention, (seq_along(cases_requiring_attention)-1) %/% (1/delta.t), sum))

  dailytime <- single.sim$time[seq(1,nrow(single.sim),(1/delta.t))]

  single.daily.data <- data.frame('days'=dailytime,'daily_cases_requiring_attention'=daily_cases_requiring_attention)
  single.daily.data$sim <- n

  if(n == 1){
    daily.data <- single.daily.data
  }else{
    daily.data <- rbind(daily.data,single.daily.data)
  }
}

daily.data %>%
  ggplot(aes(days,daily_cases_requiring_attention,color=factor(sim)))+
  geom_line()+
  xlim(0,max(time) - 1)

