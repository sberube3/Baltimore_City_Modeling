rm(list=ls())
### code adapted from Alex Becker's COVID age structured discrete time model.
### This code will incorporate age structured mixing between groups and produce estimates of the number of 
### infected individuals per day (and cumulative incident infections). We have adapted the code to add in 
### separate healthcare workers and homeless populations. 

#setwd('~/Dropbox/COVID-BaltimoreCity//')

require(socialmixr)
require(magrittr)
require(stringr)
require(reshape2)
require(dplyr)
require(ggplot2)
require(truncnorm)

### Functions: this document has 5 different functions to clean up the age structured data, make a mixing matrix, 
### rescale the mixing matrix, run the discrete time simulation, and set up all of the parameters/etc. for the discrete time simulation.
# make_age_structure_matrix(age_data_file, homeless_n, healthcare_n)
# make_polymod_matrix()
# rescale_age_matrix(Ncomp, W, BC_pop)
# sair_step(stoch = F, Ncomp, ICs, params, time, delta.t)
# setup_seir_model(stoch)

### general things
# BC_pop: the population of Baltimore by age
# W: age specific mixing rates
# C: an additional mixing matrix that can be included so mixing for healthcare workers does not decrease -- 
#    currently if you set it all to 1 then it is the same as turning it off
# Ncomp: number of compartments (ages + healthcare + homeless)
# SAIR with asymtompatics separate from symtomatics
# currently there is a bit to estimate severity that will not be needed later
# gamma <- 1/6.5 ## infectious period? make a range?
# R0 <- 2.5 ## make a range?
# prop_symptomatic <- 0.20 ## will need to update/change

### age categories (Baltimore):
# Under 5 years    41152
# 5 to 9 years    35441
# 10  to 14 years    34339
# 15 to 19 years    44278
# 20 to 24 years    56460
# 25 to 34 years   103564
# 35 to 44 years    76564
# 45 to 54 years    87445
# 55 to 59 years    37978
# 60 to 64 years    30928
# 65 to 74 years    38552
# 75+   23910 + 10350

theme_set(theme_classic(base_size=12))

## ---- loading demographic data ---- ####
  ### load [location] age structured demographic data. This code assumes that the first row will be the total 
  ### population and the following rows are the ages broken down by age groups.
  ## If you do not want homeless or healthcare workers, just set these values to zero
  ## assumes the additional of healthcare workers and homeless does not decrease from their respective age group 

  make_age_structure_matrix <- function(age_data_file, homeless_n, healthcare_n){
    ## Pull out population data
    n_age_classes <- nrow(age_data_file)-1
    BC_pop <- age_data_file[2:nrow(age_data_file),]$Estimate
    
    ## aggregate the data such that oldest class is 75 and up
    ## this is to match with the socialmixr matrix
    BC_pop[(n_age_classes-1)] <- BC_pop[(n_age_classes-1)] + BC_pop[n_age_classes]
    BC_pop <- BC_pop[1:(n_age_classes-1)]
    BC_pop <- c(BC_pop, healthcare_n) ## healthcare workers
    BC_pop <- c(BC_pop, homeless_n) ## number homeless
    
    Ncomp <- length(BC_pop) 
    return(list(Ncomp = Ncomp, BC_pop = BC_pop))  
  }

## ---- making polymod matrix ---- ####
  ### sets up polymod matrix using data from the UK (we could not find any US polymod data). It will make a 
  ### mixing matrix set up for all of the different age classes and takes the mixing patterns for the 45-55 age 
  ### group for the healthcare workers and homeless individuals

  make_polymod_matrix <- function(age.limits=c(0,5,10,15,20,25,35,45,55,60,65,75,85,90), 
                                  hcw.mix="45,55", hml.mix="45,55"){
    ## setup polymod matrix
    ## use age classes in the data
    ## for now, hard code it into the function, but can change later
    ## use the UK mixing pattern as I don't believe US is available
    data(polymod)
    age_mix <- contact_matrix(polymod, countries = "United Kingdom", age.limits = age.limits)$matrix
    
    ## assumes homeless and healthcare workers have the same mixing as 45-55 year olds
    W <- matrix(NA,nrow(age_mix)+2, ncol(age_mix)+2)
    
    rownames(W) <- c(rownames(age_mix), 'healthcare', 'homeless')
    colnames(W) <- c(colnames(age_mix), 'healthcare', 'homeless')
    W[1:nrow(age_mix),1:ncol(age_mix)] <- age_mix
    
    ## test case of uniform contact
    ## W[1:nrow(age_mix),1:ncol(age_mix)] <- 1
    
    W[nrow(age_mix)+1,] = W[grep(hcw.mix, rownames(W)),] ## healthcare
    W[nrow(age_mix)+2,] = W[grep(hml.mix, rownames(W)),] ## homeless
    W[,nrow(age_mix)+1] = W[,grep(hcw.mix, rownames(W))] ## healthcare
    W[,nrow(age_mix)+2] = W[,grep(hml.mix, rownames(W))] ## homeless
    return(W)
  }

## ---- rescaling age matrix ---- ####
  ## we also need to adjust our matrix to make sure R0 = beta/gamma
  ## we want the max eigen value of our input matrix to be 1 such that we can say R0 is just beta / gamma
  ## the following formula gives the R0 calculation for age strucutred matrix models
  ## ref : http://www.sherrytowers.com/towers_feng_2012.pdf page 242 right side (matrix is called C_ij)

  rescale_age_matrix <- function(Ncomp, W, BC_pop, c_scale_vec){
    A <- matrix(0,Ncomp,Ncomp)
    for (ii in 1:Ncomp){
      for (jj in 1:Ncomp){
        A[ii,jj] <- W[ii,jj]*BC_pop[ii]/BC_pop[jj]
      }
    }
    ## compute spectral radius / scaling parameter
    r <- eigen(A)
    lam <- r$values
    alpha <- max(Re(lam))
    W <- W / alpha
    ## now the matrix is rescaled have R0  = 1, so beta0 can be scaled to be real transmission
    C <- matrix(c_scale_vec,nrow(W),ncol(W)) ## a special contact matrix to be used to rescale health facility worker contact rates - when set to 1, it is turned off 
    return(list(W = W, C = C))
  }

## ---- main S(A)IR function ---- ####
  sair_step <- function(stoch = F, stoch.init = F, sNcomp, ICs, params, time, delta.t){
    C = params$C
    W = params$W
    beta0 = params$beta0
    beta1 = params$beta1
    phase = params$phase
    mu = params$mu
    v = params$v
    N=params$N
    sigma = params$sigma
    gamma=params$gamma
    prop_symptomatic=params$prop_symptomatic
    sd.dw = params$sd.dw
    
    ## set up a matrix to store values in by variable and time
    ## each X[it,] is the variable at one hour
    x <- matrix(NA,length(time),Ncomp * 7)
  
    ## set initial conditions
    if(stoch.init){
      Ninit <- sample(10:60, 1)
      Ninit_byage <- rmultinom(1, Ninit, prob = N/sum(N))[,1]
      Ninit_asy <- round(Ninit_byage * (1-prop_symptomatic))
      Ninit_sym <- Ninit_byage - Ninit_asy
      ICs <- c(S = N, 
               E = rep(0, Ncomp),
               A = Ninit_asy,
               I = Ninit_sym,
               R = rep(0, Ncomp),
               incid_A = rep(0, Ncomp),
               incid_I = rep(0, Ncomp))
      x[1,] <- round(ICs)
      
    }else{ x[1,] <- round(ICs) }
  
    S <- x[,1:Ncomp]; ## susceptible individuals
    E <- x[,(Ncomp+1):(2*Ncomp)]; ## exposed individuals 
    A <- x[,(2*Ncomp+1):(3*Ncomp)]; ## asymptomatic individuals
    I <- x[,(3*Ncomp+1):(4*Ncomp)];## symp individuals
    R <- x[,(4*Ncomp+1):(5*Ncomp)] ## recovered individuals
    
    ## incidence
    incid_A <- x[,(5*Ncomp+1):(6*Ncomp)];
    incid_I <- x[,(6*Ncomp+1):(7*Ncomp)];
    
    ## seasonal transmission
    seas <- beta0 * (1 + beta1 * cos(2 * pi * time/365 - phase))
    
    for(it in 1:(length(time) - 1)){
    #  WI <- C%*%W%*%(A[it,] + I[it,])
      WI <- (C*W)%*%(A[it,] + I[it,])
      
      WI[!is.finite(WI)] <- 0
      births <-rep(0,Ncomp)
      births[1] <- mu
      deaths <- rep(v,Ncomp)
      
      ## add stochasticity to FOI
      if(stoch == T){
        dw <- rtruncnorm(Ncomp, a=0, mean = 1, sd = sd.dw)
      }else{
        dw <- 1
      }
      ## declare transitions in model
      foi_prob <- 1 - exp( - seas[it] * WI/N * dw * delta.t)
      exposed_prob <- 1 - exp( - sigma * delta.t)
      inf_prob <- 1 - exp( - gamma * delta.t)
      death_prob <- 1 - exp( - deaths * delta.t)
      
      ## stochastic formulation of the model
      if(stoch == T){
        new_exp <- rbinom(n = Ncomp, size = round(S[it,]), prob = foi_prob)
        new_inf <- rbinom(n = Ncomp, size = round(E[it,]) , prob = exposed_prob)
        new_infA <- round( (1-prop_symptomatic)*new_inf )
        new_infI <- new_inf - new_infA
        new_rec_A <- rbinom(n = Ncomp, size = round(A[it,]), prob = inf_prob)
        new_rec_I <- rbinom(n = Ncomp, size = round(I[it,]), prob = inf_prob)
        
        S[it + 1, ] <- S[it,] +  births*delta.t - new_exp - rbinom(n = Ncomp, size = round(S[it,]), prob = death_prob)
        E[it + 1, ] <- E[it,] +  new_exp - new_inf - rbinom(n = Ncomp, size = round(E[it,]), prob = death_prob )
        A[it + 1, ] <- A[it,] +  new_infA - new_rec_A - rbinom(n = Ncomp, size = round(A[it,]), prob = death_prob)
        I[it + 1, ] <- I[it,] +  new_infI - new_rec_I - rbinom(n = Ncomp, size = round(I[it,]), prob = death_prob)
        R[it + 1, ] <- R[it,] +  new_rec_I + new_rec_A - rbinom(n = Ncomp, size = round(R[it,]), prob = death_prob)
        
        ## make incidence the new number of daily individuals becoming infected
        incid_A[it, ] <- new_infA
        incid_I[it, ] <- new_infI
      }
      
      ## deterministic equations
      if(stoch == F){
        S[it + 1, ] <- S[it,] + delta.t * (births - seas[it] * WI * S[it,] * dw / N  - deaths*S[it,])
        E[it + 1, ] <- E[it,] + delta.t * (seas[it] * WI * S[it,] * dw / N - deaths*E[it,] - sigma*E[it,])
        A[it + 1, ] <- A[it,] + delta.t * ( (1 - prop_symptomatic)*sigma*E[it,]   - A[it,]*(gamma - deaths))
        I[it + 1, ] <- I[it,] + delta.t * (  prop_symptomatic*sigma*E[it,] - I[it,]*(gamma - deaths) )
        R[it + 1, ] <- R[it,] + delta.t * (A[it,]*gamma+ I[it,]*gamma - R[it,]* deaths)
        incid_A[it,] <-  (1 - prop_symptomatic)*(seas[it] * WI * S[it,] * dw / N)
        incid_I[it,] <- prop_symptomatic*(seas[it] * WI * S[it,] * dw / N)
      }
    }
    out <- data.frame(cbind(time,S,E,A,I,R,incid_A, incid_I))
    names(out) <- c('time',names(ICs))
    ## output is the number in each class per time point per age-category+homeless+healthcare workers
    return(out)
  }

## ---- main SEAIR function with variable R0 input ---- ####
sair_step_variableR0 <- function(stoch = F, stoch.init = F, R0vec, Ncomp, ICs, params, time, delta.t, init.min = 10, init.max=60, init.dist = NULL, c_scale_mat=NULL){
  
  C = params$C
  W = params$W
  beta0 = params$beta0
  beta1 = params$beta1
  phase = params$phase
  mu = params$mu
  v = params$v
  N=params$N
  sigma = params$sigma
  gamma=params$gamma
  prop_symptomatic=params$prop_symptomatic
  sd.dw <- params$sd.dw
  
  ## set up a matrix to store values in by variable and time
  ## each X[it,] is the variable at one hour
  x <- matrix(NA,length(time),Ncomp * 7)
  
  ## set initial conditions
  if(stoch.init){
    Ninit <- sample(init.min:init.max, 1)
    if(is.null(init.dist)){pinit <- N / sum(N)}else{pinit <- init.dist}
    Ninit_byage <- rmultinom(1, Ninit, prob = pinit)[,1]
    Ninit_asy <- round(Ninit_byage * (1-prop_symptomatic))
    Ninit_sym <- Ninit_byage - Ninit_asy
    ICs <- c(S = N, 
             E = rep(0, Ncomp),
             A = Ninit_asy,
             I = Ninit_sym,
             R = rep(0, Ncomp),
             incid_A = rep(0, Ncomp),
             incid_I = rep(0, Ncomp))
    x[1,] <- round(ICs)
    
  }else{ x[1,] <- round(ICs) }
  
  S <- x[,1:Ncomp]; ## susceptible individuals
  E <- x[,(Ncomp+1):(2*Ncomp)]; ## exposed individuals 
  A <- x[,(2*Ncomp+1):(3*Ncomp)]; ## asymptomatic individuals
  I <- x[,(3*Ncomp+1):(4*Ncomp)];## symp individuals
  R <- x[,(4*Ncomp+1):(5*Ncomp)] ## recovered individuals
  
  ## incidence
  incid_A <- x[,(5*Ncomp+1):(6*Ncomp)]
  incid_I <- x[,(6*Ncomp+1):(7*Ncomp)]
  
  # recalculate beta0 for each R0 value
  beta0 <- R0vec * (gamma + v) * (sigma + v) / sigma 
  seas <- beta0 * (1 + beta1 * cos(2 * pi * time/365 - phase))
  R0 <- vector(length=length(time))
  Re <- vector(length=length(time))
  
  for(it in 1:(length(time) - 1)){
    
    # calculate and store the R0 value for time-specific beta0 - to confirm the new value is correct
    R0.mat <- matrix(0,Ncomp,Ncomp)
    for (i in 1:Ncomp){
      for (j in 1:Ncomp){
        R0.mat[i,j] <- W[i,j]*N[i]/N[j]* beta0[it] * sigma / ( (sigma + v) * (v + gamma))
      }
    }
    R0[it] <- Re(eigen(R0.mat)$values[1])
    
    # proceed with model step, as in sair_step
    # if c_scale_mat is provided, recalculate C matrix
    if(!is.null(c_scale_mat)){
      C <- matrix(c_scale_mat[it,],nrow(W),ncol(W))
    }
    WI <- (C*W)%*%(A[it,] + I[it,])
    WI[!is.finite(WI)] <- 0
    
    # calculate and store the effective R value (if any value of C is non-1)
    if(sum(C!=1)>0){
      Wscal <- C*W
      Re.mat <- matrix(0,Ncomp,Ncomp)
      for (i in 1:Ncomp){
        for (j in 1:Ncomp){
          Re.mat[i,j] <- Wscal[i,j]*N[i]/N[j]* beta0[it] * sigma / ( (sigma + v) * (v + gamma))
        }
      }
      Re[it] <- Re(eigen(Re.mat)$values[1])
    }else{
      Re[it] <- R0[it]
    }
    
    
    births <-rep(0,Ncomp)
    births[1] <- mu
    deaths <- rep(v,Ncomp)
    
    ## add stochasticity to FOI
    if(stoch == T){
      dw <- rtruncnorm(Ncomp, a=0, mean = 1, sd = sd.dw)
    }else{
      dw <- 1
    }
    
    ## declare transitions in model
    foi_prob <- 1 - exp( - seas[it] * WI/N * dw * delta.t)
    exposed_prob <- 1 - exp( - sigma * delta.t)
    inf_prob <- 1 - exp( - gamma * delta.t)
    death_prob <- 1 - exp( - deaths * delta.t)
    
    ## stochastic formulation of the model
    if(stoch == T){
      new_exp <- rbinom(n = Ncomp, size = round(S[it,]), prob = foi_prob)
      new_inf <- rbinom(n = Ncomp, size = round(E[it,]) , prob = exposed_prob)
      new_infA <- round( (1-prop_symptomatic)*new_inf )
      new_infI <- new_inf - new_infA
      new_rec_A <- rbinom(n = Ncomp, size = round(A[it,]), prob = inf_prob)
      new_rec_I <- rbinom(n = Ncomp, size = round(I[it,]), prob = inf_prob)
      
      S[it + 1, ] <- S[it,] +  births*delta.t - new_exp - rbinom(n = Ncomp, size = round(S[it,]), prob = death_prob)
      E[it + 1, ] <- E[it,] +  new_exp - new_inf - rbinom(n = Ncomp, size = round(E[it,]), prob = death_prob )
      A[it + 1, ] <- A[it,] +  new_infA - new_rec_A - rbinom(n = Ncomp, size = round(A[it,]), prob = death_prob)
      I[it + 1, ] <- I[it,] +  new_infI - new_rec_I - rbinom(n = Ncomp, size = round(I[it,]), prob = death_prob)
      R[it + 1, ] <- R[it,] +  new_rec_I + new_rec_A - rbinom(n = Ncomp, size = round(R[it,]), prob = death_prob)
      
      ## make incidence the new number of daily individuals becoming infected
      incid_A[it, ] <- new_infA
      incid_I[it, ] <- new_infI
    }
    
    ## deterministic equations to check -- does not currently work 
    if(stoch == F){
      S[it + 1, ] <- S[it,] + delta.t * (births - seas[it] * WI * S[it,] * dw / N  - deaths*S[it,])
      E[it + 1, ] <- E[it,] + delta.t * (seas[it] * WI * S[it,] * dw / N - deaths*E[it,] - sigma*E[it,])
      A[it + 1, ] <- A[it,] + delta.t * ( (1 - prop_symptomatic)*sigma*E[it,]   - A[it,]*(gamma - deaths))
      I[it + 1, ] <- I[it,] + delta.t * (  prop_symptomatic*sigma*E[it,] - I[it,]*(gamma - deaths) )
      R[it + 1, ] <- R[it,] + delta.t * (A[it,]*gamma+ I[it,]*gamma - R[it,]* deaths)
      incid_A[it,] <-  (1 - prop_symptomatic)*(seas[it] * WI * S[it,] * dw / N)
      incid_I[it,] <- prop_symptomatic*(seas[it] * WI * S[it,] * dw / N)
    }
  }
  
  out <- data.frame(cbind(time,S,E,A,I,R,incid_A, incid_I,R0, Re))
  names(out) <- c('time',names(ICs),"R0", "Reff")
  ## output is the number in each class per time point per age-category+homeless+healthcare workers
  return(out)
  
}

## ---- function to set up and organize the mixing data, inital conditions, parameters, etc. ---- ####

  setup_seir_model <- function(stoch, R0, c_scale_vec,
                               gamma=1/6.5, sigma=1/5.2,
                               phase=0, beta1=0, mu=0, v=0,
                               prop_symptomatic, sd.dw=0.05,
                               hcw.mix="45,55", hml.mix="45,55",
                               homeless_n=2000, healthcare_n=5000){
    
    data <- read.csv('Baltimore_AgePopulation.csv')
    age_data <- make_age_structure_matrix(data, homeless_n = homeless_n, healthcare_n = healthcare_n)
    BC_pop = age_data$BC_pop
    Ncomp = age_data$Ncomp
    W <- make_polymod_matrix(hcw.mix = hcw.mix, hml.mix = hml.mix)
    rescale_mixing <- rescale_age_matrix(Ncomp, W, BC_pop, c_scale_vec)
    W <- rescale_mixing$W
    C <- rescale_mixing$C
    
    ## set prop_symtomatic
    prop_symptomatic <- prop_symptomatic
    if(length(prop_symptomatic)!=Ncomp){warning("prop_symptomatic is the wrong length")}
  
    ## set initial conditions
    ICs <- c(S = BC_pop, 
             E = rep(0,length(BC_pop)),
             A = rep(8,length(BC_pop)),
             I = rep(2,length(BC_pop)), 
             R = BC_pop,
             incid_A = rep(0,Ncomp),
             incid_I = rep(0,Ncomp))
    
    ## set the R compartment
    ## crudely just set anything negative to be zero but we can fine tune this more depending on what we assume S0 does
    ICs[(4*Ncomp+1):(5*Ncomp)] <- BC_pop -  ICs[1:Ncomp] - ICs[(Ncomp+1):(2*Ncomp)] - ICs[(2*Ncomp+1):(3*Ncomp)] - ICs[(3*Ncomp+1):(4*Ncomp)]
    ICs[(4*Ncomp+1):(5*Ncomp)] [ICs[(4*Ncomp+1):(5*Ncomp)]  < 0 ] <- 0
    
    ## population sizes by demographic data
    N <- BC_pop
    
    ## units in days!! 
    gamma <- gamma ## infectious period
    sigma <- sigma ## latent period
    phase <- phase ## when should seasonal forcing peak?
    mu <- mu ## set births to be zero currently
    v <- v ## set natural death rate to be zero currently
    #R0 <- 2.5 ## make a range?   ## R0 = beta * sigma / ((sigma + v) * (v + gamma)) for SEAIR model with split proportion into A-I and only A and I contributing to infection
    beta0 <- R0 * (gamma + v) * (sigma + v) / sigma ## set beta based on that value
    beta1 <- beta1 ## seasonal forcing should be modest here
    sd.dw <- sd.dw
    
    ## now check to make sure the R0 we get is the R0 we put in
    ## using the same formula as above
    R0.mat <- matrix(0,Ncomp,Ncomp)
    for (i in 1:Ncomp){
      for (j in 1:Ncomp){
        R0.mat[i,j] <- W[i,j]*BC_pop[i]/BC_pop[j]* beta0 * sigma / ( (sigma + v) * (v + gamma))
      }
    }
    # for (i in 1:Ncomp){
    #   for (j in 1:Ncomp){
    #     R0.mat[i,j] <- W[i,j]*BC_pop[i]/BC_pop[j]* beta0 * (gamma + v) * (sigma +v )/sigma
    #   }
    # }
    print(eigen(R0.mat)$values[1]) ## just a check
    return(list(C = C, W = W, beta0 = beta0, beta1 = beta1, 
                phase = phase, mu = mu, v = v, ICs = ICs, 
                Ncomp = Ncomp, N=N, gamma=gamma, sigma = sigma,
                prop_symptomatic=prop_symptomatic, sd.dw=sd.dw))
  }

## ---- Setting up a run ---- ####
  ## how long to run the model for?
  ## currently set to be 300 days integrated at the day
  # delta.t <- 1/1
  # time <- seq(1,300,by = delta.t)
  # 
  # prop_symptomatic <- c(0.141, 0.106, 0.074, 0.184, 0.293, 0.387, 0.438, 
  #                       0.535, 0.693, 0.816, 0.765, 0.749, 0.535, 0.535)
  # 
  # all_prelim_info <- setup_seir_model(stoch = TRUE, R0 = 2.0, c_scale_vec = 0.3)
  # ## what is printed should be R0 
  # Ncomp = all_prelim_info$Ncomp
  # ICs = all_prelim_info$ICs
  # 
  # params = list(C = all_prelim_info$C, 
  #               W = all_prelim_info$W, 
  #               beta0 = all_prelim_info$beta0, 
  #               beta1 = all_prelim_info$beta1, 
  #               phase = all_prelim_info$phase, 
  #               mu = all_prelim_info$mu, 
  #               v = all_prelim_info$v, 
  #               N=all_prelim_info$N, 
  #               gamma=all_prelim_info$gamma, 
  #               sigma = all_prelim_info$sigma, 
  #               prop_symptomatic=all_prelim_info$prop_symptomatic, 
  #               sd.dw=all_prelim_info$sd.dw)

  ## running some different simulations and making some basic plots
  # single.sim <- sair_step(stoch = TRUE, Ncomp, ICs, params, time, delta.t)
  # single.sim %>% ggplot(aes(time,I5))+geom_line()
  # single.sim %>% dplyr::select(paste0('I',1:Ncomp)) %>% rowSums() -> totalI
  # single.sim %>% dplyr::select(paste0('incid_I',1:Ncomp)) %>% rowSums() -> totalincid_I
  # single.sim %>% dplyr::select(paste0('incid_A',1:Ncomp)) %>% rowSums() -> totalincid_A
  # single.sim <- cbind(single.sim, totalincid_A)
  # single.sim <- cbind(single.sim, totalincid_I)
  # par(mfrow=c(1,2))
  # plot(single.sim$totalincid_A,type='l')
  # plot(single.sim$totalincid_I,type='l')
  # 
  # test_sim <- sair_step(stoch = TRUE, Ncomp, ICs, params, time, delta.t)
  # run_index = rep(n, nrow(test_sim))
  # test_sim <- cbind(run_index, test_sim)
  # all_sim <- rbind(all_sim, test_sim)
  # 
  # run_index = rep(0, nrow(test_sim))
  # test_sim <- cbind(run_index, test_sim)
  
## ---- running multiple simulations ---- ####     
  # nsim <- 500
  # start_index <- seq(1, nsim*length(time)+1, by = length(time))
  # all_sim <- matrix(NA,1,(Ncomp*7)+2)
  # colnames(all_sim) <- c('run_index', colnames(test_sim))#matrix(,nsim*length(time),(Ncomp*5)+2) ## +2 -> time step, run index
  # Sys.time()
  # for(n in 1:nsim){
  # #  prop_serious <- 0.2 ### will update later
  #   ## run the simulation one time
  #   single.sim <- sair_step(stoch = TRUE, Ncomp, ICs, params, time, delta.t)
  #   run_index = rep(n, nrow(single.sim))
  #   single.sim <- cbind(run_index, single.sim)
  #   all_sim <- rbind(all_sim, single.sim)
  #   all_sim[start_index[n]:(start_index[n+1]-1),] = single.sim
  # }
  # write.csv(all_sim, file = 'Output_20200322/TEST_SEIR_results__n500__r02_2__c1.csv')
  # Sys.time()

## ---- loop over r0 and c values ---- ####
  # write_output_files = TRUE
  # r0_values <- c(1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5)
  # c_values <- c(1, 0.75, 0.5, 0.25, 0.1, 0.01)
  # 
  # delta.t <- 1/1
  # time <- seq(1,365,by = delta.t)
  # Ncomp = all_prelim_info$Ncomp
  # ICs = all_prelim_info$ICs
  # 
  # for(ii in 1:length(r0_values)){
  #   for(jj in 1:length(c_values)){
  #     R0_test = r0_values[ii]
  #     c_test = c_values[jj]
  #     print(c(R0_test, c_test))
  #     nsim <- 500
  #     start_index <- seq(1, nsim*length(time)+1, by = length(time))
  #     all_sim <- matrix(NA,1,(Ncomp*7)+2)
  #     colnames(all_sim) <- colnames(test_sim)#matrix(,nsim*length(time),(Ncomp*5)+2) ## +2 -> time step, run index
  #     all_prelim_info <- setup_seir_model(stoch = TRUE, R0 = R0_test, c_scale_vec = 1)
  #     params = list(C = all_prelim_info$C, 
  #                   W = all_prelim_info$W, 
  #                   beta0 = all_prelim_info$beta0, 
  #                   beta1 = all_prelim_info$beta1, 
  #                   phase = all_prelim_info$phase, 
  #                   mu = all_prelim_info$mu, 
  #                   v = all_prelim_info$v, 
  #                   N=all_prelim_info$N, 
  #                   gamma=all_prelim_info$gamma,
  #                   prop_symptomatic=all_prelim_info$prop_symptomatic, 
  #                   sigma = all_prelim_info$sigma)
  #     for(n in 1:nsim){
  #       #  prop_serious <- 0.2 ### will update later
  #       ## run the simulation one time
  #       single.sim <- sair_step(stoch = TRUE, Ncomp, ICs, params, time, delta.t)
  #       run_index = rep(n, nrow(single.sim))
  #       single.sim <- cbind(run_index, single.sim)
  #       all_sim <- rbind(all_sim, single.sim)
  #       all_sim[start_index[n]:(start_index[n+1]-1),] = single.sim
  #     }
  #     if(write_output_files == TRUE){
  #       write.csv(all_sim, file = paste(paste(paste('Output_20200322_v2/SEIR_results__n500__r0', R0_test*10, sep = ''), c_test*100, sep = '__'), 'csv', sep = '.'))
  #     }
  # }
  # }
  # 
  # 
  # library(tidyverse)
  # ## c(0.9797704, 9.487135e-05, 5.578693e-05, 8.364770e-05, 0.01999531)
  # 
  # write_summary_files = TRUE
  # if(write_summary_files == TRUE){
  #   ### loop over r0 and c values to write summary incidence files will produce 
  #   ### new files that have the mean, median, IQR and 95% quantile per time step across each incidence class. 
  #   r0_values <- c(1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5)
  #   c_values <- c(1, 0.75, 0.5, 0.25, 0.1, 0.01)
  #   first_incid_name = 'incid_A1'
  #   last_incid_name = 'incid_I14'
  #   n_incid = 14
  #   ## mean, median, IQR, 2.5%, 97.5%
  #   for(ii in 1:length(r0_values)){
  #     for(jj in 1:length(c_values)){
  #       R0_test = r0_values[ii]
  #       c_test = c_values[jj]
  #       print(c(R0_test, c_test))
  #       file_name = paste(paste(paste('Output_20200322_v2/SEIR_results__n500__r0', R0_test*10, sep = ''), c_test*100, sep = '__'), 'csv', sep = '.')
  #       data_file <- read.csv(file_name)
  #       inc_data <- data_file[,c(grep('time', colnames(data_file)), min(grep(first_incid_name, colnames(data_file))):max(grep(last_incid_name, colnames(data_file))))]
  #       inc_data[is.na(inc_data)] <- 0
  #       summary_inc_data <- inc_data %>% 
  #                           group_by(time) %>% 
  #                           summarise_all(.funs = list(mean = mean, median = median, IQR = IQR, Q1 =~quantile(x=.,probs = 0.025), Q4 = ~quantile(x=., probs = 0.975)))
  #       output_file_name <- paste(paste(paste('Output_20200322_v2/SUMMARY_INCIDENCE_SEIR_results__n500__r0', R0_test*10, sep = ''), c_test*100, sep = '__'), 'csv', sep = '.')
  #       write.csv(summary_inc_data, output_file_name, row.names = FALSE)
  #     }
  #   }
  # }

# 
  #   daily_incid <- unname(tapply(totalincid, (seq_along(totalincid)-1) %/% (1/delta.t), sum))
# ## will need to be changed



# 
# 
# 
# 
#   single.sim %>%
#     ggplot(aes(time,I5))+geom_line()
# 
#   single.sim %>%
#     dplyr::select(paste0('I',1:Ncomp)) %>%
#     rowSums() -> totalI
# 
#   single.sim %>%
#     dplyr::select(paste0('incid',1:Ncomp)) %>%
#     rowSums() -> totalincid
# 
#   daily_incid <- unname(tapply(totalincid, (seq_along(totalincid)-1) %/% (1/delta.t), sum))
# ## will need to be changed
#   cases_requiring_attention <- totalI * prop_serious
# 
#   daily_cases_requiring_attention <- unname(tapply(cases_requiring_attention, (seq_along(cases_requiring_attention)-1) %/% (1/delta.t), sum))
# 
#   dailytime <- single.sim$time[seq(1,nrow(single.sim),(1/delta.t))]
# 
#   single.daily.data <- data.frame('days'=dailytime,'daily_cases_requiring_attention'=daily_cases_requiring_attention)
#   single.daily.data$sim <- n
# 
#   if(n == 1){
#     daily.data <- single.daily.data
#   }else{
#     daily.data <- rbind(daily.data,single.daily.data)
#   }
# }
# 
# daily.data %>%
#   ggplot(aes(days,daily_cases_requiring_attention,color=factor(sim)))+
#   geom_line()+
#   xlim(0,max(time) - 1)


