

## Setup ####

  require(socialmixr)
  require(magrittr)
  require(stringr)
  require(reshape2)
  require(dplyr)
  require(ggplot2)
  require(truncnorm)

  source("asyptomatic_age.R")

  age.limits <- c(0,5,10,15,20,25,35,45,55,60,65,75,85,90)
  prop_symptomatic <- c(0.141, 0.106, 0.074, 0.184, 0.293, 0.387, 0.438, 
                        0.535, 0.693, 0.816, 0.765, 0.749, 0.535, 0.535)
  
  delta.t <- 1/1
  time <- seq(1,300,by = delta.t)
  t_March19 <- as.numeric(as.Date("2020-03-19") - as.Date("2020-03-01"))
  t_May14 <- as.numeric(as.Date("2020-05-14") - as.Date("2020-03-01"))
  
  nsim <- 500
  start_index <- seq(1, nsim*length(time)+1, by = length(time))
  
  all_prelim_info <- setup_seir_model(stoch = TRUE, 
                                      R0 = 2, 
                                      c_scale_vec = 1, 
                                      prop_symptomatic = prop_symptomatic,
                                      sd.dw = 0.1)
  Ncomp = all_prelim_info$Ncomp
  ICs = all_prelim_info$ICs
  params = list(C = all_prelim_info$C, 
                W = all_prelim_info$W, 
                beta0 = all_prelim_info$beta0, 
                beta1 = all_prelim_info$beta1, 
                phase = all_prelim_info$phase, 
                mu = all_prelim_info$mu, 
                v = all_prelim_info$v, 
                N=all_prelim_info$N, 
                gamma=all_prelim_info$gamma, 
                sigma = all_prelim_info$sigma, 
                prop_symptomatic=all_prelim_info$prop_symptomatic,
                sd.dw = all_prelim_info$sd.dw)
  
  
  cnames.allsim <- c('run_index', 'time',
                     paste0("S", 1:Ncomp),
                     paste0("E", 1:Ncomp),
                     paste0("A", 1:Ncomp),
                     paste0("I", 1:Ncomp),
                     paste0("R", 1:Ncomp),
                     paste0("incid_A", 1:Ncomp),
                     paste0("incid_I", 1:Ncomp),
                     "R0")
  
## Scenario 1A: R0 2 - 3, uncontrolled ####
  all_sim <- matrix(NA,1,(Ncomp*7)+3)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- runif(length(time), min = 2, max = 3)
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec, Ncomp, ICs, params, time, delta.t)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200408/scenario1A_lowR0_uncontrolled.csv")
  
  
## Scenario 1B: R0 2 - 3, 30-40% reduction March 19 - May 14 ####
  all_sim <- matrix(NA,1,(Ncomp*7)+3)
  colnames(all_sim) <- cnames.allsim
    
  for(n in 1:nsim){
    R0vec <- runif(length(time), min = 2, max = 3)
    R0vec[t_March19:t_May14] <- runif(length(t_March19:t_May14), 2*(1-0.4), 3*(1-0.3))
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec, Ncomp, ICs, params, time, delta.t)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200408/scenario1B_lowR0_MildDistancing.csv")
  
  
## Scenario 1C: R0 2 - 3, 45-65% reduction March 19 - May 14, 48-76% May 15 - Dec 26 ####
  all_sim <- matrix(NA,1,(Ncomp*7)+3)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- runif(length(time), min = 2, max = 3)
    R0vec[t_March19:t_May14] <- runif(length(t_March19:t_May14), 2*(1-0.65), 3*(1-0.45))
    R0vec[(t_May14+1):length(time)] <- runif(length((t_May14+1):length(time)), 2*(1-0.76), 3*(1-0.48))
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec, Ncomp, ICs, params, time, delta.t)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200408/scenario1C_lowR0_ModDistancing.csv")

  
## Scenario 2A: R0 3.5-4, uncontrolled ####
  all_sim <- matrix(NA,1,(Ncomp*7)+3)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- runif(length(time), min = 3.5, max = 4)
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec, Ncomp, ICs, params, time, delta.t)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200408/scenario2A_highR0_uncontrolled.csv")
  
  
## Scenario 2B: R0 3.5-4, 30-40% reduction March 19 - May 14 ####
  all_sim <- matrix(NA,1,(Ncomp*7)+3)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- runif(length(time), min = 3.5, max = 4)
    R0vec[t_March19:t_May14] <- runif(length(t_March19:t_May14), 3.5*(1-0.4), 4*(1-0.3))
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec, Ncomp, ICs, params, time, delta.t)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200408/scenario2B_highR0_MildDistancing.csv")
  
  
## Scenario 2C: R0 3.5 - 4, 45-65% reduction March 19 - May 14, 48-76% May 15 - Dec 26 ####
  all_sim <- matrix(NA,1,(Ncomp*7)+3)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- runif(length(time), min = 3.5, max = 4)
    R0vec[t_March19:t_May14] <- runif(length(t_March19:t_May14), 3.5*(1-0.65), 4*(1-0.45))
    R0vec[(t_May14+1):length(time)] <- runif(length((t_May14+1):length(time)), 3.5*(1-0.76), 4*(1-0.48))
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec, Ncomp, ICs, params, time, delta.t)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200408/scenario2C_highR0_ModDistancing.csv")
  
  
  
  

  
