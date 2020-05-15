

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
  t_March13 <- as.numeric(as.Date("2020-03-13") - as.Date("2020-03-01"))
  t_March30 <- as.numeric(as.Date("2020-03-30") - as.Date("2020-03-01"))
  t_April19 <- as.numeric(as.Date("2020-04-19") - as.Date("2020-03-01"))
  t_April26 <- as.numeric(as.Date("2020-04-26") - as.Date("2020-03-01"))
  t_May15 <- as.numeric(as.Date("2020-05-15") - as.Date("2020-03-01"))
  t_May30 <- as.numeric(as.Date("2020-05-30") - as.Date("2020-03-01"))
  t_June15 <- as.numeric(as.Date("2020-06-15") - as.Date("2020-03-01"))
  t_July1 <- as.numeric(as.Date("2020-07-01") - as.Date("2020-03-01"))
  t_July31 <- as.numeric(as.Date("2020-07-31") - as.Date("2020-03-01"))
  t_Aug1<- as.numeric(as.Date("2020-08-01") - as.Date("2020-03-01"))
  t_Aug15<- as.numeric(as.Date("2020-08-15") - as.Date("2020-03-01"))
  t_Oct1<- as.numeric(as.Date("2020-10-01") - as.Date("2020-03-01"))
  
  Rt_stay_at_home_phase3 <- c(0.8, 1.08) #Baltimore average April 27 - May 4
  Rt_stay_at_home_phase2 <- c(1.016, 1.107) #Baltimore average April 20 - April 26
  Rt_stay_at_home_phase1 <- c(1.016, 1.288) #Baltimore average March 30 - April 19
  Rt_safer_at_home <- c(1.148, 1.198) # georgia average May 1 - 10
  Rt_mod_sd <- c(1.337, 1.990) # Baltimore average March 20 - 29
  Rt_unc <- c(2, 3)
  
  nsim <- 1000
  start_index <- seq(1, nsim*length(time)+1, by = length(time))
  
  all_prelim_info <- setup_seir_model(stoch = TRUE, 
                                      R0 = 2, 
                                      c_scale_vec = 1, 
                                      prop_symptomatic = prop_symptomatic,
                                      sd.dw = 0.1, 
                                      healthcare_n = 26890)
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
                     "R0",
                     "Reff")
  
## ---- Scenario 1 ---- ####
#  stay at home March 29 - June 15
#  safer at home June 16 - Aug 15
#  moderate social distancing Aug 16 - Dec 31
  
  all_sim <- matrix(NA,1,(Ncomp*7)+4)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- rep(runif(1, min = 2, max = 3), length(time))
    R0vec[t_March13:t_March30] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    R0vec[(t_March30+1):t_April19] <- runif(1, min = Rt_stay_at_home_phase1[1], max = Rt_stay_at_home_phase1[2])
    R0vec[(t_April19+1):t_April26] <- runif(1, min = Rt_stay_at_home_phase2[1], max = Rt_stay_at_home_phase2[2])
    R0vec[(t_April26+1):t_June15] <- runif(1, min = Rt_stay_at_home_phase3[1], max = Rt_stay_at_home_phase3[2])
    R0vec[(t_June15+1):t_Aug15] <- runif(1, min = Rt_safer_at_home[1], max = Rt_safer_at_home[2])
    R0vec[(t_Aug15+1):length(time)] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    
    c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
    
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec = R0vec, Ncomp = Ncomp, 
                                ICs = ICs, params = params, 
                                time = time, delta.t = delta.t, 
                                c_scale_mat = c_scale_mat)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200514/scenario1.csv")
  
  
## ---- Scenario 2 --- ####
# stay at home March 29 - May 30
# safer at home June 1 - July 31
# moderate social distancing Aug 1 - Dec 31

  all_sim <- matrix(NA,1,(Ncomp*7)+4)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- rep(runif(1, min = 2, max = 3), length(time))
    R0vec[t_March13:t_March30] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    R0vec[(t_March30+1):t_April19] <- runif(1, min = Rt_stay_at_home_phase1[1], max = Rt_stay_at_home_phase1[2])
    R0vec[(t_April19+1):t_April26] <- runif(1, min = Rt_stay_at_home_phase2[1], max = Rt_stay_at_home_phase2[2])
    R0vec[(t_April26+1):t_May30] <- runif(1, min = Rt_stay_at_home_phase3[1], max = Rt_stay_at_home_phase3[2])
    R0vec[(t_May30+1):t_July31] <- runif(1, min = Rt_safer_at_home[1], max = Rt_safer_at_home[2])
    R0vec[(t_Aug1+1):length(time)] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    
    c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
    
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec = R0vec, Ncomp = Ncomp, 
                                ICs = ICs, params = params, 
                                time = time, delta.t = delta.t, 
                                c_scale_mat = c_scale_mat)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200514/scenario2.csv")
  
## ---- Scenario 3 ---- ####
#  stay at home March 29 - May 30
#  safer at home June 1 - June 30
#  moderate social distancing July 1 - Sept 30
#  uncontrolled Oct 1 - Dec 31

  all_sim <- matrix(NA,1,(Ncomp*7)+4)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- rep(runif(1, min = 2, max = 3), length(time))
    R0vec[t_March13:t_March30] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    R0vec[(t_March30+1):t_April19] <- runif(1, min = Rt_stay_at_home_phase1[1], max = Rt_stay_at_home_phase1[2])
    R0vec[(t_April19+1):t_April26] <- runif(1, min = Rt_stay_at_home_phase2[1], max = Rt_stay_at_home_phase2[2])
    R0vec[(t_April26+1):t_May30] <- runif(1, min = Rt_stay_at_home_phase3[1], max = Rt_stay_at_home_phase3[2])
    R0vec[(t_May30+1):(t_July1-1)] <- runif(1, min = Rt_safer_at_home[1], max = Rt_safer_at_home[2])
    R0vec[t_July1:(t_Oct1-1)] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    
    c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)
    
    
    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec = R0vec, Ncomp = Ncomp, 
                                ICs = ICs, params = params, 
                                time = time, delta.t = delta.t, 
                                c_scale_mat = c_scale_mat)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200514/scenario3.csv")
  
  
## ---- Scenario 4 --- ####
#  stay at home March 29 - May 15
#  safer at home May 16 - June 15
#  moderate social distancing June 16 - Aug 15
#  uncontrolled Aug 16 - Dec 31

  all_sim <- matrix(NA,1,(Ncomp*7)+4)
  colnames(all_sim) <- cnames.allsim
  
  for(n in 1:nsim){
    R0vec <- rep(runif(1, min = 2, max = 3), length(time))
    R0vec[t_March13:t_March30] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    R0vec[(t_March30+1):t_April19] <- runif(1, min = Rt_stay_at_home_phase1[1], max = Rt_stay_at_home_phase1[2])
    R0vec[(t_April19+1):t_April26] <- runif(1, min = Rt_stay_at_home_phase2[1], max = Rt_stay_at_home_phase2[2])
    R0vec[(t_April26+1):t_May15] <- runif(1, min = Rt_stay_at_home_phase3[1], max = Rt_stay_at_home_phase3[2])
    R0vec[(t_May15+1):t_June15] <- runif(1, min = Rt_safer_at_home[1], max = Rt_safer_at_home[2])
    R0vec[(t_June15+1):t_Aug15] <- runif(1, min = Rt_mod_sd[1], max = Rt_mod_sd[2])
    
    c_scale_mat <- matrix(1, nrow = length(time), ncol=Ncomp)

    tmp <- sair_step_variableR0(stoch = TRUE, stoch.init = TRUE, 
                                R0vec = R0vec, Ncomp = Ncomp, 
                                ICs = ICs, params = params, 
                                time = time, delta.t = delta.t, 
                                c_scale_mat = c_scale_mat)
    run_index = rep(n, nrow(tmp))
    tmp <- cbind(run_index, tmp)
    all_sim <- rbind(all_sim, tmp)
  }
  all_sim <- all_sim[-1,]
  write.csv(all_sim, file="output_20200514/scenario4.csv")
  
  
