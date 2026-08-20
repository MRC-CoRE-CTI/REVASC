# MAIN POWER SIMULATION
# simrun.R
# IW 31mar2026 building on Wenyue's code

library(Rlab)
library(survival)
# library(SurvRegCensCov)
library(foreach)
library(doParallel)
library(purrr)
library(bayesmeta)

setwd("C:/temp/Revasc")
source("gendata.R")
source("prepdata.R")
source("anadata.R")
source("manyreps.R")

# common DGM parameters
target_hr = 0.7
crossover = 0.08
prevalence_set = c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
acc_period = 36
time_censor = 48
K = length(prevalence_set)
common_shape = 0.48 # CORRECTED 31MAR2026
cont_cuminc_5y = c(.623,.490,.398,.609,.392,.720,.689) # Weiqi's 5-year cum inc for PCI group
scale_con = 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape) # Weibull calculation


simrun <- function(nonnull,N_tot,nSims){
  results <- manyreps(
    # specific DGM parameters
    N_tot = N_tot, nonnull=nonnull,
    # common DGM parameters
    target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
    # analysis parameters
    sep = TRUE, scale_hn = 0.26, 
    # sim parameters
    nSims = nSims, seed = 100
  )
  alpha <- 0.025 # 1-sided nominal sig level
  powers <- colMeans(results>1-alpha) # cluster-wise T1ER/power
  fwer <- mean((results>1-alpha)%*%(1-nonnull)>0) # FWER
  avepower <- sum(nonnull*powers)/sum(nonnull) # average power
  wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
  ressep <- cbind(t(powers),fwer,avepower,wtavepower)

  results <- manyreps(
    # specific DGM parameters
    N_tot = N_tot, nonnull=nonnull,
    # common DGM parameters
    target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
    # analysis parameters
    sep = FALSE, scale_hn = 0.26,
    # sim parameters
    nSims = nSims, seed = 100
  )
  alpha <- 0.025 # 1-sided nominal sig level
  powers <- colMeans(results>1-alpha) # cluster-wise T1ER/power
  fwer <- mean((results>1-alpha)%*%(1-nonnull)>0) # FWER
  avepower <- sum(nonnull*powers)/sum(nonnull) # average power
  wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
  resborrow <- cbind(t(powers),fwer,avepower,wtavepower)
  
  res <- rbind(ressep, resborrow)
  rownames(res)<-c("separate","borrow")
  return(res)
}

nSims=5000
### GLOBAL NULL ###
simrun(nonnull=rep(0,K),N_tot=1500,nSims=nSims)

### GLOBAL ALTERNATIVE ###
simrun(nonnull=rep(1,K),N_tot=1500,nSims=nSims)
simrun(nonnull=rep(1,K),N_tot=2000,nSims=nSims)
simrun(nonnull=rep(1,K),N_tot=2500,nSims=nSims)
simrun(nonnull=rep(1,K),N_tot=3000,nSims=nSims)

N_tot=3000
### PARTIAL ALTERNATIVE 1: ONLY SMALLEST CLUSTER (2) IS NULL
simrun(nonnull=c(1,0,1,1,1,1,1),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 2: ONLY LARGEST CLUSTER (4) IS NULL
simrun(nonnull=c(1,1,1,0,1,1,1),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 3: ONLY SMALLEST CLUSTER (2) IS NON-NULL
simrun(nonnull=c(0,1,0,0,0,0,0),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 4: ONLY LARGEST CLUSTER (4) IS NON-NULL
simrun(nonnull=c(0,0,0,1,0,0,0),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 5: THREE LARGEST CLUSTERS (5,1,4) ARE NULL
simrun(nonnull=c(0,1,1,0,0,1,1),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 6: THREE LARGEST CLUSTERS (5,1,4) ARE NON-NULL
simrun(nonnull=c(1,0,0,1,1,0,0),N_tot=N_tot,nSims=nSims)


### SAVE AND END ###

save.image(file="simrun.Rdata")
# took about 12h