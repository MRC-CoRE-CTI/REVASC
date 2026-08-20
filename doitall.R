# MAIN POWER SIMULATION
# doitall.R
# IW 27mar2026 building on Wenyue's code

library(Rlab)
library(survival)
# library(SurvRegCensCov)
library(foreach)
library(doParallel)
library(purrr)
library(bayesmeta)

source("S:/MRCCTU_Methodology/1_Design/Projects/High-risk REVASC trial design/Ian code/gendata.R")
source("S:/MRCCTU_Methodology/1_Design/Projects/High-risk REVASC trial design/Ian code/prepdata.R")
source("S:/MRCCTU_Methodology/1_Design/Projects/High-risk REVASC trial design/Ian code/anadata.R")
source("S:/MRCCTU_Methodology/1_Design/Projects/High-risk REVASC trial design/Ian code/manyreps.R")

# common DGM parameters
target_hr = 0.7
crossover = 0.08
prevalence_set = c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
acc_period = 36
time_censor = 48
K = length(prevalence_set)
common_shape = 1.1
cont_cuminc_5y = c(.623,.490,.398,.609,.392,.720,.689) # Weiqi's 5-year cum inc for PCI group
scale_con = 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape) # Weibull calculation
# scale_con = c(50, 68, 86, 50, 90, 38, 42) # Wenyue's numbers



### GLOBAL NULL, SEPARATE ANALYSES ###

start.time <- Sys.time()
nullsep <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=rep(0,K),
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = TRUE, scale_hn = 0.26, 
  # sim parameters
  nSims = 5000, seed = 100
)
T_manyreps_sep <- Sys.time() - start.time
T_manyreps_sep

alpha<-0.025 # 1-sided nominal sig level

colMeans(nullsep>1-alpha) # cluster-wise T1ER/power
mean(rowSums(nullsep>1-alpha)>0) # family-wise T1ER/power



### GLOBAL NULL, BORROWING ANALYSES ###
start.time <- Sys.time()
nullborrow <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=rep(0,K),
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = FALSE, scale_hn = 0.26, 
  # sim parameters
  nSims = 500, seed = 100
)
T_manyreps_borrow <- Sys.time() - start.time
T_manyreps_borrow

alpha<-0.025 # 1-sided nominal sig level

colMeans(nullborrow>1-alpha) # cluster-wise T1ER
mean(rowSums(nullborrow>1-alpha)>0) # family-wise T1ER



### GLOBAL ALTERNATIVE, SEPARATE ANALYSIS ###

galtsep <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=rep(1,K),
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = TRUE, scale_hn = 0.26, 
  # sim parameters
  nSims = 500, seed = 100
)

alpha<-0.025 # 1-sided nominal sig level

colMeans(galtsep>1-alpha) # cluster-wise power
mean(galtsep>1-alpha) # average power



### GLOBAL ALTERNATIVE, BORROWING ANALYSIS ###
nonnull = rep(1,K)
galtborrow <- manyreps(
  # specific DGM parameters
  N_tot = 2650, nonnull=nonnull,
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = FALSE, scale_hn = 0.26, 
  # sim parameters
  nSims = 500, seed = 100
)
alpha<-0.025 # 1-sided nominal sig level
powers<-colMeans(galtborrow>1-alpha) # cluster-wise T1ER/power
avepower <- sum(nonnull*powers)/sum(nonnull) # average power
wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
cbind(t(powers),avepower,wtavepower)


### PARTIAL ALTERNATIVES, SEPARATE ANALYSIS ###

### PALT1: ONLY SMALLEST CLUSTER (2) IS NULL
nonnull=c(1,0,1,1,1,1,1)
palt1sep <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=nonnull,
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = TRUE, scale_hn = 0.26, 
  # sim parameters
  nSims = 5000, seed = 100
)

alpha<-0.025 # 1-sided nominal sig level
powers<-colMeans(palt1sep>1-alpha) # cluster-wise T1ER/power
avepower <- sum(nonnull*powers)/sum(nonnull) # average power
wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
cbind(t(powers),avepower,wtavepower)


### PALT2: ONLY LARGEST CLUSTER (4) IS NULL
nonnull=c(1,1,1,0,1,1,1)
palt2sep <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=nonnull,
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = TRUE, scale_hn = 0.26, 
  # sim parameters
  nSims = 5000, seed = 100
)

alpha<-0.025 # 1-sided nominal sig level
powers<-colMeans(palt2sep>1-alpha) # cluster-wise T1ER/power
avepower <- sum(nonnull*powers)/sum(nonnull) # average power
wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
cbind(t(powers),avepower,wtavepower)


### PALT3: THREE LARGEST CLUSTERS (5,1,4) ARE NULL
nonnull=c(0,1,1,0,0,1,1)
palt3sep <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=nonnull,
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = TRUE, scale_hn = 0.26, 
  # sim parameters
  nSims = 5000, seed = 100
)

alpha<-0.025 # 1-sided nominal sig level
powers<-colMeans(palt3sep>1-alpha) # cluster-wise T1ER/power
avepower <- sum(nonnull*powers)/sum(nonnull) # average power
wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
cbind(t(powers),avepower,wtavepower)


### PALT4: THREE LARGEST CLUSTERS (5,1,4) ARE NON-NULL
nonnull=c(1,0,0,1,1,0,0)
palt4sep <- manyreps(
  # specific DGM parameters
  N_tot = 1500, nonnull=nonnull,
  # common DGM parameters
  target_hr = target_hr, crossover = crossover, prevalence_set = prevalence_set, common_shape = common_shape, scale_con = scale_con, acc_period = acc_period, time_censor = time_censor,
  # analysis parameters
  sep = TRUE, scale_hn = 0.26, 
  # sim parameters
  nSims = 5000, seed = 100
)

alpha<-0.025 # 1-sided nominal sig level
powers<-colMeans(palt4sep>1-alpha) # cluster-wise T1ER/power
avepower <- sum(nonnull*powers)/sum(nonnull) # average power
wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
cbind(t(powers),avepower,wtavepower)

### SAVE AND END ###

save.image(file="doitall.Rdata")
