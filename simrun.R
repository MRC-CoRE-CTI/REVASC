# MAIN POWER SIMULATION
# simrun.R
# IW 31mar2026 building on Wenyue's code
# IW 20aug2026 tidied up
# IW 24aug2026 manyreps returns both sep and BB

# clear workspace
rm(list = ls())
options(width = 120)

# load all simulation programs
source("setup.R")

# DGM settings for this run
target_hr <- 0.7
crossover <- 0.08
common_shape <- 1.1
acc_period <- 36
time_censor <- 48
run <- 2 # 1: 7 classes incl. HF; 2: 6 classes excl. HF
if(run==1) {
  # parameters with HF patients (7 subgroups)
  prevalence_set <- c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
  cont_cuminc_5y <- c(.623, .490, .398, .609, .392, .720, .689) # Weiqi's 5-year cum inc for PCI group
}
if(run==2){
  # parameters without HF patients (6 subgroups)
  # figures from Weiqi 24/8/2026
  nobs <- c(21878, 6477, 24506, 23179, 6153, 4323)
  nobspci <- c(8002, 2138, 9336, 16331, 1523, 875)
  neventpci <- c(4857, 850, 5605, 6350, 1017, 591) # all-cause mortality or CV hospitalisation to 5 years
  prevalence_set <- nobs / sum(nobs)
  cont_cuminc_5y <- neventpci / nobspci
}

# derived parameters
K <- length(prevalence_set)
scale_con <- 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape) # Weibull calculation

# simulation settings
N_reps <- 5000
# I'll set the seed before each run, so that the script can be run in parts

# define simulation program
sim_run <- function(nonnull, N_tot, N_reps) {
  K <- length(prevalence_set)
  alpha <- 0.025 # 1-sided nominal sig level
  # do simulation
  resultsall <- manyreps(
    N_tot = N_tot,
    target_hr = target_hr,
    crossover = crossover,
    nonnull = nonnull,
    prevalence_set = prevalence_set,
    common_shape = common_shape,
    scale_con = scale_con,
    acc_period = acc_period,
    time_censor = time_censor,
    # analysis parameters
    scale_hn = 0.26,
    # sim parameters
    N_reps = N_reps,
    seed = 100
  )
  # compute T1ERs and powers for separate and BB analyses
  for(i in 1:2) {
    if(i==1) results <- resultsall[,1:K]
    else results <- resultsall[,(K+1):(2*K)]
    powers <- colMeans(results > 1 - alpha) # subgroup-wise T1ER/power
    fwer <- mean((results > 1 - alpha) %*% (1 - nonnull) > 0) # FWER
    avepower <- sum(nonnull * powers) / sum(nonnull) # average power
    wtavepower <- sum(prevalence_set * nonnull * powers) /
      sum(nonnull * prevalence_set) # weighted average power
    if(i==1) ressep <- cbind(t(powers), fwer, avepower, wtavepower)
    else  resborrow <- cbind(t(powers), fwer, avepower, wtavepower)
  }
  # tidy up and output
  res <- rbind(c(nonnull, rep(NA, 3)), ressep, resborrow)
  rownames(res) <- c("nonnull", "separate", "borrow")
  res <- cbind(N_tot, N_reps, res)
  colnames(res)[1] <- "N_tot"
  colnames(res)[2] <- "N_reps"
  colnames(res)[3:(K+2)] <- c(paste0("Cl", 1:K))
  return(res)
}
# now run the simulations

### GLOBAL NULL ###
set.seed(101)
resgn <- sim_run(nonnull = rep(0, K), N_tot = 1500, N_reps = N_reps)
resgn

### GLOBAL ALTERNATIVE, VARY N ###
for(i in 1:4) {
  set.seed(101)
  N_tot=1000+500*i # 1500 to 3000
  print("Global alternative, N_tot =", N_tot)
  thisresult = paste0("resga",i)
  assign(thisresult, sim_run(nonnull = rep(1, K), N_tot = N_tot, N_reps = N_reps))
  print(get(thisresult))
}

### PARTIAL ALTERNATIVES, VARY NONNULL ###
if(run==1) {
  N_tot <- 3000
  ### PARTIAL ALTERNATIVE 1: ONLY SMALLEST SUBGROUP (2) IS NULL
  nonnull_pa1 <- c(1, 0, 1, 1, 1, 1, 1)
  ### PARTIAL ALTERNATIVE 2: ONLY LARGEST SUBGROUP (4) IS NULL
  nonnull_pa2 <- c(1, 1, 1, 0, 1, 1, 1)
  ### PARTIAL ALTERNATIVE 3: ONLY SMALLEST SUBGROUP (2) IS NON-NULL
  nonnull_pa3 <- c(0, 1, 0, 0, 0, 0, 0)
  ### PARTIAL ALTERNATIVE 4: ONLY LARGEST SUBGROUP (4) IS NON-NULL
  nonnull_pa4 <- c(0, 0, 0, 1, 0, 0, 0)
  ### PARTIAL ALTERNATIVE 5: THREE LARGEST SUBGROUPS (5,1,4) ARE NULL
  nonnull_pa5 <- c(0, 1, 1, 0, 0, 1, 1)
  ### PARTIAL ALTERNATIVE 6: THREE LARGEST SUBGROUPS (5,1,4) ARE NON-NULL
  nonnull_pa6 <- c(1, 0, 0, 1, 1, 0, 0)
}
if(run==2) {
  N_tot <- 2500
  ### PARTIAL ALTERNATIVE 1: ONLY SMALLEST SUBGROUP (6) IS NULL
  nonnull_pa1 <- c(1, 1, 1, 1, 1, 0)
  ### PARTIAL ALTERNATIVE 2: ONLY LARGEST SUBGROUP (3) IS NULL
  nonnull_pa2 <- c(1, 1, 0, 1, 1, 1)
  ### PARTIAL ALTERNATIVE 3: ONLY SMALLEST SUBGROUP (6) IS NON-NULL
  nonnull_pa3 <- c(0, 0, 0, 0, 0, 1)
  ### PARTIAL ALTERNATIVE 4: ONLY LARGEST SUBGROUP (3) IS NON-NULL
  nonnull_pa4 <- c(0, 0, 1, 0, 0, 0)
  ### PARTIAL ALTERNATIVE 5: THREE LARGEST SUBGROUPS (1,3,4) ARE NULL
  nonnull_pa5 <- c(0, 1, 0, 0, 1, 1)
  ### PARTIAL ALTERNATIVE 6: THREE LARGEST SUBGROUPS (1,3,4) ARE NON-NULL
  nonnull_pa6 <- c(1, 0, 1, 1, 0, 0)
}
for(i in 1:6){
  set.seed(101)
  print("Partial alternative", i)
  thisresult <- paste0("respa",i) 
  thisnonnull <- paste0("nonnull_pa",i)
  assign(thisresult, sim_run(get(thisnonnull), N_tot = N_tot, N_reps = N_reps))
  print(get(thisresult))
}

### SAVE AND END ###
allresults <- rbind(resgn, resga1, resga2, resga3, resga4, respa1, respa2, respa3, respa4, respa5, respa6)
allresults
library(writexl)
outfile <- paste0("simrun_results_run", run)
write_xlsx(as.data.frame(allresults), paste0(outfile,".xlsx"))
save.image(file = paste0(outfile,".Rdata"))
# 5000 reps took about 12h