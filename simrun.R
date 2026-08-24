# MAIN POWER SIMULATION
# simrun.R
# IW 31mar2026 building on Wenyue's code
# IW 20aug2026 tidied up
# IW 24aug2026 manyreps returns both sep and BB

source("setup.R")

# settings for this run
target_hr <- 0.7
crossover <- 0.08
common_shape <- 1.1
acc_period <- 36
time_censor <- 48

# parameters with HF patients (7 subgroups)
# prevalence_set <- c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
# cont_cuminc_5y <- c(.623, .490, .398, .609, .392, .720, .689) # Weiqi's 5-year cum inc for PCI group

# parameters without HF patients (6 subgroups)
# figures from Weiqi 24/8/2026
nobs <- c(21878, 6477, 24506, 23179, 6153, 4323)
nobspci <- c(8002, 2138, 9336, 16331, 1523, 875)
neventpci <- c(4857, 850, 5605, 6350, 1017, 591) # all-cause mortality or CV hospitalisation to 5 years
prevalence_set <- nobs / sum(nobs)
cont_cuminc_5y <- neventpci / nobspci

scale_con <- 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape) # Weibull calculation

set.seed(101)

K <- length(prevalence_set)


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
    if (i==1) ressep <- cbind(t(powers), fwer, avepower, wtavepower)
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

N_reps <- 1000

### GLOBAL NULL ###
resgn <- sim_run(nonnull = rep(0, K), N_tot = 1500, N_reps = N_reps)
resgn

### GLOBAL ALTERNATIVE, VARY N ###
resga1 <- sim_run(nonnull = rep(1, K), N_tot = 1500, N_reps = N_reps)
resga1
resga2 <- sim_run(nonnull = rep(1, K), N_tot = 2000, N_reps = N_reps)
resga2
resga3 <- sim_run(nonnull = rep(1, K), N_tot = 2500, N_reps = N_reps)
resga3
resga4 <- sim_run(nonnull = rep(1, K), N_tot = 3000, N_reps = N_reps)
resga4

N_tot <- 3000
### PARTIAL ALTERNATIVE 1: ONLY SMALLEST SUBGROUP (2) IS NULL
respa1 <- sim_run(nonnull = c(1, 0, 1, 1, 1, 1, 1), N_tot = N_tot, N_reps = N_reps)
respa1
### PARTIAL ALTERNATIVE 2: ONLY LARGEST SUBGROUP (4) IS NULL
respa2 <- sim_run(nonnull = c(1, 1, 1, 0, 1, 1, 1), N_tot = N_tot, N_reps = N_reps)
respa2
### PARTIAL ALTERNATIVE 3: ONLY SMALLEST SUBGROUP (2) IS NON-NULL
respa3 <- sim_run(nonnull = c(0, 1, 0, 0, 0, 0, 0), N_tot = N_tot, N_reps = N_reps)
respa3
### PARTIAL ALTERNATIVE 4: ONLY LARGEST SUBGROUP (4) IS NON-NULL
respa4 <- sim_run(nonnull = c(0, 0, 0, 1, 0, 0, 0), N_tot = N_tot, N_reps = N_reps)
respa4
### PARTIAL ALTERNATIVE 5: THREE LARGEST SUBGROUPS (5,1,4) ARE NULL
respa5 <- sim_run(nonnull = c(0, 1, 1, 0, 0, 1, 1), N_tot = N_tot, N_reps = N_reps)
respa5
### PARTIAL ALTERNATIVE 6: THREE LARGEST SUBGROUPS (5,1,4) ARE NON-NULL
respa6 <- sim_run(nonnull = c(1, 0, 0, 1, 1, 0, 0), N_tot = N_tot, N_reps = N_reps)
respa6
### SAVE AND END ###

allresults <- rbind(resgn, resga1, resga2, resga3, resga4, respa1, respa2, respa3, respa4, respa5, respa6)
library(writexl)
write_xlsx(as.data.frame(allresults), "simrun_results_v2.xlsx")
save.image(file = "simrun.Rdata")
# 5000 reps took about 12h