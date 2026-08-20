# MAIN POWER SIMULATION
# simrun.R
# IW 31mar2026 building on Wenyue's code
# IW 20aug2026 tidied up

source("setup.R")

simrun <- function(nonnull,N_tot,nSims){
  for (sep in c(TRUE, FALSE)) {
  results <- manyreps(
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
    sep = sep, 
    scale_hn = 0.26, 
    # sim parameters
    nSims = nSims, 
    seed = 100
  )
  alpha <- 0.025 # 1-sided nominal sig level
  powers <- colMeans(results>1-alpha) # subgroup-wise T1ER/power
  fwer <- mean((results>1-alpha)%*%(1-nonnull)>0) # FWER
  avepower <- sum(nonnull*powers)/sum(nonnull) # average power
  wtavepower <- sum(prevalence_set*nonnull*powers)/sum(nonnull*prevalence_set) # weighted average power
  if (sep) ressep <- cbind(t(powers),fwer,avepower,wtavepower)
  else  resborrow <- cbind(t(powers),fwer,avepower,wtavepower)
  }
  res <- rbind(c(nonnull,rep(NA,3)), ressep, resborrow)
  rownames(res) <- c("nonnull","separate","borrow")
  res <- cbind(N_tot, res)
  colnames(res2)[1] <- "N_tot"
  return(res)
}

nSims=500

### GLOBAL NULL ###
resgn <- simrun(nonnull=rep(0,K),N_tot=1500,nSims=nSims)

### GLOBAL ALTERNATIVE, VARY N ###
resga1 <- simrun(nonnull=rep(1,K),N_tot=1500,nSims=nSims)
resga2 <- simrun(nonnull=rep(1,K),N_tot=2000,nSims=nSims)
resga3 <- simrun(nonnull=rep(1,K),N_tot=2500,nSims=nSims)
resga4 <- simrun(nonnull=rep(1,K),N_tot=3000,nSims=nSims)

N_tot=3000
### PARTIAL ALTERNATIVE 1: ONLY SMALLEST subgroup (2) IS NULL
respa1 <- simrun(nonnull=c(1,0,1,1,1,1,1),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 2: ONLY LARGEST subgroup (4) IS NULL
respa2 <- simrun(nonnull=c(1,1,1,0,1,1,1),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 3: ONLY SMALLEST subgroup (2) IS NON-NULL
respa3 <- simrun(nonnull=c(0,1,0,0,0,0,0),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 4: ONLY LARGEST subgroup (4) IS NON-NULL
respa4 <- simrun(nonnull=c(0,0,0,1,0,0,0),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 5: THREE LARGEST subgroupS (5,1,4) ARE NULL
respa5 <- simrun(nonnull=c(0,1,1,0,0,1,1),N_tot=N_tot,nSims=nSims)

### PARTIAL ALTERNATIVE 6: THREE LARGEST subgroupS (5,1,4) ARE NON-NULL
respa6 <- simrun(nonnull=c(1,0,0,1,1,0,0),N_tot=N_tot,nSims=nSims)

### SAVE AND END ###

rbind(resgn, resga1, resga2, resga3, respa1, respa2, respa3, respa4, respa5, respa6)

save.image(file="simrun.Rdata")
# 5000 reps took about 12h