# check_manyreps.R

source("setup.R")

manyreps(
    N_tot=1500, 
    target_hr=target_hr, 
    crossover=crossover, # anticipated cross-over 
    nonnull=rep(0,7), # which subgroups are non-null
    prevalence_set=prevalence_set, # prevalence set
    common_shape=common_shape, # Weibull shape parameter (common for all subgroups)
    scale_con=scale_con, # control arm scale parameter (vary by subgroup)
    acc_period=acc_period, # recruitment period
    time_censor=time_censor, # time to censoring
# analysis parameters
    sep = TRUE,
    scale_hn=0.26, # scale parameter of the halfnormal prior
# sim parameters
    nSims=50, # number of simulations
    seed=10 # for generating the fixed sample size set
)
