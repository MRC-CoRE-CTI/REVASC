# IW manyreps.R built from Wenyue's ma_bayesmeta.R
# Runs the simulation many times
# revised 20aug2026
# 24aug2026 returns both sep and BB

manyreps <- function(
	# DGM parameters
		N_tot, # total sample size 
		target_hr, # target HR 
		crossover, # anticipated cross-over 
		nonnull, # which subgroups are non-null
		prevalence_set, # prevalence set
		common_shape, # Weibull shape parameter (common for all subgroups)
		scale_con, # control arm scale parameter (vary by subgroup)
		acc_period, # recruitment period
		time_censor, # time to censoring
	# analysis parameters
		scale_hn, # scale parameter of the halfnormal prior
	# sim parameters
		N_reps, # number of simulations
		seed # for generating the fixed sample size set
)
{
	cat("Running", N_reps, "repetitions with K =", K, "\n")
	set.seed(seed)

	K <- length(prevalence_set)

	# Ensure cluster exists & register it
	cl <- parallel::makeCluster(4)
	doParallel::registerDoParallel(cl)

	out <- foreach(
		sim = 1:N_reps,
		.packages = c('Rlab','survival','foreach','purrr','bayesmeta'),
		.export = c("gendata", "prepdata", "anadata")
	) %dopar% {
		# use %do% for debugging and change to %dopar% once code works
		mygendata <- gendata(
			N_tot = N_tot,
			target_hr = target_hr, 
			crossover = crossover, 
			prevalence_set 	= prevalence_set,
			common_shape = common_shape,
			scale_con =	 scale_con,
			nonnull = nonnull,
			acc_period = acc_period,
			time_censor = time_censor
		)
		myprepdata <- prepdata(mygendata)
		beta_hat <- myprepdata[[1]]
		sd.within <- myprepdata[[2]]

		myanadata.sep <- pnorm(-beta_hat/sd.within)
		myanadata.bb <- anadata(myprepdata, scale_hn = scale_hn)

		return(c(myanadata.sep, myanadata.bb))
	}
	parallel::stopCluster(cl)

	df <- as.data.frame(do.call(rbind, out))
	colnames(df) <- c(paste0("psep", 1:K), paste0("pbb", 1:K))
	return(df)
}
