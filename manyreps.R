# IW manyreps.R built from Wenyue's ma_bayesmeta.R
# Runs the simulation many times

manyreps <- function(
				# DGM parameters
							N_tot, # total sample size 
							target_hr, # target HR 
							crossover, # anticipated cross-over 
							nonnull, # which clusters are non-null
							prevalence_set, # prevalence set
							common_shape, # Weibull shape parameter (common for all clusters)
							scale_con, # control arm scale parameter (vary by cluster)
							acc_period, # recruitment period
							time_censor, # time to censoring
				# analysis parameters
							sep = TRUE,
							scale_hn, # scale parameter of the halfnormal prior
				# sim parameters
							nSims, # number of simulations
							seed # for generating the fixed sample size set
)
{  
	set.seed(seed)

	K <- length(prevalence_set)

	# Ensure cluster exists & register it
	cl <- parallel::makeCluster(4)
	doParallel::registerDoParallel(cl)

#	# Export ALL required user-defined functions / objects to workers
#	clusterExport(cl, varlist = c(
#		"N_tot", "target_hr", "crossover", "nonnull",
#		"prevalence_set", "common_shape", "scale_con",
#		"acc_period", "time_censor", "sep", "scale_hn",
#		"nSims", "seed", 
#		"gendata", "prepdata", "anadata"
#	))

	out <- foreach(
    sim = 1:nSims,
    .packages = c('Rlab','survival','foreach','purrr','bayesmeta'),
    .export = c("gendata", "prepdata", "anadata")
	) %dopar% {
		# use %do% for debugging and change to %dopar% once code works

		mygendata <- gendata(
			N_tot = N_tot,
			target_hr = target_hr, 
			crossover = crossover, 
			nonnull = nonnull,
			prevalence_set 	= prevalence_set,
			common_shape = common_shape,
			scale_con =	 scale_con,
			acc_period = acc_period,
			time_censor = time_censor
		)

		myprepdata <- prepdata(mygendata)
		beta_hat<-myprepdata[[1]]
		sd.within<-myprepdata[[2]]

		if(sep) {
			myanadata <- pnorm(-beta_hat/sd.within)
		}
		else {
			myanadata <- anadata(myprepdata,
				scale_hn = scale_hn)
		}

		return(myanadata)
	}

	df <- as.data.frame(do.call(rbind, out))
	colnames(df) <- paste0("Cl", 1:K)
	return(df)
}


# colMeans(df>1-crit) 
	