# IW anadata.R built from Wenyue's ma_bayesmeta.R
# takes subgroup-specific logHRs and SEs and returns posterior probs of theta_k<0

anadata <- function(prepdata, scale_hn=0.26) {
	beta.hat <- prepdata[[1]]
	sd.within <- prepdata[[2]]
	K <- length(beta.hat)

  # Second-stage analysis
  mybayesmeta <- bayesmeta(y = beta.hat,
                        sigma = sd.within,
                        tau.prior = function(t){dhalfnormal(t, scale=scale_hn)},
                        mu.prior = c("mean"=0,"sd"=sqrt(1000)),
                        interval.type = "central",
                        delta = 0.01, epsilon = 0.0001)

	# delta, epsilon are tuning parameters for the bayesmeta computational method

  # Subgroup-specific HRs
  # mean_beta[1, ] <- as.numeric(results$theta['mean', ])
  
  # Average treatment effect
  # average_beta[1] <- results$summary['mean', 'mu']
  
  # Between-subgroup heterogeneity
  # between_sd[1] <- results$summary['mean', 'tau']
  
  # Power/Type I error
  # this is p(theta1<0|data) - must be >1-alpha/2 to reject: 
  pprob <- vapply(1:K,
    function(id) mybayesmeta$pposterior(theta = 0, individual = id),
    numeric(1)
  )
  # to get posterior interval for study 1: 
      # bayesmeta$post.interval(mu.level=.95,individual=1)

  return(pprob)
}
