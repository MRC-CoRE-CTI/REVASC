# IW anadata.R built from Wenyue's ma_bayesmeta.R
# takes cluster-specific logHRs and SEs and returns posterior probs of theta_k<0
# EITHER from bayesmeta OR from separate analyses


anadata <- function(prepdata, scale_hn=0.26)
{
	beta_hat <- prepdata[[1]]    
	sd.within <- prepdata[[2]]    
	K <- length(beta_hat)

		  # Second-stage analysis
      mybayesmeta <- bayesmeta(y = beta_hat,
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
      pprob <- c(mybayesmeta$pposterior(theta = 0, individual = 1),
                       mybayesmeta$pposterior(theta = 0, individual = 2),
                       mybayesmeta$pposterior(theta = 0, individual = 3),
                       mybayesmeta$pposterior(theta = 0, individual = 4),
                       mybayesmeta$pposterior(theta = 0, individual = 5),
                       mybayesmeta$pposterior(theta = 0, individual = 6),
                       mybayesmeta$pposterior(theta = 0, individual = 7)
      ) 
	  # bayesmeta$pposterior(theta = 0, individual = 1) is p(theta1<0|data)
	  # must be >1-alpha/2 to reject
      # to get posterior interval for study 1: bayesmeta$post.interval(mu.level=.95,individual=1)
	  
      
      results <- list(
        # adj_N_tot,
        # mean_beta,
        # average_beta,
        # between_sd,
        # coverage,
        pprob
        # var,
        # squared_error,
        # low,
        # upp,
        # wid_ci
      )

    return(pprob)
  }
