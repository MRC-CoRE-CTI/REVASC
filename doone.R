N_tot = 1500
nSims=5
target_hr = 0.7
crossover = 0.08
nonnull = rep(1,K)
prevalence_set = c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
common_shape = 1.1
cont_cuminc_5y = c(.623,.490,.398,.609,.392,.720,.689) # Weiqi's numbers
scale_con = 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape)
# scale_con = c(50, 68, 86, 50, 90, 38, 42) # Wenyue's numbers
sep=TRUE
scale_hn = 0.26
acc_period = 36
time_censor = 48
seed = 100



set.seed(seed)

K <- length(prevalence_set)



mygendata <- gendata(N_tot = N_tot,
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
}else {
  myanadata <- anadata(myprepdata,
                       scale_hn = scale_hn)
}

myanadata
