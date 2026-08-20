# check_gendata.R
# create big data sets that we can analyse in Stata to confirm control S(5)



N_tot = 100000
nSims=5
target_hr = 0.7
crossover = 0.08
nonnullclusters = 7
prevalence_set = c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
common_shape = 1.1
acc_period = 36
time_censor = 72

cont_cuminc_5y = c(.623,.490,.398,.609,.392,.720,.689) # Weiqi's numbers
scale_con = 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape)
# scale_con = c(50, 68, 86, 50, 90, 38, 42) # Wenyue's numbers

set.seed(101)

K <- length(prevalence_set)



mygendata <- gendata(N_tot = N_tot,
                     target_hr = target_hr, 
                     crossover = crossover, 
                     nonnullclusters = nonnullclusters,
                     prevalence_set 	= prevalence_set,
                     common_shape = common_shape,
                     scale_con =	 scale_con,
                     acc_period = acc_period,
                     time_censor = time_censor
)

for(i in 1:7){write.dta(mygendata[[i]], paste("c:/temp/mydata", i, ".dta", sep = ""))}
