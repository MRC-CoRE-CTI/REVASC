# check_gendata.R
# create big data sets that we can analyse in Stata to confirm control S(5)

K <- length(prevalence_set)
mygendata <- gendata(
    N_tot = 3000,
    target_hr = target_hr, 
    crossover = crossover, 
    prevalence_set = prevalence_set,
    common_shape = common_shape,
    scale_con = scale_con,
    nonnull = rep(1, K),
    acc_period = acc_period,
    time_censor = time_censor
)

# Export to Stata
library(foreign)
write.dta(mygendata, "gendata.dta")

# various analyses and descriptive statistics
table(mygendata$subgroupid, mygendata$rand)
# check proportions are correct

coxph(Surv(survt,status)~rand+factor(subgroupid), data=mygendata)
# check log HR is correct

fit <- survfit(Surv(survt, status) ~ subgroupid, data=mygendata)
plot(fit, xlab = "Time", ylab = "Survival probability")
# check shapes

myprepdata <- prepdata(mygendata) 
myprepdata

myanadata<- anadata(myprepdata)
myanadata
