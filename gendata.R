# Alternative to gendata.R
# Generate data in one data frame
# Allows for staggered entry
# Assumes entry times are uniformly distributed over the accrual period
# IW 30apr2026
# Checked by check_gendata_R: does give correct answers

gendata <- function(K, # subgroups
                    N_tot, # total sample size 
                    target_hr, # target HR 
                    crossover, # anticipated cross-over 
                    prevalence_set, # prevalence set
                    common_shape, # Weibull shape parameter (common for all subgroups)
                    scale_con, # control arm scale parameter (vary by subgroup)
                    nonnull, # which subgroups are non-null
                    acc_period, # recruitment period
                    time_censor, # time to censoring
                    palloc=0.5 # p(rand to treatment)
)
{
# Calculations by subgroup
subgroupid = c(1:K)

# Generate a fixed sample size set
N_subgroup <- round(N_tot * prevalence_set)
adj_N_tot <- sum(N_subgroup) # Actual total sample size

adj_target_loghr <- log(target_hr)*(1-crossover)^2 # Conversion to the adjusted target log(HR)
beta <- adj_target_loghr * nonnull

# Create data frame of individuals
subgroups = as.data.frame(t(rbind(subgroupid, N_subgroup, scale_con, beta)))
df <- subgroups[rep(1:nrow(subgroups), subgroups$N_subgroup), ]
df$N_subgroup <- NULL

# Time of accrual
df$start = runif(adj_N_tot,min=0,max=acc_period)

# Randomise
df$rand = rbinom(adj_N_tot,size=1,prob=palloc)

# Survival time
df$end = df$start+rweibull(adj_N_tot, shape = common_shape, scale = df$scale_con*exp(-df$beta*df$rand/common_shape))

# Censoring
df$status = df$end < acc_period + time_censor
df$survt = pmin(df$end, acc_period + time_censor) - df$start

# Tidy up
df$scale_con <- NULL
df$beta <- NULL
df$end <- NULL

return(df)
}

