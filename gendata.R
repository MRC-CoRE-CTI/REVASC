# IW gendata.R built from Wenyue's ma_bayesmeta.R
# Notes: seed arg removed - should only set seed once
#   added nonnullclusters, removed nSims
#  removed analysis and returning data

gendata <- function(N_tot, # total sample size 
                         target_hr, # target HR 
                         crossover, # anticipated cross-over 
                         nonnull, # which clusters are non-null
                         prevalence_set, # prevalence set
                         common_shape, # Weibull shape parameter (common for all clusters)
                         scale_con, # control arm scale parameter (vary by cluster)
                         acc_period, # recruitment period
                         time_censor # time to censoring
)
{
  K <- length(prevalence_set)
  adj_target_loghr <- log(target_hr)*(1-crossover)^2 # Conversion to the adjusted target log(HR)
  
  # Generate a fixed sample size set
  N_cluster <- round(N_tot * prevalence_set)

  adj_N_tot <- sum(N_cluster) # Actual total sample size
  acc_rate <- N_cluster/acc_period # Mean number of occurrences per month
  
  # True log(HR) set
  beta <- adj_target_loghr * nonnull

  # Weibull model parameters
  scale_exp <- as.numeric(K) # experimental arm scale parameter given control arm scale parameter
  for (i in 1:K){
    scale_exp[i] <- exp(log(scale_con[i]) - (1/common_shape) * beta[i]) 
  }
  
    # mean_beta <- data.frame(matrix(nrow = 1, ncol = K))
    # weight <- data.frame(matrix(nrow = 1, ncol = K))
    # average_beta <- as.numeric(1)
    # between_sd <- as.numeric(1)
    # coverage <- data.frame(matrix(nrow = 1, ncol = K))
    reject <- data.frame(matrix(nrow = 1, ncol = K))
    # var <- data.frame(matrix(nrow = 1, ncol = K))
    # squared_error <- data.frame(matrix(nrow = 1, ncol = K))
    # low <- data.frame(matrix(nrow = 1, ncol = K))
    # upp <- data.frame(matrix(nrow = 1, ncol = K))
    # wid_ci <- data.frame(matrix(nrow = 1, ncol = K))
    
    # Poisson process with rate parameter
    # (also a set of random arrival times as Poisson point process of intensity at the rate)
    inter_1 <- rexp(N_cluster[1], rate = acc_rate[1])
    arr_1 <- cumsum(inter_1)
    
    inter_2 <- rexp(N_cluster[2], rate = acc_rate[2])
    arr_2 <- cumsum(inter_2) 
    
    inter_3 <- rexp(N_cluster[3], rate = acc_rate[3])
    arr_3 <- cumsum(inter_3) 
    
    inter_4 <- rexp(N_cluster[4], rate = acc_rate[4])
    arr_4 <- cumsum(inter_4)
    
    inter_5 <- rexp(N_cluster[5], rate = acc_rate[5])
    arr_5 <- cumsum(inter_5)
    
    inter_6 <- rexp(N_cluster[6], rate = acc_rate[6])
    arr_6 <- cumsum(inter_6)
    
    inter_7 <- rexp(N_cluster[7], rate = acc_rate[7])
    arr_7 <- cumsum(inter_7)
    
    # Number of arrivals for each month (can be either earlier or delayed by 1-2 months)
    arr_month_1 <- numeric(as.integer(max(arr_1)) + 1)
    arr_month_2 <- numeric(as.integer(max(arr_2)) + 1)
    arr_month_3 <- numeric(as.integer(max(arr_3)) + 1)
    arr_month_4 <- numeric(as.integer(max(arr_4)) + 1)
    arr_month_5 <- numeric(as.integer(max(arr_5)) + 1)
    arr_month_6 <- numeric(as.integer(max(arr_6)) + 1)
    arr_month_7 <- numeric(as.integer(max(arr_7)) + 1)
    
    for (i in 1:length(arr_month_1)){
      arr_month_1[i] <- sum(arr_1 > i-1 & arr_1 <= i)
    }
    
    for (i in 1:length(arr_month_2)){
      arr_month_2[i] <- sum(arr_2 > i-1 & arr_2 <= i)
    }
    
    for (i in 1:length(arr_month_3)){
      arr_month_3[i] <- sum(arr_3 > i-1 & arr_3 <= i)
    }
    
    for (i in 1:length(arr_month_4)){
      arr_month_4[i] <- sum(arr_4 > i-1 & arr_4 <= i)
    }
    
    for (i in 1:length(arr_month_5)){
      arr_month_5[i] <- sum(arr_5 > i-1 & arr_5 <= i)
    }
    
    for (i in 1:length(arr_month_6)){
      arr_month_6[i] <- sum(arr_6 > i-1 & arr_6 <= i)
    }
    
    for (i in 1:length(arr_month_7)){
      arr_month_7[i] <- sum(arr_7 > i-1 & arr_7 <= i)
    }
    
    # Treatment arm allocation
    arm_1 <- NULL
    arm_2 <- NULL
    arm_3 <- NULL
    arm_4 <- NULL
    arm_5 <- NULL
    arm_6 <- NULL
    arm_7 <- NULL
    con_month_1 <- as.numeric(length(arr_month_1))
    exp_month_1 <- as.numeric(length(arr_month_1))
    con_month_2 <- as.numeric(length(arr_month_2))
    exp_month_2 <- as.numeric(length(arr_month_2))
    con_month_3 <- as.numeric(length(arr_month_3))
    exp_month_3 <- as.numeric(length(arr_month_3))
    con_month_4 <- as.numeric(length(arr_month_4))
    exp_month_4 <- as.numeric(length(arr_month_4))
    con_month_5 <- as.numeric(length(arr_month_5))
    exp_month_5 <- as.numeric(length(arr_month_5))
    con_month_6 <- as.numeric(length(arr_month_6))
    exp_month_6 <- as.numeric(length(arr_month_6))
    con_month_7 <- as.numeric(length(arr_month_7))
    exp_month_7 <- as.numeric(length(arr_month_7))
    
    for (i in 1:length(arr_month_1)){
      new_arm_1 <- rbinom(n = arr_month_1[i], size = 1, prob = 0.5)
      arm_1 <- append(arm_1, new_arm_1)
      con_month_1[i] <- sum(new_arm_1 == '0')
      exp_month_1[i] <- arr_month_1[i] - con_month_1[i]
    }
    
    for (i in 1:length(arr_month_2)){
      new_arm_2 <- rbinom(n = arr_month_2[i], size = 1, prob = 0.5)
      arm_2 <- append(arm_2, new_arm_2)
      con_month_2[i] <- sum(new_arm_2 == '0')
      exp_month_2[i] <- arr_month_2[i] - con_month_2[i]
    }
    
    for (i in 1:length(arr_month_3)){
      new_arm_3 <- rbinom(n = arr_month_3[i], size = 1, prob = 0.5)
      arm_3 <- append(arm_3, new_arm_3)
      con_month_3[i] <- sum(new_arm_3 == '0')
      exp_month_3[i] <- arr_month_3[i] - con_month_3[i]
    }
    
    for (i in 1:length(arr_month_4)){
      new_arm_4 <- rbinom(n = arr_month_4[i], size = 1, prob = 0.5)
      arm_4 <- append(arm_4, new_arm_4)
      con_month_4[i] <- sum(new_arm_4 == '0')
      exp_month_4[i] <- arr_month_4[i] - con_month_4[i]
    }
    
    for (i in 1:length(arr_month_5)){
      new_arm_5 <- rbinom(n = arr_month_5[i], size = 1, prob = 0.5)
      arm_5 <- append(arm_5, new_arm_5)
      con_month_5[i] <- sum(new_arm_5 == '0')
      exp_month_5[i] <- arr_month_5[i] - con_month_5[i]
    }
    
    for (i in 1:length(arr_month_6)){
      new_arm_6 <- rbinom(n = arr_month_6[i], size = 1, prob = 0.5)
      arm_6 <- append(arm_6, new_arm_6)
      con_month_6[i] <- sum(new_arm_6 == '0')
      exp_month_6[i] <- arr_month_6[i] - con_month_6[i]
    }
    
    for (i in 1:length(arr_month_7)){
      new_arm_7 <- rbinom(n = arr_month_7[i], size = 1, prob = 0.5)
      arm_7 <- append(arm_7, new_arm_7)
      con_month_7[i] <- sum(new_arm_7 == '0')
      exp_month_7[i] <- arr_month_7[i] - con_month_7[i]
    }
    
    # Generating Weibull survival times for the control patients
    Survtime_con_1 <- rweibull(sum(con_month_1), shape = common_shape, scale = scale_con[1])
    Survtime_con_2 <- rweibull(sum(con_month_2), shape = common_shape, scale = scale_con[2])
    Survtime_con_3 <- rweibull(sum(con_month_3), shape = common_shape, scale = scale_con[3])
    Survtime_con_4 <- rweibull(sum(con_month_4), shape = common_shape, scale = scale_con[4])
    Survtime_con_5 <- rweibull(sum(con_month_5), shape = common_shape, scale = scale_con[5])
    Survtime_con_6 <- rweibull(sum(con_month_6), shape = common_shape, scale = scale_con[6])
    Survtime_con_7 <- rweibull(sum(con_month_7), shape = common_shape, scale = scale_con[7])
    
    # And for the experimental...
    Survtime_exp_1 <- rweibull(sum(exp_month_1), shape = common_shape, scale = scale_exp[1])
    Survtime_exp_2 <- rweibull(sum(exp_month_2), shape = common_shape, scale = scale_exp[2])
    Survtime_exp_3 <- rweibull(sum(exp_month_3), shape = common_shape, scale = scale_exp[3])
    Survtime_exp_4 <- rweibull(sum(exp_month_4), shape = common_shape, scale = scale_exp[4])
    Survtime_exp_5 <- rweibull(sum(exp_month_5), shape = common_shape, scale = scale_exp[5])
    Survtime_exp_6 <- rweibull(sum(exp_month_6), shape = common_shape, scale = scale_exp[6])
    Survtime_exp_7 <- rweibull(sum(exp_month_7), shape = common_shape, scale = scale_exp[7])
    
    # Survival status (based on administrative censoring)
    if (!is.null(time_censor) == TRUE) {
      test_con_1 <- (Survtime_con_1 > time_censor) #Identify which times need right censoring.
      test_exp_1 <- (Survtime_exp_1 > time_censor)
      test_con_2 <- (Survtime_con_2 > time_censor)
      test_exp_2 <- (Survtime_exp_2 > time_censor)
      test_con_3 <- (Survtime_con_3 > time_censor)
      test_exp_3 <- (Survtime_exp_3 > time_censor)
      test_con_4 <- (Survtime_con_4 > time_censor)
      test_exp_4 <- (Survtime_exp_4 > time_censor)
      test_con_5 <- (Survtime_con_5 > time_censor)
      test_exp_5 <- (Survtime_exp_5 > time_censor)
      test_con_6 <- (Survtime_con_6 > time_censor)
      test_exp_6 <- (Survtime_exp_6 > time_censor)
      test_con_7 <- (Survtime_con_7 > time_censor)
      test_exp_7 <- (Survtime_exp_7 > time_censor)
      status_con_1 <- rep(1, length(Survtime_con_1)) #Initialize the Status variable.
      status_exp_1 <- rep(1, length(Survtime_exp_1))
      status_con_2 <- rep(1, length(Survtime_con_2)) 
      status_exp_2 <- rep(1, length(Survtime_exp_2))
      status_con_3 <- rep(1, length(Survtime_con_3)) 
      status_exp_3 <- rep(1, length(Survtime_exp_3))
      status_con_4 <- rep(1, length(Survtime_con_4)) 
      status_exp_4 <- rep(1, length(Survtime_exp_4))
      status_con_5 <- rep(1, length(Survtime_con_5)) 
      status_exp_5 <- rep(1, length(Survtime_exp_5))
      status_con_6 <- rep(1, length(Survtime_con_6)) 
      status_exp_6 <- rep(1, length(Survtime_exp_6))
      status_con_7 <- rep(1, length(Survtime_con_7)) 
      status_exp_7 <- rep(1, length(Survtime_exp_7))
      #For all observations that need to be right censored, set the time
      # value to the right censor value and set the status value to 0
      # (indicating right censoring).  Status=1 implies observed event.
      # Status=0 implies right censored event.
      Survtime_con_1[test_con_1] <- time_censor
      Survtime_exp_1[test_exp_1] <- time_censor
      status_con_1[test_con_1] <- 0
      status_exp_1[test_exp_1] <- 0
      Survtime_con_2[test_con_2] <- time_censor
      Survtime_exp_2[test_exp_2] <- time_censor
      status_con_2[test_con_2] <- 0
      status_exp_2[test_exp_2] <- 0
      Survtime_con_3[test_con_3] <- time_censor
      Survtime_exp_3[test_exp_3] <- time_censor
      status_con_3[test_con_3] <- 0
      status_exp_3[test_exp_3] <- 0
      Survtime_con_4[test_con_4] <- time_censor
      Survtime_exp_4[test_exp_4] <- time_censor
      status_con_4[test_con_4] <- 0
      status_exp_4[test_exp_4] <- 0
      Survtime_con_5[test_con_5] <- time_censor
      Survtime_exp_5[test_exp_5] <- time_censor
      status_con_5[test_con_5] <- 0
      status_exp_5[test_exp_5] <- 0
      Survtime_con_6[test_con_6] <- time_censor
      Survtime_exp_6[test_exp_6] <- time_censor
      status_con_6[test_con_6] <- 0
      status_exp_6[test_exp_6] <- 0
      Survtime_con_7[test_con_7] <- time_censor
      Survtime_exp_7[test_exp_7] <- time_censor
      status_con_7[test_con_7] <- 0
      status_exp_7[test_exp_7] <- 0
    }
    #Create status variable if censor_value is not specified (in such a case
    # status = 1 for all observations.)
    if (is.null(time_censor) == TRUE) {
      status_con_1 <- rep(1, length(Survtime_con_1))
      status_exp_1 <- rep(1, length(Survtime_exp_1))
      status_con_2 <- rep(1, length(Survtime_con_2))
      status_exp_2 <- rep(1, length(Survtime_exp_2))
      status_con_3 <- rep(1, length(Survtime_con_3))
      status_exp_3 <- rep(1, length(Survtime_exp_3))
      status_con_4 <- rep(1, length(Survtime_con_4))
      status_exp_4 <- rep(1, length(Survtime_exp_4))
      status_con_5 <- rep(1, length(Survtime_con_5))
      status_exp_5 <- rep(1, length(Survtime_exp_5))
      status_con_6 <- rep(1, length(Survtime_con_6))
      status_exp_6 <- rep(1, length(Survtime_exp_6))
      status_con_7 <- rep(1, length(Survtime_con_7))
      status_exp_7 <- rep(1, length(Survtime_exp_7))
    }
    
    # Construct survival datasets (overall)
    data_1 <- as.data.frame(cbind(arr_1, arm_1))
    data_1$id <- seq(1, N_cluster[1], 1)
    data_1[data_1$arm_1 == '0', "survt"] <- Survtime_con_1
    data_1[data_1$arm_1 == '1', "survt"] <- Survtime_exp_1
    data_1[data_1$arm_1 == '0', "status"] <- status_con_1
    data_1[data_1$arm_1 == '1', "status"] <- status_exp_1
    
    data_2 <- as.data.frame(cbind(arr_2, arm_2))
    data_2$id <- seq(1, N_cluster[2], 1)
    data_2[data_2$arm_2 == '0', "survt"] <- Survtime_con_2
    data_2[data_2$arm_2 == '1', "survt"] <- Survtime_exp_2
    data_2[data_2$arm_2 == '0', "status"] <- status_con_2
    data_2[data_2$arm_2 == '1', "status"] <- status_exp_2
    
    data_3 <- as.data.frame(cbind(arr_3, arm_3))
    data_3$id <- seq(1, N_cluster[3], 1)
    data_3[data_3$arm_3 == '0', "survt"] <- Survtime_con_3
    data_3[data_3$arm_3 == '1', "survt"] <- Survtime_exp_3
    data_3[data_3$arm_3 == '0', "status"] <- status_con_3
    data_3[data_3$arm_3 == '1', "status"] <- status_exp_3
    
    data_4 <- as.data.frame(cbind(arr_4, arm_4))
    data_4$id <- seq(1, N_cluster[4], 1)
    data_4[data_4$arm_4 == '0', "survt"] <- Survtime_con_4
    data_4[data_4$arm_4 == '1', "survt"] <- Survtime_exp_4
    data_4[data_4$arm_4 == '0', "status"] <- status_con_4
    data_4[data_4$arm_4 == '1', "status"] <- status_exp_4
    
    data_5 <- as.data.frame(cbind(arr_5, arm_5))
    data_5$id <- seq(1, N_cluster[5], 1)
    data_5[data_5$arm_5 == '0', "survt"] <- Survtime_con_5
    data_5[data_5$arm_5 == '1', "survt"] <- Survtime_exp_5
    data_5[data_5$arm_5 == '0', "status"] <- status_con_5
    data_5[data_5$arm_5 == '1', "status"] <- status_exp_5
    
    data_6 <- as.data.frame(cbind(arr_6, arm_6))
    data_6$id <- seq(1, N_cluster[6], 1)
    data_6[data_6$arm_6 == '0', "survt"] <- Survtime_con_6
    data_6[data_6$arm_6 == '1', "survt"] <- Survtime_exp_6
    data_6[data_6$arm_6 == '0', "status"] <- status_con_6
    data_6[data_6$arm_6 == '1', "status"] <- status_exp_6
    
    data_7 <- as.data.frame(cbind(arr_7, arm_7))
    data_7$id <- seq(1, N_cluster[7], 1)
    data_7[data_7$arm_7 == '0', "survt"] <- Survtime_con_7
    data_7[data_7$arm_7 == '1', "survt"] <- Survtime_exp_7
    data_7[data_7$arm_7 == '0', "status"] <- status_con_7
    data_7[data_7$arm_7 == '1', "status"] <- status_exp_7

	return(list(data_1,data_2,data_3,data_4,data_5,data_6,data_7))

}
