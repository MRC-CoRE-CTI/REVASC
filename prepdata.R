# IW prepdata.R built from Wenyue's ma_bayesmeta.R
# takes cluster-specific datasets and returns logHRs and SEs

prepdata<-function(gendata) {
    # First-stage analysis
    # Cox regression
    model_1 <- tryCatch(model_1 <- coxph(Surv(survt, status) ~ arm_1, data = gendata[[1]]),
                        error = function(e) NULL)
    model_2 <- tryCatch(model_2 <- coxph(Surv(survt, status) ~ arm_2, data = gendata[[2]]),
                        error = function(e) NULL)
    model_3 <- tryCatch(model_3 <- coxph(Surv(survt, status) ~ arm_3, data = gendata[[3]]),
                        error = function(e) NULL)
    model_4 <- tryCatch(model_4 <- coxph(Surv(survt, status) ~ arm_4, data = gendata[[4]]),
                        error = function(e) NULL)
    model_5 <- tryCatch(model_5 <- coxph(Surv(survt, status) ~ arm_5, data = gendata[[5]]),
                        error = function(e) NULL)
    model_6 <- tryCatch(model_6 <- coxph(Surv(survt, status) ~ arm_6, data = gendata[[6]]),
                        error = function(e) NULL)
    model_7 <- tryCatch(model_7 <- coxph(Surv(survt, status) ~ arm_7, data = gendata[[7]]),
                        error = function(e) NULL)
    
    # Weibull regression
    # model_1 <- tryCatch(model_1 <- WeibullReg(Surv(survt, status) ~ arm_1, data = data_1), 
    #                     error = function(e) NULL)
    # model_2 <- tryCatch(model_2 <- WeibullReg(Surv(survt, status) ~ arm_2, data = data_2), 
    #                     error = function(e) NULL)
    # model_3 <- tryCatch(model_3 <- WeibullReg(Surv(survt, status) ~ arm_3, data = data_3), 
    #                     error = function(e) NULL)
    # model_4 <- tryCatch(model_4 <- WeibullReg(Surv(survt, status) ~ arm_4, data = data_4), 
    #                     error = function(e) NULL)
    # model_5 <- tryCatch(model_5 <- WeibullReg(Surv(survt, status) ~ arm_5, data = data_5), 
    #                     error = function(e) NULL)
    # model_6 <- tryCatch(model_6 <- WeibullReg(Surv(survt, status) ~ arm_6, data = data_6), 
    #                     error = function(e) NULL)
    # model_7 <- tryCatch(model_7 <- WeibullReg(Surv(survt, status) ~ arm_7, data = data_7), 
    #                     error = function(e) NULL)
    
    if(!is.null(model_1) | !is.null(model_2) | !is.null(model_3) | !is.null(model_4) | !is.null(model_5) | !is.null(model_6) | !is.null(model_7)){
      # Summary stats from Cox regression
      beta_1 <- as.numeric(model_1$coefficients)
      se_1 <- as.numeric(sqrt(model_1$var))
      
      beta_2 <- as.numeric(model_2$coefficients)
      se_2 <- as.numeric(sqrt(model_2$var))
      
      beta_3 <- as.numeric(model_3$coefficients)
      se_3 <- as.numeric(sqrt(model_3$var))
      
      beta_4 <- as.numeric(model_4$coefficients)
      se_4 <- as.numeric(sqrt(model_4$var))
 
      beta_5 <- as.numeric(model_5$coefficients)
      se_5 <- as.numeric(sqrt(model_5$var))
      
      beta_6 <- as.numeric(model_6$coefficients)
      se_6 <- as.numeric(sqrt(model_6$var))
      
      beta_7 <- as.numeric(model_7$coefficients)
      se_7 <- as.numeric(sqrt(model_7$var))
      
      # Summary stats from Weibull regression
      # beta_1 <- model_1$coef['arm_1', 'Estimate']
      # se_1 <- model_1$coef['arm_1', 'SE']
 
      # beta_2 <- model_2$coef['arm_2', 'Estimate']
      # se_2 <- model_2$coef['arm_2', 'SE']

      # beta_3 <- model_3$coef['arm_3', 'Estimate']
      # se_3 <- model_3$coef['arm_3', 'SE']

      # beta_4 <- model_4$coef['arm_4', 'Estimate']
      # se_4 <- model_4$coef['arm_4', 'SE']

      # beta_5 <- model_5$coef['arm_5', 'Estimate']
      # se_5 <- model_5$coef['arm_5', 'SE']

      # beta_6 <- model_6$coef['arm_6', 'Estimate']
      # se_6 <- model_6$coef['arm_6', 'SE']

      # beta_7 <- model_7$coef['arm_7', 'Estimate']
      # se_7 <- model_7$coef['arm_7', 'SE']
      
      beta_hat <- c(beta_1, beta_2, beta_3, beta_4, beta_5, beta_6, beta_7)
      sd.within <- c(se_1, se_2, se_3, se_4, se_5, se_6, se_7)
      var.within <- sd.within^2
      
	  results <- list(beta_hat,sd.within)
}
    else{
      results <- list(
        rep(NA, 1)
      )
    }

return(results)

}
