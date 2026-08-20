# IW prepdata.R built from Wenyue's ma_bayesmeta.R
# takes subgroup-specific datasets and returns logHRs and SEs
# revised 20aug2026 for single dataset

prepdata<-function(gendata) {
  # First-stage analysis
  # Cox regression
  K <- max(gendata$subgroupid)
  beta.hat <- rep(NA,K)
  sd.within <- rep(NA,K)
  for(i in 1:K){
    print(i)
    model <- tryCatch(
      coxph(Surv(survt, status) ~ rand, data = subset(mygendata,subgroupid==i)), 
      error = function(e) NULL)
    if(!is.null(model)) {
      # Summary stats from Cox regression
      beta.hat[i] <- as.numeric(model$coefficients)
      sd.within[i] <- as.numeric(sqrt(model$var))
    }
  }
  results <- list(beta.hat=beta.hat,sd.within=sd.within)
  return(results)
}
