# setup.R - read files and set parms

# User-specific setup
setwd("C:/ian/git/REVASC")
# End of user-specific setup

library(Rlab)
library(survival)
# library(SurvRegCensCov)
library(foreach)
library(doParallel)
library(purrr)
library(bayesmeta)

source("gendata.R")
source("prepdata.R")
source("anadata.R")
source("manyreps.R")

target_hr = 0.7
crossover = 0.08
prevalence_set = c(0.25, 0.04, 0.06, 0.27, 0.24, 0.09, 0.05)
common_shape = 1.1
cont_cuminc_5y = c(.623,.490,.398,.609,.392,.720,.689) # Weiqi's 5-year cum inc for PCI group
scale_con = 60 * (-log(1-cont_cuminc_5y))^(-1/common_shape) # Weibull calculation
acc_period = 36
time_censor = 48

set.seed(101)

K <- length(prevalence_set)
