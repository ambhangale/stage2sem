## Aditi M. Bhangale
## Last updated: 2 July 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Incorporating a level-specific covariate

## PREPARE STAGE1 `lavaan.srm` RESULTS FOR STAGE 2 + STAGE 2 `lavaan.srm`

# rm(list = ls())

setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem/sim_code/sim2_covariate/")
# getwd()

#################################
# FUNCTIONS FOR SIMULATIONS ----
#################################

# function 0: generate level-specific (co)variance matrices----
getSigma <- function(level) {
  
}
#----

# function 1: simulate data for a single group----

#----

# function 2: simulate data for multiple groups----

#----

# function 3: minimum loss function to optimise over to choose hyperparameters for beta priors----

minLoss <- function(par = c(1.5, 1.5), targetCorr, accept.dev) {
  # targetCorr = population correlation value we want to estimate
  # accept.dev = average deviation from targetCorr we are willing to accept
  
  # default mean (location) of beta distribution (by default, is 0.5)
  defaultM <- par[1] / (par[1] + par[2])
  # default SD of beta distribution (by default, is 0.25)
  defaultSD <- sqrt(par[1]*par[2] / ((par[1]+par[2])^2 * (par[1] + par[2] + 1))) 
  
  
  # target correlation (location of beta) and accepted deviation around it (accept.dev)
  # convert to beta scale
  targetM <- (targetCorr + 1) / 2 # convert targetCorr to beta scale (location of beta)
  targetSD <- accept.dev / 2 # convert accept.dev to beta scale (SD of beta)
  
  mSE <- (defaultM - targetM)^2 # squared error for mean
  sdSE <- (defaultSD - targetSD)^2 # squared error for sd
  
  mSE + sdSE # aim: minimise the sum of squared error (loss function)-- thus, use optim()
}

# optim(par = c(1.5, 1.5), fn = minLoss, targetCorr = 0.3, accept.dev = 0.1,
#       method = "L-BFGS-B", lower = 0)

#----

# function 4: prophetic priors

#----

# function 5: ANOVA-based priors

ANOVA_priors <- function(dat_rated, dat_casecov) {
  
  
}

#----

# function 6: set customised priors for MCMC stage----

#----

# function 7: stage 1 of lavaan.srm----
#TODO save stage1 files, so that T doesn't have to run them again when trying different
# methods for test stats, robust ACOV, etc
#----

# function 8: stage 2 of lavaan.srm----

#----

## priors: default, ANOVA, prophetic
### ANOVA priors: the prior locations among the case-level covaraites are easy. just run `cor()`
### then run a univariate srm for each case-dyad level pair, treating the case-level covariate
### as self ratings