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

ANOVA_priors <- function(dat_rated, dat_casecov, 
                         rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
                         case.covs = c("self.iq1", "self.iq5", "self.iq6",
                                       "grade1", "grade2", "grade4"),
                         IDout, IDin, IDgroup, precision, default_prior) {
  # make sure names(dat_casecov) <- c("id", "group.id", "self.iq1", "self.iq5",
  # "self.iq6", "grade1", "grade2", "grade4")
  
  priors <- default_prior
  library(TripleR)
  
  matNames_c <- c(paste0(rep(rr.vars, each = 2), c("_out", "_in")), case.covs)
  covMat_c <- matrix(0, length(rr.vars)*2 + length(case.covs), length(rr.vars)*2 + length(case.covs),
                     dimnames = list(matNames_c, matNames_c))
  matNames_d <- paste0(rep(rr.vars, each = 2), c("_ij", "_ji"))
  covMat_d <- matrix(0, length(rr.vars)*2, length(rr.vars)*2,
                     dimnames = list(matNames_d, matNames_d))
  
  for (r in 1:length(rr.vars)) {
    uv.formula <- paste0(rr.vars[r], "~", IDout, "*", IDin)
    if (!missing(IDgroup)) uv.formula <- paste0(uv.formula, "|", IDgroup, collapse = "")
    uvSRM <- RR(as.formula(uv.formula), data = dat_rated, na.rm = T)
    uv.gEsts <- cbind(uvSRM$varComp.groups[uvSRM$varComp.groups$type != "error variance",], 
                      rr.var = rr.vars[r])
    uv.varNames <- unique(uv.gEsts$type[grep(" variance", uv.gEsts$type)])
    uv.covNames <- unique(uv.gEsts$type[grep("covariance", uv.gEsts$type)])
    
    # rescale negative variances (inadmissable cases) to 0
    uv.gEsts[uv.gEsts$type %in% uv.varNames,]$estimate <- 
      ifelse(uv.gEsts[uv.gEsts$type %in% uv.varNames,]$estimate < 0, 0, 
             uv.gEsts[uv.gEsts$type %in% uv.varNames,]$estimate)
    
    # save variances per variable (weighted mean)
    uv.var <- sapply(uv.varNames, function(type) {
      n <- uv.gEsts[uv.gEsts$type == type,]$group.size
      weighted.mean(uv.gEsts[uv.gEsts$type == type,]$estimate, w = n-1)
    })
    
    # save covariances per variable (weighted mean)
    uv.cov <- sapply(uv.covNames, function(type) {
      n <- uv.gEsts[uv.gEsts$type == type,]$group.size
      weighted.mean(uv.gEsts[uv.gEsts$type == type,]$estimate, w = n-1)
    })
    
    # construct covariance matrices (to compute correlation matrices)
    covMat_c[paste0(rr.vars[r], "_out"), paste0(rr.vars[r], "_out")] <- uv.var["actor variance"]
    covMat_c[paste0(rr.vars[r], "_in"), paste0(rr.vars[r], "_in")] <- uv.var["partner variance"]
    covMat_c[paste0(rr.vars[r], "_out"), paste0(rr.vars[r], "_in")] <- covMat_c[paste0(rr.vars[r], "_in"), paste0(rr.vars[r], "_out")] <- 
      uv.cov["actor-partner covariance"]
    
    covMat_d[paste0(rr.vars[r], "_ij"), paste0(rr.vars[r], "_ij")] <- 
      covMat_d[paste0(rr.vars[r], "_ji"), paste0(rr.vars[r], "_ji")] <- uv.var["relationship variance"]
    covMat_d[paste0(rr.vars[r], "_ij"), paste0(rr.vars[r], "_ji")] <- 
      covMat_d[paste0(rr.vars[r], "_ji"), paste0(rr.vars[r], "_ij")] <- uv.cov["relationship covariance"]
    
    ## SD priors
    priors$rr_in_t[rr.vars[r],  "m"] <- sqrt(uv.var["partner variance"])
    priors$rr_out_t[rr.vars[r], "m"] <- sqrt(uv.var["actor variance"])
    priors$rr_rel_t[rr.vars[r], "m"] <- sqrt(uv.var["relationship variance"])
    
    priors$rr_in_t[rr.vars[r],  "sd"] <- priors$rr_out_t[rr.vars[r], "sd"] <- priors$rr_rel_t[rr.vars[r], "sd"] <- precision
    
    
    combi_dat <- merge(dat_casecov, uvSRM$effects, by = c("id", "group.id"))
    names(combi_dat) <- gsub(pattern = "[.]a", replacement = "_out", 
                             gsub(pattern = "[.]p", replacement = "_in", names(combi_dat)))
    
    ## correlations with case-level covariates
    for (c in 1:length(case.covs)) {
      covMat_c[paste0(rr.vars[r], "_out"), case.covs[c]] <- covMat_c[case.covs[c], paste0(rr.vars[r], "_out")] <- 
        parCor(combi_dat[,case.covs[c]], combi_dat[,paste0(rr.vars[r], "_out")], combi_dat$group.id)$par.cor # covariance with outgoing effects
      covMat_c[paste0(rr.vars[r], "_out"), case.covs[c]] <- covMat_c[case.covs[c], paste0(rr.vars[r], "_out")] <- 
        parCor(combi_dat[,case.covs[c]], combi_dat[,paste0(rr.vars[r], "_out")], combi_dat$group.id)$par.cor # covariance with incoming effects
      
      # SDs for each covariate in diagonal
      covMat_c[case.covs[c], case.covs[c]] <- sd(dat_casecov[,case.covs[c]])
    }
  }
  
  
  
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