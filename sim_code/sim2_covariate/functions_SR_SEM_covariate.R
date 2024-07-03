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

# iq.data = readRDS("rated_iq.rds")
# covariate.data = readRDS("covariate_dat.rds")
# default_prior <- lavaan.srm::srm_priors(data = iq.data[4:6], case_data = covariate.data[3:8])

ANOVA_priors <- function(rr.data, case.data, 
                         rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
                         case.covs = c("self.iq1", "self.iq5", "self.iq6",
                                       "grade1", "grade2", "grade4"),
                         IDout, IDin, IDgroup, 
                         precision = 0.1, default_prior) {
  # make sure names(case.data) <- c("id", "group.id", "self.iq1", "self.iq5",
  # "self.iq6", "grade1", "grade2", "grade4")
  
  priors <- default_prior
  library(TripleR)
  
  matNames_c <- c(paste0(rep(rr.vars, each = 2), c("_out", "_in")), case.covs)
  corMat_c <- matrix(0, length(rr.vars)*2 + length(case.covs), length(rr.vars)*2 + length(case.covs),
                     dimnames = list(matNames_c, matNames_c)) # correlation matrix at case level
  diag(corMat_c) <- 1
  matNames_d <- paste0(rep(rr.vars, each = 2), c("_ij", "_ji"))
  covMat_d <- matrix(0, length(rr.vars)*2, length(rr.vars)*2,
                     dimnames = list(matNames_d, matNames_d)) # covariance matrix at dyad level

  uv.SDs <- vector() # to use in bivariate SRMs
  
  for (rr in 1:length(rr.vars)) {
    uv.formula <- paste0(rr.vars[rr], "~", IDout, "*", IDin)
    if (!missing(IDgroup)) uv.formula <- paste0(uv.formula, "|", IDgroup, collapse = "")
    uvSRM <- RR(as.formula(uv.formula), data = rr.data, na.rm = T)
    uv.gEsts <- cbind(uvSRM$varComp.groups[uvSRM$varComp.groups$type != "error variance",], 
                      rr.var = rr.vars[rr])
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
    corMat_c[paste0(rr.vars[rr], "_out"), paste0(rr.vars[rr], "_in")] <- corMat_c[paste0(rr.vars[rr], "_in"), paste0(rr.vars[rr], "_out")] <- 
      uv.cov["actor-partner covariance"] / sqrt(uv.var["actor variance"]*uv.var["partner variance"])
    
    covMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[rr], "_ij")] <- 
      covMat_d[paste0(rr.vars[rr], "_ji"), paste0(rr.vars[rr], "_ji")] <- uv.var["relationship variance"]
    covMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[rr], "_ji")] <- 
      covMat_d[paste0(rr.vars[rr], "_ji"), paste0(rr.vars[rr], "_ij")] <- uv.cov["relationship covariance"]
    
    ## SD priors
    priors$rr_in_t[rr.vars[rr],  "m"] <- sqrt(uv.var["partner variance"])
    priors$rr_out_t[rr.vars[rr], "m"] <- sqrt(uv.var["actor variance"])
    priors$rr_rel_t[rr.vars[rr], "m"] <- sqrt(uv.var["relationship variance"])
    
    priors$rr_in_t[rr.vars[rr],  "sd"] <- priors$rr_out_t[rr.vars[rr], "sd"] <- priors$rr_rel_t[rr.vars[rr], "sd"] <- precision
    
    
    combi_dat <- merge(case.data, uvSRM$effects, by = c("id", "group.id"))
    names(combi_dat) <- gsub(pattern = "[.]a", replacement = "_out", 
                             gsub(pattern = "[.]p", replacement = "_in", names(combi_dat)))
    
    for (c in 1:length(case.covs)) {
      corMat_c[paste0(rr.vars[rr], "_out"), case.covs[c]] <- corMat_c[case.covs[c], paste0(rr.vars[rr], "_out")] <- 
        parCor(combi_dat[,case.covs[c]], combi_dat[,paste0(rr.vars[rr], "_out")], combi_dat$group.id)$par.cor # correlation with outgoing effects
      corMat_c[paste0(rr.vars[rr], "_in"), case.covs[c]] <- corMat_c[case.covs[c], paste0(rr.vars[rr], "_in")] <- 
        parCor(combi_dat[,case.covs[c]], combi_dat[,paste0(rr.vars[rr], "_in")], combi_dat$group.id)$par.cor # correlation with incoming effects
      #TODO change this to IDgroup --- make sure group IDs are "group.id"
      
      # SD priors for each covariate
      priors$case_cov_t[case.covs[c], "m"] <- sd(case.data[,case.covs[c]])
      priors$case_cov_t[case.covs[c], "sd"] <- precision
    }
    
    # saving to use for the bivariate SRMs
    uv.SDs[paste0(rr.vars[rr], "_out")] <- sqrt(uv.var["actor variance"])
    uv.SDs[paste0(rr.vars[rr], "_in")] <- sqrt(uv.var["partner variance"])
  }
  
  for (rr in 2:length(rr.vars)) {
    for (kk in 1:(rr-1)) {
      bv.formula <- paste0(rr.vars[rr], "+", rr.vars[kk], "~", IDout, "*", IDin)
      if (!missing(IDgroup)) bv.formula <- paste0(bv.formula, "|", IDgroup, collapse = "")
      bvSRM <- RR(as.formula(paste0(bv.formula)), data = rr.data, na.rm = T)
      bv.gEsts <- do.call("rbind", lapply(1:length(bvSRM$groups), 
                                          function (g) cbind(bvSRM$groups[[g]]$bivariate, 
                                                             Group = g, 
                                                             rr.var = paste0(rr.vars[rr],",",rr.vars[kk]))))
      bv.covNames <- unique(bv.gEsts$type[grep("covariance", bv.gEsts$type)])
      
      ## weighted mean of the covariances
      bv.cov <- sapply(bv.covNames, function(type) {
        #FIXME: TDJ: Observed "n" in data might be too large if listwise deleted by TripleR.
        #       Can TripleR add $group.size to result? Why is it missing for bivariate?
        #TODO: TDJ: weight dyad level by number of dyads?
        n <- sapply(unique(rr.data[, IDgroup]),
                    function(x) length(unique(rr.data[rr.data[,IDgroup] == x,][,IDout])), 
                    simplify = TRUE)
        weighted.mean(bv.gEsts[bv.gEsts$type == type,]$estimate, w = n-1)
      })
      
      corMat_c[paste0(rr.vars[rr], "_out"), paste0(rr.vars[kk], "_out")] <- 
        corMat_c[paste0(rr.vars[kk], "_out"), paste0(rr.vars[rr], "_out")] <- 
        bv.cov["actor-actor covariance"] / (uv.SDs[paste0(rr.vars[rr], "_out")]* uv.SDs[paste0(rr.vars[kk], "_out")])
      corMat_c[paste0(rr.vars[rr], "_in"), paste0(rr.vars[kk], "_in")] <- 
        corMat_c[paste0(rr.vars[kk], "_in"), paste0(rr.vars[rr], "_in")] <- 
        bv.cov["partner-partner covariance"] / (uv.SDs[paste0(rr.vars[rr], "_in")]* uv.SDs[paste0(rr.vars[kk], "_in")])
      corMat_c[paste0(rr.vars[rr], "_out"), paste0(rr.vars[kk], "_in")] <- 
        corMat_c[paste0(rr.vars[kk], "_in"), paste0(rr.vars[rr], "_out")] <- 
        bv.cov["actor-partner covariance"] / (uv.SDs[paste0(rr.vars[rr], "_out")]* uv.SDs[paste0(rr.vars[kk], "_in")])
      corMat_c[paste0(rr.vars[rr], "_in"), paste0(rr.vars[kk], "_out")] <- 
        corMat_c[paste0(rr.vars[kk], "_out"), paste0(rr.vars[rr], "_in")] <- 
        bv.cov["partner-actor covariance"] / (uv.SDs[paste0(rr.vars[rr], "_in")]* uv.SDs[paste0(rr.vars[kk], "_out")])
      
      covMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[kk], "_ij")] <- 
        covMat_d[paste0(rr.vars[kk], "_ij"), paste0(rr.vars[rr], "_ij")] <- 
        covMat_d[paste0(rr.vars[rr], "_ji"), paste0(rr.vars[kk], "_ji")] <- 
        covMat_d[paste0(rr.vars[kk], "_ji"), paste0(rr.vars[rr], "_ji")] <- 
        bv.cov["intrapersonal relationship covariance"]
      
      covMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[kk], "_ji")] <- 
        covMat_d[paste0(rr.vars[kk], "_ji"), paste0(rr.vars[rr], "_ij")] <- 
        covMat_d[paste0(rr.vars[rr], "_ji"), paste0(rr.vars[kk], "_ij")] <- 
        covMat_d[paste0(rr.vars[kk], "_ij"), paste0(rr.vars[rr], "_ji")] <- 
        bv.cov["interpersonal relationship covariance"]
    }
  }
  
  corMat_d <- cov2cor(covMat_d)
  
  for (c in 2:length(case.covs)) { # insert correlations between case-level covariates
    for (cc in 1:(c-1)) {
      corMat_c[case.covs[c], case.covs[cc]] <- corMat_c[case.covs[cc], case.covs[c]] <-
        parCor(combi_dat[,case.covs[c]], combi_dat[,case.covs[cc]], combi_dat$group.id)$par.cor
    }
  }
  
  # if case-level correlation > abs(0.9) then rescale to +/-0.9 and populate srm_priors matrices
  for (x in 2:nrow(corMat_c)) {
    for (y in 1:(x-1)) {
      corMat_c[x,y] <- ifelse(corMat_c[x,y]> 0.9, 0.9, 
                              ifelse(corMat_c[x,y] < -0.9, -0.9, corMat_c[x,y]))
      hyperpars_c <- optim(par = c(1.5, 1.5), fn = minLoss,
                           targetCorr = corMat_c[x,y], accept.dev = precision,
                           method = "L-BFGS-B", lower = 0)$par
      priors$case_beta[x,y] <- hyperpars_c[1] # lower triangle
      priors$case_beta[y,x] <- hyperpars_c[2] # upper triangle
    }
  }
  
  # if dyad-level correlation > abs(0.9) then rescale to +/-0.9, then populate srm_priors matrices 
  corMat_d[row(corMat_d)!=col(corMat_d)] <- ifelse(corMat_d[row(corMat_d)!=col(corMat_d)] > 0.9, 0.9,
                                                   ifelse(corMat_d[row(corMat_d)!=col(corMat_d)] < -0.9,
                                                          -0.9, corMat_d[row(corMat_d)!=col(corMat_d)]))
  
  for (rr in 1:length(rr.vars)) {
    
    ### dyadic-- diagonal
    dyadic_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                              targetCorr = corMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[rr], "_ji")],
                              accept.dev = precision, 
                              method = "L-BFGS-B", lower = 0)$par
    priors$rr_beta_a[rr,rr] <- dyadic_hyperpars[1]
    priors$rr_beta_b[rr,rr] <- dyadic_hyperpars[2]
    
    if (rr > 1L) for (kk in 1:(rr-1)) { 
      ### intra-- below diagonal
      intra_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                               targetCorr = corMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[kk], "_ij")],
                               accept.dev = precision, 
                               method = "L-BFGS-B", lower = 0)$par
      priors$rr_beta_a[rr,kk] <- intra_hyperpars[1]
      priors$rr_beta_b[rr,kk] <- intra_hyperpars[2]
      ### inter-- above diagonal
      inter_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                               targetCorr = corMat_d[paste0(rr.vars[rr], "_ij"), paste0(rr.vars[kk], "_ji")],
                               accept.dev = precision, 
                               method = "L-BFGS-B", lower = 0)$par
      priors$rr_beta_a[kk,rr] <- inter_hyperpars[1]
      priors$rr_beta_b[kk,rr] <- inter_hyperpars[2]
    }
  }
  
  priors
} #TODO check again for bugs

ANOVA_priors(rr.data = iq.data, case.data = covariate.data, 
             IDout = "ego", IDin = "alter", IDgroup = "group",
             default_prior = default_prior)

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