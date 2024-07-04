## Aditi M. Bhangale
## Last updated: 4 July 2024

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
getSigma <- function(return_mats = TRUE) {
    popcov.vec_c <- c(0.367178156002207,0.0282638172371223,0.232566187414696,-0.00042303152370159,
                      0.145404838183355,0.00256135404013844,0.19650649315067,0.141654940954547,
                      0.0832506443945378,0.0179926118180312,0.0516975118928625,-0.00815220571888204,
                      0.0282638172371223,0.118625041049478,0.00436854440779768,0.0994599593219533,
                      0.00611922017653478,0.0940880968540291,0.0772761368510366,0.0918692560924219,
                      0.0809180729246042,0.0826428880708735,0.0796935097239754,0.06850886275001,
                      0.232566187414696,0.00436854440779768,0.470529300000656,-0.0119247565252019,
                      0.309544649467953,-0.0430818983370867,0.122650544584262,0.248888862834933,
                      0.158902022440876,0.00522078558607705,0.0818404324127288,-0.0155742030238611,
                      -0.00042303152370159,0.0994599593219533,-0.0119247565252019,0.197907713215963,
                      0.0107240331607381,0.178920190736849,0.0914084704107476,0.0991671583731612,
                      0.114559778372169,0.192815472860378,0.168808200401105,0.140288023680691,
                      0.145404838183355,0.00611922017653478,0.309544649467953,0.0107240331607381,
                      0.419707598143518,0.0021001712597054,0.0906398691171177,0.16010599825563,
                      0.293721236120613,0.0234391069284459,0.0570586359951145,0.0202568439692539,
                      0.00256135404013844,0.0940880968540291,-0.0430818983370867,0.178920190736849,
                      0.0021001712597054,0.189996325669525,0.0820405997800706,0.0897308819246081,
                      0.103607444797212,0.179600581343282,0.168834134277284,0.144804682212415,
                      0.19650649315067,0.0772761368510366,0.122650544584262,0.0914084704107476,
                      0.0906398691171177,0.0820405997800706,0.73648636471313,0.364561502339407,
                      0.346920016798681,0.102077086997001,0.0810390897699736,0.0658465393457167,
                      0.141654940954547,0.0918692560924219,0.248888862834933,0.0991671583731612,
                      0.16010599825563,0.0897308819246081,0.364561502339407,0.650560216188368,
                      0.40072530226852,0.123389906869327,0.117031365329952,0.0602534410885451,
                      0.0832506443945378,0.0809180729246042,0.158902022440876,0.114559778372169,
                      0.293721236120613,0.103607444797212,0.346920016798681,0.40072530226852,
                      0.714613137577745,0.135069195417365,0.1605460671607,0.134148244398022,
                      0.0179926118180312,0.0826428880708735,0.00522078558607705,0.192815472860378,
                      0.0234391069284459,0.179600581343282,0.102077086997001,0.123389906869327,
                      0.135069195417365,0.508493900461932,0.287147572145552,0.208055233197341,
                      0.0516975118928625,0.0796935097239754,0.0818404324127288,0.168808200401105,
                      0.0570586359951145,0.168834134277284,0.0810390897699736,0.117031365329952,
                      0.1605460671607,0.287147572145552,0.516170394363591,0.280945879685804,
                      -0.00815220571888204,0.06850886275001,-0.0155742030238611,0.140288023680691,
                      0.0202568439692539,0.144804682212415,0.0658465393457167,0.0602534410885451,
                      0.134148244398022,0.208055233197341,0.280945879685804,0.475466332681621)
      
    pop.names_c <- c(paste0(rep(c("peer.iq1","peer.iq5","peer.iq6"), each = 2), c("_out", "_in")),
                     "self.iq1", "self.iq5", "self.iq6",
                     "grade1", "grade2", "grade4")
    
    popcov.mat_c <- matrix(popcov.vec_c, nrow = 12, ncol = 12, 
                           dimnames = list(pop.names_c, pop.names_c))
    
    popcor.mat_c <- cov2cor(popcov.mat_c)
    popSD.vec_c <- sqrt(diag(popcov.mat_c))
    
    
    popcov.vec_d <- c(0.56852540197962,0.0216108061369262,0.0964086730430072,0.0468845156295315,
                      0.0842288129678118,0.0423410677086223,0.0216108061369262,0.56852540197962,
                      0.0468845156295315,0.0964086730430072,0.0423410677086223,0.0842288129678118,
                      0.0964086730430072,0.0468845156295315,0.601841259106791,0.0445854813706347,
                      0.263483539276719,0.0497165258128205,0.0468845156295315,0.0964086730430072,
                      0.0445854813706347,0.601841259106791,0.0497165258128205,0.263483539276719,
                      0.0842288129678118,0.0423410677086223,0.263483539276719,0.0497165258128205,
                      0.689620794901647,0.0383268533702981,0.0423410677086223,0.0842288129678118,
                      0.0497165258128205,0.263483539276719,0.0383268533702981,0.689620794901647)
    
    pop.names_d <- paste0(rep(c("peer.iq1","peer.iq5","peer.iq6"), each = 2), c("_ij", "_ji"))
    
    popcov.mat_d <- matrix(popcov.vec_d, nrow = 6, ncol = 6, 
                           dimnames = list(pop.names_d, pop.names_d))
    
    popcor.mat_d <- cov2cor(popcov.mat_d)
    popSD.vec_d <- sqrt(diag(popcov.mat_d))
    
    if (return_mats) {
      return(list(pop.cov_c = popcov.mat_c, pop.cor_c = popcor.mat_c, pop.SD_c = popSD.vec_c,
                  pop.cov_d = popcov.mat_d, pop.cor_d = popcor.mat_d, pop.SD_d = popSD.vec_d))
    } else {
      popcov.df_c <- as.data.frame(as.table(popcov.mat_c))
      popcov.df_c$par_names <- paste0(popcov.df_c$Var1, "~~", popcov.df_c$Var2)
      names(popcov.df_c)[names(popcov.df_c) == "Freq"] <- "pop.cov"
      popcov.df_c <- popcov.df_c[, c("par_names", "pop.cov")]
      
      popcor.df_c <- as.data.frame(as.table(popcor.mat_c))
      popcor.df_c$par_names <- paste0(popcor.df_c$Var1, "~~", popcor.df_c$Var2)
      names(popcor.df_c)[names(popcor.df_c) == "Freq"] <- "pop.cor"
      popcor.df_c <- popcor.df_c[, c("par_names", "pop.cor")]
      
      popSD.df_c <- as.data.frame(as.table(popSD.vec_c))
      popSD.df_c$par_names <- paste0(popSD.df_c$Var1, "~~", popSD.df_c$Var1)
      names(popSD.df_c)[names(popSD.df_c) == "Freq"] <- "pop.SD"
      popSD.df_c <- popSD.df_c[, c("par_names", "pop.SD")]
      
      popcov.df_d <- as.data.frame(as.table(popcov.mat_d))
      popcov.df_d$par_names <- paste0(popcov.df_d$Var1, "~~", popcov.df_d$Var2)
      names(popcov.df_d)[names(popcov.df_d) == "Freq"] <- "pop.cov"
      popcov.df_d <- popcov.df_d[, c("par_names", "pop.cov")]
      
      popcor.df_d <- as.data.frame(as.table(popcor.mat_d))
      popcor.df_d$par_names <- paste0(popcor.df_d$Var1, "~~", popcor.df_d$Var2)
      names(popcor.df_d)[names(popcor.df_d) == "Freq"] <- "pop.cor"
      popcor.df_d <- popcor.df_d[, c("par_names", "pop.cor")]
      
      popSD.df_d <- as.data.frame(as.table(popSD.vec_d))
      popSD.df_d$par_names <- paste0(popSD.df_d$Var1, "~~", popSD.df_d$Var1)
      names(popSD.df_d)[names(popSD.df_d) == "Freq"] <- "pop.SD"
      popSD.df_d <- popSD.df_d[, c("par_names", "pop.SD")]
    
      return(list(pop.cov = rbind(popcov.df_c, popcov.df_d),
                  pop.cor = rbind(popcov.df_c, popcov.df_d),
                  pop.SD = rbind(popSD.df_c, popSD.df_d)))
      
      }
} # TODO if bugs in the ANOVA_priors function() these valeus will have to be updated

# getSigma()
# getSigma(return_mats = FALSE)

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

# ANOVA_priors(rr.data = iq.data, case.data = covariate.data, 
#              IDout = "ego", IDin = "alter", IDgroup = "group",
#              default_prior = default_prior)


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