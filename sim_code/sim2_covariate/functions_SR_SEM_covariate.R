## Aditi M. Bhangale
## Last updated: 6 July 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Incorporating a level-specific covariate

## PREPARE STAGE1 `lavaan.srm` RESULTS FOR STAGE 2 + STAGE 2 `lavaan.srm`

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem/sim_code/sim2_covariate/")
# getwd()

# MCSampID = 1; n = 5; G = 3
# dat <- genGroups(1, 5, 3)
# rr.data = dat$rr.dat
# case.data = dat$covariate.dat
# rr.vars = c("peer.iq1","peer.iq5","peer.iq6")
# case.covs = c("self.iq1", "self.iq5", "self.iq6",
#               "grade1", "grade2", "grade4")
# IDout = "ego"; IDin = "alter"; IDgroup = "group"
# precision = 0.1
# iter = 10


#################################
# FUNCTIONS FOR SIMULATIONS ----
#################################

# function 0: generate level-specific (co)variance matrices----
getSigma <- function(return_mats = TRUE) {
    popcov.vec_c <- c(0.367178156002207,0.0282638172371223,0.232566187414696,
                      -0.00042303152370159,0.145404838183355,0.00256135404013844,
                      0.19650649315067,0.141654940954547,0.0832506443945378,
                      0.0179926118180312,0.0516975118928625,-0.00815220571888204,
                      0.0282638172371223,0.118625041049478,0.00436854440779768,
                      0.0994599593219533,0.00611922017653478,0.0940880968540291,
                      0.0772761368510366,0.0918692560924219,0.0809180729246042,
                      0.0826428880708735,0.0796935097239754,0.06850886275001,
                      0.232566187414696,0.00436854440779768,0.470529300000656,
                      -0.0119247565252019,0.309544649467953,-0.0430818983370867,
                      0.122650544584262,0.248888862834933,0.158902022440876,
                      0.00522078558607705,0.0818404324127288,-0.0155742030238611,
                      -0.00042303152370159,0.0994599593219533,-0.0119247565252019,
                      0.197907713215963,0.0107240331607381,0.178920190736849,
                      0.0914084704107476,0.0991671583731612,0.114559778372169,
                      0.192815472860378,0.168808200401105,0.140288023680691,
                      0.145404838183355,0.00611922017653478,0.309544649467953,
                      0.0107240331607381,0.419707598143518,0.0021001712597054,
                      0.0906398691171177,0.16010599825563,0.293721236120613,
                      0.0234391069284459,0.0570586359951145,0.0202568439692539,
                      0.00256135404013844,0.0940880968540291,-0.0430818983370867,
                      0.178920190736849,0.0021001712597054,0.189996325669525,
                      0.0820405997800706,0.0897308819246081,0.103607444797212,
                      0.179600581343282,0.168834134277284,0.144804682212415,
                      0.19650649315067,0.0772761368510366,0.122650544584262,
                      0.0914084704107476,0.0906398691171177,0.0820405997800706,
                      0.73648636471313,0.364561502339407,0.346920016798681,
                      0.102077086997001,0.0810390897699736,0.0658465393457167,
                      0.141654940954547,0.0918692560924219,0.248888862834933,
                      0.0991671583731612,0.16010599825563,0.0897308819246081,
                      0.364561502339407,0.650560216188368,0.40072530226852,
                      0.123389906869327,0.117031365329952,0.0602534410885451,
                      0.0832506443945378,0.0809180729246042,0.158902022440876,
                      0.114559778372169,0.293721236120613,0.103607444797212,
                      0.346920016798681,0.40072530226852,0.714613137577745,
                      0.135069195417365,0.1605460671607,0.134148244398022,
                      0.0179926118180312,0.0826428880708735,0.00522078558607705,
                      0.192815472860378,0.0234391069284459,0.179600581343282,
                      0.102077086997001,0.123389906869327,0.135069195417365,
                      0.508493900461932,0.287147572145552,0.208055233197341,
                      0.0516975118928625,0.0796935097239754,0.0818404324127288,
                      0.168808200401105,0.0570586359951145,0.168834134277284,
                      0.0810390897699736,0.117031365329952,0.1605460671607,
                      0.287147572145552,0.516170394363591,0.280945879685804,
                      -0.00815220571888204,0.06850886275001,-0.0155742030238611,
                      0.140288023680691,0.0202568439692539,0.144804682212415,
                      0.0658465393457167,0.0602534410885451,0.134148244398022,
                      0.208055233197341,0.280945879685804,0.475466332681621)
      
    pop.names_c <- c(paste0(rep(c("peer.iq1","peer.iq5","peer.iq6"), each = 2), c("_out", "_in")),
                     "self.iq1", "self.iq5", "self.iq6",
                     "grade1", "grade2", "grade4")
    
    popcov.mat_c <- matrix(popcov.vec_c, nrow = 12, ncol = 12, 
                           dimnames = list(pop.names_c, pop.names_c))
    
    popcor.mat_c <- cov2cor(popcov.mat_c)
    popSD.vec_c <- sqrt(diag(popcov.mat_c))
    
    
    popcov.vec_d <- c(0.56852540197962,0.0216108061369262,0.0964086730430072,
                      0.0468845156295315,0.0842288129678118,0.0423410677086223,
                      0.0216108061369262,0.56852540197962,0.0468845156295315,
                      0.0964086730430072,0.0423410677086223,0.0842288129678118,
                      0.0964086730430072,0.0468845156295315,0.601841259106791,
                      0.0445854813706347,0.263483539276719,0.0497165258128205,
                      0.0468845156295315,0.0964086730430072,0.0445854813706347,
                      0.601841259106791,0.0497165258128205,0.263483539276719,
                      0.0842288129678118,0.0423410677086223,0.263483539276719,
                      0.0497165258128205,0.689620794901647,0.0383268533702981,
                      0.0423410677086223,0.0842288129678118,0.0497165258128205,
                      0.263483539276719,0.0383268533702981,0.689620794901647)
    
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
      
      popsem.names_c <- c('grade4~grade2','grade4~peer.iq6_in','grade2~grade1',
                          'grade2~peer.iq5_in','grade1~peer.iq1_in',
                          'peer.iq6_in~peer.iq5_in','peer.iq6_in~grade2',
                          'peer.iq5_in~peer.iq1_in','peer.iq5_in~self.iq1',
                          'peer.iq5_in~grade1','self.iq6~self.iq5','self.iq6~peer.iq5_in',
                          'self.iq5~self.iq1','self.iq5~peer.iq1_in','self.iq5~grade1',
                          'grade4~~grade4','grade2~~grade2','grade1~~grade1',
                          'peer.iq6_in~~peer.iq6_in','peer.iq5_in~~peer.iq5_in',
                          'self.iq6~~self.iq6','self.iq5~~self.iq5','grade4~~self.iq6',
                          'peer.iq1_in~~peer.iq1_in','peer.iq1_in~~self.iq1',
                          'self.iq1~~self.iq1')
      
      popsem.ests_c <- c(0.405120883710625,0.403981609684353,0.378260262751128,
                         0.486959275968029,0.696672632598519,0.86641850652361,
                         0.0438158459260351,0.62232343063585,0.0387921375357394,
                         0.273939241174995,0.551687603067031,0.303990860771671,
                         0.438132160259862,0.429938689227192,0.0848295430112434,
                         0.300669089338049,0.324763370086136,0.450918912821721,
                         0.0275314576601663,0.0786205514346074,0.451860569400425,
                         0.440868881402277,0.032534353460116,0.118625139274703,
                         0.0772761368510366,0.736486364496517)
      
      pop.sem_c <- data.frame(par_names = popsem.names_c, 
                                pop.sem = popsem.ests_c)
    
      return(list(pop.cov = rbind(popcov.df_c, popcov.df_d),
                  pop.cor = rbind(popcov.df_c, popcov.df_d),
                  pop.SD = rbind(popSD.df_c, popSD.df_d),
                  pop.sem = pop.sem_c))
      
      }
} 

# getSigma()
# getSigma(return_mats = FALSE)

#----

# function 1: simulate data for a single group----
genData <- function(n) {
  popMats <- getSigma()
  
  pop.cov_c <- popMats$pop.cov_c
  pop.cov_d <- popMats$pop.cov_d
  
  pop.mu_c <- rep(0, 12)
  pop.mu_d <- rep(0, 6)
  
  library(mnormt) # for rmnorm()
  
  dat_c <- rmnorm(n = n, mean = pop.mu_c, varcov = pop.cov_c)
  rr.dat_d <- rmnorm(n = (n*(n-1))/2, mean = pop.mu_d, varcov = pop.cov_d)
  
  covariate.dat_c <- as.data.frame(cbind(ID = 1:n, dat_c[, c("self.iq1", "self.iq5", "self.iq6",
                               "grade1", "grade2", "grade4")]))
  
  rr.dat_c <- dat_c[, paste0(rep(c("peer.iq1","peer.iq5","peer.iq6"), 
                                 each = 2), c("_out", "_in"))]
  
  Y.dat <- NULL
  rr.names <- c("peer.iq1","peer.iq5","peer.iq6") 
  
  for (rr in rr.names) {
    ego.name <- paste0(rr, "_out")
    ego.mat <- matrix(rr.dat_c[,ego.name], nrow = n, ncol = n, byrow = F) # each row = new ego
    diag(ego.mat) <- NA
    
    alter.name <- paste0(rr, "_in")
    alter.mat <- matrix(rr.dat_c[,alter.name], nrow = n, ncol = n, byrow = T) # each col = new alter
    diag(alter.mat) <- NA
    
    R.mat <- matrix(NA, n, n)
    ij.name <- paste0(rr, "_ij")
    ji.name <- paste0(rr, "_ji")
    
    lt.indices <- which(lower.tri(R.mat, diag = F), arr.ind = T) # lower triangle indices
    
    R.mat[lt.indices] <- rr.dat_d[,ij.name]
    R.mat[lt.indices[,2:1]] <- rr.dat_d[,ji.name]
    
    Y.adj <- ego.mat + alter.mat + R.mat # adjacency matrix for Y 
    
    # convert to long format
    tempY.low <- cbind(lt.indices, Y.adj[lt.indices])
    colnames(tempY.low) <- c("ego", "alter", rr)
    
    tempY.up <- cbind(lt.indices[,2:1], Y.adj[lt.indices[,2:1]])
    colnames(tempY.up) <- c("ego", "alter", rr)
    
    tempY <- rbind(tempY.low, tempY.up)
  
  if (is.null(Y.dat)) {
    Y.dat <- tempY
  } else {
    Y.dat <- merge(Y.dat, tempY) # necessary to generate multigroup data (next function)
  }
  }
  return(list(rr.dat = Y.dat, covariate.dat = covariate.dat_c))
}

# genData(n = 5)

#----

# function 2: simulate data for multiple groups----
genGroups <- function(MCSampID, n, G) {
  # set the seed
  library(parallel)
  library(portableParallelSeeds)
  mySeeds <- seedCreator(nReps = 5000, streamsPerRep = 1, seed = 1974)
  setSeeds(projSeeds = mySeeds, run = MCSampID)
  
  rr.dat_multi <- lapply(1:G, function(g){
    rr.dat_single <- genData(n = n)$rr.dat
    
    rr.dat_single$ego <- rr.dat_single$ego + g*100
    rr.dat_single$alter <- rr.dat_single$alter + g*100
    
    cbind(group = g, rr.dat_single)
  })
  
  rr.dat_multi <- do.call("rbind", rr.dat_multi)
  
  covariate.dat_multi <- lapply(1:G, function(g) {
    covariate.dat_single <- genData(n = n)$covariate.dat
    
    covariate.dat_single$ID <- covariate.dat_single$ID + g*100
    
    cbind(group = g, covariate.dat_single)
  })
  
  covariate.dat_multi <- do.call("rbind", covariate.dat_multi)
  
  return(list(rr.dat = rr.dat_multi, covariate.dat = covariate.dat_multi))
}

# genGroups(MCSampID = 1, n = 5, G = 3)

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
prophetic_priors <- function(pop.corMat, pop.SDvec, rr.vars, case.covs, 
                             precision, default_prior) {
  
  priors <- default_prior
  
  # priors for SDs---t priors
  pop.SD_c <- pop.SDvec$pop.SD_c
  pop.SD_d <- pop.SDvec$pop.SD_d
  
  priors$rr_in_t$m <- pop.SD_c[paste0(rr.vars, "_in")]
  priors$rr_out_t$m <- pop.SD_c[paste0(rr.vars, "_out")]
  priors$rr_rel_t$m <- pop.SD_d[paste0(rr.vars, "_ij")]
  priors$case_cov_t$m <- pop.SD_c[case.covs]
  
  priors$rr_in_t$sd <- priors$rr_out_t$sd <- priors$rr_rel_t$sd <- rep(precision, length(rr.vars))
  priors$case_cov_t$sd <- rep(precision, length(case.covs))
  
  
  # priors for correlations---beta priors
  
  ## begin: case-level hyperpars----
  pop.cor_c <- pop.corMat$pop.cor_c
  
  targetVals_c <- as.list(pop.cor_c[lower.tri(pop.cor_c)])
  newHyperpars_c <- mapply(function(tc, ...) optim(par = c(1.5, 1.5), fn = minLoss, ...)$par, 
                           targetCorr = targetVals_c, accept.dev = precision,
                           method = "L-BFGS-B", lower = 0, SIMPLIFY = FALSE)
  
  # save the new (strongly informative) alpha and beta parameters
  alpha_c <- sapply(newHyperpars_c, "[[", 1)
  beta_c <- sapply(newHyperpars_c, "[[", 2)
  
  # populate alpha hyperpars
  priors$case_beta[lower.tri(priors$case_beta, diag = FALSE)] <- alpha_c
  
  # populate beta parameters
  betamat_c <- matrix(NA, 12, 12)
  betamat_c[lower.tri(betamat_c, diag = FALSE)] <- beta_c
  betamat_c <- t(betamat_c)
  priors$case_beta[upper.tri(priors$case_beta, diag = FALSE)] <- 
  betamat_c[upper.tri(betamat_c, diag = FALSE)]
  ## end: case-level hyperpars----
  
  ## begin: dyad-level hyperpars----
  pop.cor_d <- pop.corMat$pop.cor_d
  
  # saving the population correlation values
  ## DIAG == dyadic reciprocity -- [2,1]; [4,3]; [6,5]
  ## ABOVE = inter -- [3,2]; [5,2]; [5,4]
  ## BELOW = intra -- [4,2]; [6,2]; [6,4]
  targetVals_d <- list(dyad1 = pop.cor_d[2,1], 
                       dyad2 = pop.cor_d[4,3], 
                       dyad3 = pop.cor_d[6,5],
                       inter21 = pop.cor_d[3,2], 
                       inter31 = pop.cor_d[5,2], 
                       inter32 = pop.cor_d[5,4],
                       intra21 = pop.cor_d[4,2], 
                       intra31 = pop.cor_d[6,2], 
                       intra32 = pop.cor_d[6,4])
  
  newHyperpars_d <- mapply(function(tc, ...) optim(par = c(1.5, 1.5), fn = minLoss, ...)$par, 
                           targetCorr = targetVals_d, accept.dev = precision,
                           method = "L-BFGS-B", lower = 0, SIMPLIFY = FALSE)
  
  # new hyperpars for dyadic reciprocity, interpersonal corr, and intrapersonal corr
  dyad_alpha <- c(newHyperpars_d$dyad1[1], newHyperpars_d$dyad2[1], newHyperpars_d$dyad3[1])
  dyad_beta <- c(newHyperpars_d$dyad1[2], newHyperpars_d$dyad2[2], newHyperpars_d$dyad3[2])
  
  inter_alpha <- c(newHyperpars_d$inter21[1], newHyperpars_d$inter31[1], newHyperpars_d$inter32[1])
  inter_beta <- c(newHyperpars_d$inter21[2], newHyperpars_d$inter31[2], newHyperpars_d$inter32[2])
  
  intra_alpha <- c(newHyperpars_d$intra21[1], newHyperpars_d$intra31[1], newHyperpars_d$intra32[1])
  intra_beta <- c(newHyperpars_d$intra21[2], newHyperpars_d$intra31[2], newHyperpars_d$intra32[2])
  
  # replace in priors tables
  diag(priors$rr_beta_a) <- dyad_alpha
  diag(priors$rr_beta_b) <- dyad_beta # dyadic reciprocities on diagonal
  
  priors$rr_beta_a[upper.tri(priors$rr_beta_a)] <- inter_alpha
  priors$rr_beta_b[upper.tri(priors$rr_beta_b)] <- inter_beta # inter = above
  
  priors$rr_beta_a[lower.tri(priors$rr_beta_a)] <- intra_alpha
  priors$rr_beta_b[lower.tri(priors$rr_beta_b)] <- intra_beta # intra = below
  
  ## end: dyad-level hyperpars----
  
  return(priors)
}

# library(lavaan.srm)
# dat <- genGroups(1, 5, 3)
# popVals <- getSigma()
# prophetic_priors(pop.corMat = list(pop.cor_c = popVals$pop.cor_c, pop.cor_d = popVals$pop.cor_d),
#                  pop.SDvec = list(pop.SD_c = popVals$pop.SD_c, pop.SD_d = popVals$pop.SD_d),
#                  rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
#                  case.covs = c("self.iq1", "self.iq5", "self.iq6", "grade1", "grade2", "grade4"),
#                  precision = 0.1,
#                  default_prior = srm_priors(data = dat$rr.dat[,4:6], case_data = dat$covariate.dat[,3:8]))

#----

# function 5: ANOVA-based priors----
ANOVA_priors <- function(rr.data, case.data, rr.vars, case.covs, IDout, IDin, 
                         IDgroup, precision = 0.1, default_prior) {
  
  # to make compatible with TripleR output
  names(case.data)[names(case.data) == "group"] <- "group.id"
  names(case.data)[names(case.data) == "ID"] <- "id"
  
  # names(case.data) <- c("group.id", "id", case.covs) # to make compatible with TripleR output
  
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
    
    # construct correlation matrix at case level 
    # and covariance matrix at dyad level (to compute correlation matrix)
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
    
    for (cc in 1:length(case.covs)) {
      corMat_c[paste0(rr.vars[rr], "_out"), case.covs[cc]] <- corMat_c[case.covs[cc], paste0(rr.vars[rr], "_out")] <- 
        parCor(combi_dat[,case.covs[cc]], combi_dat[,paste0(rr.vars[rr], "_out")], combi_dat$group.id)$par.cor # correlation with outgoing effects
      corMat_c[paste0(rr.vars[rr], "_in"), case.covs[cc]] <- corMat_c[case.covs[cc], paste0(rr.vars[rr], "_in")] <- 
        parCor(combi_dat[,case.covs[cc]], combi_dat[,paste0(rr.vars[rr], "_in")], combi_dat$group.id)$par.cor # correlation with incoming effects
      
      # # SD priors for each covariate
      # priors$case_cov_t[case.covs[cc], "m"] <- sd(case.data[,case.covs[cc]])
      # priors$case_cov_t[case.covs[cc], "sd"] <- precision
      ## leave the SD priors for each covariate as the default (discussed this with T)
      
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
  
  corMat_d <- cov2cor(covMat_d) # convert dyad-level covariances to correlations
  
  for (cc in 2:length(case.covs)) { # insert correlations between case-level covariates
    for (bb in 1:(cc-1)) {
      corMat_c[case.covs[cc], case.covs[bb]] <- corMat_c[case.covs[bb], case.covs[cc]] <-
        parCor(combi_dat[,case.covs[cc]], combi_dat[,case.covs[bb]], combi_dat$group.id)$par.cor
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
  
  return(priors)
}

# library(lavaan.srm)
# dat <- genGroups(1, 5, 3)
# ANOVA_priors(rr.data = dat$rr.dat, case.data = dat$covariate.dat,
#              rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
#              case.covs = c("self.iq1", "self.iq5", "self.iq6",
#                            "grade1", "grade2", "grade4"),
#              IDout = "ego", IDin = "alter", IDgroup = "group", precision = 0.1,
#              default_prior = srm_priors(data = dat$rr.dat[,4:6], case_data = dat$covariate.dat[,3:8]))


#----

# function 6: set customised priors for MCMC stage----
set_priors <- function(rr.data, case.data, rr.vars, case.covs,
                       IDout = NULL, IDin = NULL, IDgroup = NULL, priorType, 
                       precision = NULL) {
  prior_env <- new.env()
  prior_env$default_prior <- srm_priors(data = rr.data[rr.vars], 
                                        case_data = case.data[case.covs])
  prior_env$popValues <- getSigma()
  prior_env$pop.corMat <- list(pop.cor_c = get("popValues", envir = prior_env)$pop.cor_c,
                               pop.cor_d = get("popValues", envir = prior_env)$pop.cor_d)
  prior_env$pop.SDvec <- list(pop.SD_c = get("popValues", envir = prior_env)$pop.SD_c, 
                              pop.SD_d = get("popValues", envir = prior_env)$pop.SD_d)
  
  if (priorType == "default") {
    srmPriors <- get("default_prior", prior_env)
  } else if (priorType == "prophetic") {
    srmPriors <- prophetic_priors(pop.corMat = get("pop.corMat", envir = prior_env),
                                  pop.SDvec = get("pop.SDvec", envir = prior_env),
                                  rr.vars = rr.vars, case.covs = case.covs, 
                                  precision = precision,
                                  default_prior = get("default_prior", prior_env))
  } else if (priorType == "ANOVA") {
    srmPriors <- ANOVA_priors(rr.data = rr.data, case.data = case.data, 
                              rr.vars = rr.vars, case.covs = case.covs,
                              IDout = IDout, IDin = IDin, IDgroup = IDgroup, 
                              precision = precision, 
                              default_prior = get("default_prior", prior_env))
  } else {
    stop("specify valid priorType (default, prophetic or ANOVA)")
  }
  return(srmPriors)
}

# set_priors(rr.data = rr.data, case.data = case.data,
#            rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
#            case.covs = c("self.iq1", "self.iq5", "self.iq6",
#                          "grade1", "grade2", "grade4"),
#            priorType = "default")
# set_priors(rr.data = rr.data, case.data = case.data,
#            rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
#            case.covs = c("self.iq1", "self.iq5", "self.iq6",
#                          "grade1", "grade2", "grade4"),
#            priorType = "prophetic", precision = 0.1)
# set_priors(rr.data = rr.data, case.data = case.data,
#            rr.vars = c("peer.iq1","peer.iq5","peer.iq6"),
#            case.covs = c("self.iq1", "self.iq5", "self.iq6",
#                          "grade1", "grade2", "grade4"),
#            IDout = "ego", IDin = "alter", IDgroup = "group",
#            priorType = "ANOVA", precision = 0.1)


#----

# function 7: stage 1 of lavaan.srm----
lavs1 <- function(MCSampID, n, G, rr.vars = c("peer.iq1","peer.iq5","peer.iq6"), 
                  case.covs = c("self.iq1", "self.iq5", "self.iq6", 
                                "grade1", "grade2", "grade4"), 
                  IDout = "ego", IDin = "alter", IDgroup = "group", 
                  priorType, precision = 0.1, iter = 2000, savefile = F) {
  library(lavaan.srm)
  library(coda) # for gelman.diag()
  library(rstan) # for As.mcmc.list()
  
  s1_env <- new.env()
  s1_env$dat <- genGroups(MCSampID = MCSampID, n = n, G = G)
  s1_env$rr.data <- get("dat", envir = s1_env)$rr.dat
  s1_env$covariate.data <- get("dat", envir = s1_env)$covariate.dat
  s1_env$MCMC_pars <- c(paste0("pSigma[", outer(1:12, 1:12, FUN = paste, sep = ",")[upper.tri(diag(12), diag = TRUE)], "]"), # case level
                       paste0("dSigma[", c(paste0(rep(1, each = 6), ",", 1:6), # relvar1; dyadcov11; intra/inter12; intra/inter13
                                           paste0(rep(3, each = 4), ",", 3:6), #relvar2; dyadcov22; intra/inter23
                                           paste0(rep(5, each = 2), ",", 5:6)), "]")) # relvar3; dyadcov33
  
  t0 <- Sys.time()
  
  if (priorType == "default") {
    
    s1_priors <- set_priors(rr.data = get("rr.data", envir = s1_env),
                            case.data = get("covariate.data", envir = s1_env),
                            rr.vars = rr.vars, case.covs = case.covs,
                            priorType = priorType)
    s1ests <- mvsrm(data = get("rr.data", envir = s1_env),
                    rr.vars = rr.vars, IDout = IDout, IDin = IDin, 
                    IDgroup = IDgroup, fixed.groups = T, 
                    case_data = get("covariate.data", envir = s1_env), init_r = 0.5, 
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                 probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
    
    if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
      iter = iter*2
      
      s1ests <- mvsrm(data = get("rr.data", envir = s1_env),
                      rr.vars = rr.vars, IDout = IDout, IDin = IDin, 
                      IDgroup = IDgroup, fixed.groups = T, 
                      case_data = get("covariate.data", envir = s1_env), init_r = 0.5, 
                      iter = iter, priors = s1_priors, seed = 1512, verbose = F)
      
      s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                   probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
      
    }
  } else if (priorType == "prophetic") {
    
    s1_priors <- set_priors(rr.data = get("rr.data", envir = s1_env),
                            case.data = get("covariate.data", envir = s1_env),
                            rr.vars = rr.vars, case.covs = case.covs,
                            priorType = priorType, precision = precision)
    s1ests <- mvsrm(data = get("rr.data", envir = s1_env),
                    rr.vars = rr.vars, IDout = IDout, IDin = IDin, 
                    IDgroup = IDgroup, fixed.groups = T, 
                    case_data = get("covariate.data", envir = s1_env), init_r = 0.5, 
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                 probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
    
    if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
      iter = iter*2
      
      s1ests <- mvsrm(data = get("rr.data", envir = s1_env),
                      rr.vars = rr.vars, IDout = IDout, IDin = IDin, 
                      IDgroup = IDgroup, fixed.groups = T, 
                      case_data = get("covariate.data", envir = s1_env), init_r = 0.5, 
                      iter = iter, priors = s1_priors, seed = 1512, verbose = F)
      
      s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                   probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
      
    }
    
  } else if (priorType == "ANOVA") {
    
    s1_priors <- set_priors(rr.data = get("rr.data", envir = s1_env),
                            case.data = get("covariate.data", envir = s1_env),
                            rr.vars = rr.vars, case.covs = case.covs,
                            IDout = IDout, IDin = IDin, IDgroup = IDgroup,
                            priorType = priorType, precision = precision)
    s1ests <- mvsrm(data = get("rr.data", envir = s1_env),
                    rr.vars = rr.vars, IDout = IDout, IDin = IDin, 
                    IDgroup = IDgroup, fixed.groups = T, 
                    case_data = get("covariate.data", envir = s1_env), init_r = 0.5, 
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                 probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
    
    if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
      iter = iter*2
      
      s1ests <- mvsrm(data = get("rr.data", envir = s1_env),
                      rr.vars = rr.vars, IDout = IDout, IDin = IDin, 
                      IDgroup = IDgroup, fixed.groups = T, 
                      case_data = get("covariate.data", envir = s1_env), init_r = 0.5, 
                      iter = iter, priors = s1_priors, seed = 1512, verbose = F)
      
      s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                   probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
      
    }
  }
  t1 <- Sys.time()
  
  # create outlier flag 
  attr(s1ests, "mPSRF") <- gelman.diag(As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env)), 
                                       autoburnin = FALSE)$mpsrf
  if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
    attr(s1ests, "Reff_outlier") <- TRUE # if ANY Rhats or effective sample sizes indicate non-convergence
  } else {
    attr(s1ests, "Reff_outlier") <- FALSE
  }
  
  if (attr(s1ests, "mPSRF") > 1.05) {
    attr(s1ests, "mPSRF_outlier") <- TRUE # if the overall mPSRF indicates non-convergence
  } else {
    attr(s1ests, "mPSRF_outlier") <- FALSE
  }
  
  # assign simulation conditions as attributes to s1ests to use for stage2
  attr(s1ests, "MCSampID") <- MCSampID
  attr(s1ests, "n") <- n
  attr(s1ests, "G") <- G
  attr(s1ests, "priorType") <- priorType
  attr(s1ests, "precision") <- ifelse(priorType == "default", "NA", precision)
  attr(s1ests, "iter") <- iter
  attr(s1ests, "RunTime") <- difftime(t1, t0, units = "mins")
  
  if (savefile) saveRDS(s1ests, paste0("s1covariate_ID", MCSampID, ".nG", G, ".n", n, "_", 
                                       priorType,  "_", ifelse(priorType == "default", "NA", precision), ".rds"))
  return(s1ests)
  
}

# foo <- lavs1(MCSampID = 1, n = 5, G = 3, priorType = "default", precision = NA, iter = 100)
# s1ests_prophetic <- lavs1(MCSampID = 1, n = 15, G = 10, priorType = "prophetic", precision = 0.1, iter = 500)
# baz <- lavs1(MCSampID = 1, n = 5, G = 3, priorType = "ANOVA", precision = 0.1, iter = 100)

#----

# function 8: stage 2 of lavaan.srm----
lavs2 <- function(s1ests, savefile = F) {
  library(lavaan.srm)
  
  if (attr(s1ests, "Reff_outlier") && attr(s1ests, "mPSRF_outlier")) {
    s2result <- NULL
  } else { 
  t0 <- Sys.time()  
  
  ## stage2
  mod_covariate <- '# regressions
    grade4 ~ grade2 + peer.iq6_in
    grade2 ~ grade1 + peer.iq5_in
    grade1 ~ peer.iq1_in
  
    peer.iq6_in ~ peer.iq5_in + grade2
    peer.iq5_in ~ peer.iq1_in + self.iq1 + grade1

    self.iq6 ~ self.iq5 + peer.iq5_in
    self.iq5 ~ self.iq1 + peer.iq1_in + grade1
  '
  
  fit_covariate <- sem.srm(mod_covariate, data = s1ests, component = "case",
                           posterior.est = "mean", 
                           test = c("satorra.bentler","scaled.shifted","browne.residual.adf"))
  
  #TODO talk to T about this: why does it by default correlate grade4 and self.iq6? i don't get it

  if (lavInspect(fit_covariate, "converged")) {
    ## save parameter estimates
    s2ests <- parameterEstimates(fit_covariate)
    s2ests$par_names <- paste0(s2ests$lhs, s2ests$op, s2ests$rhs)
    popVals <- getSigma(return_mats = F)$pop.sem 
    s2ests <- merge(s2ests, popVals, by = "par_names")
    s2ests$MCSampID <- attr(s1ests, "MCSampID")
    s2ests$n <- attr(s1ests, "n")
    s2ests$G <- attr(s1ests, "G")
    s2ests$condition <- paste0(attr(s1ests, "n"), "-", attr(s1ests, "G"))
    s2ests$analType <- "2SMLE"
    s2ests$s1priorType <- paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision"))
    s2ests$s1iter <- attr(s1ests, "iter")
    s2ests$s1mPSRF <- attr(s1ests, "mPSRF")
    s2ests$coverage <- s2ests$ci.lower < s2ests$pop & s2ests$pop < s2ests$ci.upper
    s2ests$bias <- s2ests$est - s2ests$pop
    s2ests$RB <- s2ests$bias / s2ests$pop
    s2ests$s1Reff_outlier <- attr(s1ests, "Reff_outlier")
    s2ests$s1mPSRF_outlier <- attr(s1ests, "mPSRF_outlier")
    
    s2ests <- s2ests[, c("MCSampID", "n", "G", "condition", "analType", "s1priorType", 
                         "s1iter", "s1mPSRF", "s1Reff_outlier", "s1mPSRF_outlier", 
                         "par_names", "pop.sem", "est", "se", "ci.lower", 
                         "ci.upper", "coverage", "bias", "RB")] # remove non-redundant columns and reorder
  
    ## evaluate model fit
    standard.LRT_covariate <- lavTestLRT(fit_covariate, method = "standard")
    standard.stat_covariate <- standard.LRT_covariate$`Chisq diff`[2]
    df_covariate <- standard.LRT_covariate$`Df diff`[2]
    standard.p_covariate <- standard.LRT_covariate$`Pr(>Chisq)`[2]
    
    adf.LRT_covariate <- lavTestLRT(fit_covariate, type = "browne.residual.adf")
    adf.stat_covariate <- adf.LRT_covariate$`Chisq diff`[2]
    adf.p_covariate <- adf.LRT_covariate$`Pr(>Chisq)`[2]
    
    sb.LRT_covariate <- lavTestLRT(fit_covariate, test = "satorra.bentler")
    sb.stat_covariate <- sb.LRT_covariate$`Chisq diff`[2]
    sb.p_covariate <- sb.LRT_covariate$`Pr(>Chisq)`[2]
    
    ss.LRT_covariate <- lavTestLRT(fit_covariate, test = "scaled.shifted")
    ss.stat_covariate <- ss.LRT_covariate$`Chisq diff`[2]
    ss.p_covariate <- ss.LRT_covariate$`Pr(>Chisq)`[2]
    
    N_covariate <- attr(s1ests, "nobs")["case"] + attr(s1ests, "nobs")["dyad"]
    yb.corrected.stat_covariate <- adf.stat_covariate / (1 + (adf.stat_covariate/N_covariate))
    yb.corrected.p_covariate <- pchisq(yb.corrected.stat_covariate, df = df_covariate, lower = F)
    
    yb.F.stat_covariate <- ((N_covariate - df_covariate)/(N_covariate * df_covariate))*adf.stat_covariate
    yb.F.df2_covariate <- N_covariate - df_covariate
    yb.F.p_covariate <- pf(yb.F.stat_covariate, df1 = df_covariate, df2 = yb.F.df2_covariate, lower.tail = F)
    
  } else {
      s2ests <- NULL
      
      standard.stat_covariate <- standard.p_covariate <- adf.stat_covariate <- adf.p_covariate <-
        sb.stat_covariate <- sb.p_covariate <- ss.stat_covariate <- ss.p_covariate <- 
        yb.corrected.stat_covariate <- yb.corrected.p_covariate <- yb.F.stat_covariate <- 
        yb.F.p_covariate <- NULL
      
      N_covariate <- attr(s1ests, "nobs")["covariate"]
      yb.F.df2_covariate <- N_covariate - df_covariate
    }
  
  standard <- c(fitStat.type = "standard", 
                covariate.stat = ifelse(!is.null(standard.stat_covariate), standard.stat_covariate, "NA"), 
                covariate.df = df_covariate, 
                covariate.p = ifelse(!is.null(standard.p_covariate), standard.p_covariate, "NA")) # Standard chi-square test
  
  adf <- c(fitStat.type = "browne.residual.adf", 
           covariate.stat = ifelse(!is.null(adf.stat_covariate), adf.stat_covariate, "NA"), 
           covariate.df = df_covariate, 
           covariate.p = ifelse(!is.null(adf.p_covariate), adf.p_covariate, "NA")) # Browne's residual-based ADF
  
  sb <- c(fitStat.type = "satorra.bentler.2001", 
          covariate.stat = ifelse(!is.null(sb.stat_covariate), sb.stat_covariate, "NA"), 
          covariate.df = df_covariate, 
          covariate.p = ifelse(!is.null(sb.p_covariate), sb.p_covariate, "NA")) # Satorra and Bentler's corrected statistic
  
  ss <- c(fitStat.type = "scaled.shifted", 
          covariate.stat = ifelse(!is.null(ss.stat_covariate), ss.stat_covariate, "NA"), 
          covariate.df = df_covariate, 
          covariate.p = ifelse(!is.null(ss.p_covariate), ss.p_covariate, "NA")) # Scaled-shifted statistic
  
  yb.corrected <- c(fitStat.type = "yuan.bentler.corrected.adf", 
                    covariate.stat = ifelse(!is.null(yb.corrected.stat_covariate), yb.corrected.stat_covariate, "NA"), 
                    covariate.df = df_covariate, 
                    covariate.p = ifelse(!is.null(yb.corrected.p_covariate), yb.corrected.p_covariate, "NA")) # Yuan and Bentler's small-sample correction for ADF
  
  yb.F <- c(fitStat.type = "yuan.bentler.F",
            covariate.stat = ifelse(!is.null(yb.F.stat_covariate), yb.F.stat_covariate, "NA"), 
            covariate.df = paste0(df_covariate, ",", yb.F.df2_covariate), 
            covariate.p = ifelse(!is.null(yb.F.p_covariate), yb.F.p_covariate, "NA"))
  
  s2mod <- data.frame(cbind(MCSampID = attr(s1ests, "MCSampID"),
                            n = attr(s1ests, "n"),
                            G = attr(s1ests, "G"),
                            condition = paste0(attr(s1ests, "n"), "-", attr(s1ests, "G")),
                            analType = "2SMLE",
                            s1priorType = paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision")),
                            s1iter = attr(s1ests, "iter"),
                            s1mPSRF = attr(s1ests, "mPSRF"), rbind(standard, adf, sb, ss, yb.corrected, yb.F)))
  
  t1 <- Sys.time()
  s2ests$RunTime <- difftime(t1, t0, units = "mins")
  
  s2result <- list(s2ests = s2ests, s2mod = s2mod)
  }
  if (savefile) saveRDS(s2result, file = paste0("s2covariate_ID", attr(s1ests, "MCSampID"),
                                                ".nG", attr(s1ests, "G"), ".n", 
                                                attr(s1ests, "n"), "_", attr(s1ests, "priorType"),
                                                "_", attr(s1ests, "precision"), ".rds")) #TODO check me
  return(s2result)
}

# lavs2(s1ests_prophetic)

#----

# function 9: create runsim files----
makeRunsim <- function(nSamps, n, G, priorType, precision = NULL, sim) {
  runsimfile <- paste0('## Aditi M. Bhangale
## Last updated:', Sys.Date(), 
                       
'\n# Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data

# covariate simulation

# runsim_s1covariate_',priorType, ifelse(!is.null(precision), paste0("_", precision), ""), '_n', n, '_G', G, '_',sim,'

source("functions_SR_SEM_covariate.R")

# specify conditions\n',
                       
                       priorType,'_grid <- expand.grid(MCSampID = 1:', nSamps, ', n = ', n, ', G = ', G, ', priorType = "', 
                       priorType, '",',
                       ifelse(!is.null(precision), paste0('precision = ', precision), 
                              paste0('precision = NA')), ',
                              stringsAsFactors = F)\n',
                       
                       priorType,'_grid$row_num <- 1:nrow(', priorType,'_grid)

# prepare parallel processing\n
library(doSNOW)

nClus <- 124
cl <- makeCluster(nClus)
registerDoSNOW(cl)

# run simulation\n',
                       
                       paste0('s1Result <- foreach(row_num = 1:nrow(',priorType,'_grid),
                    .packages = c("mnormt", "parallel", "portableParallelSeeds",
                    "lavaan.srm", "coda",
                    "modeest", "HDInterval", "rstan")) %dopar% {
                                    
out <- try(lavs1(MCSampID = ',priorType,'_grid[row_num, ]$MCSampID, 
                    n = ',priorType,'_grid[row_num, ]$n, G = ',priorType,'_grid[row_num, ]$G,
                    priorType = ',priorType,'_grid[row_num, ]$priorType, 
                    precision = ',priorType,'_grid[row_num, ]$precision), silent = T)
if(inherits(out, "try-error")) out <- NULL
                                    
return(out)
  }
         
# close cluster\n
stopCluster(cl)
         
saveRDS(s1Result, paste0("results_s1covariate_', priorType, 
                              ifelse(!is.null(precision), paste0("_", precision), ""), 
                              '_n', n, '_G', G, '_', sim,'_", Sys.Date(),".rds"))
         
         ')
  )
  
  cat(runsimfile, file = paste0("runsim_s1covariate_", priorType, ifelse(!is.null(precision), paste0("_", precision), ""), 
                                "_n", n, "_G", G, "_", sim, ".R"))
  invisible(NULL)
}

# makeRunsim(nSamps = 1, n = 5, G = 3, priorType = "prophetic", precision = 0.1, sim = "sim2")
# makeRunsim(nSamps = 5, n = 5, G = 3, priorType = "ANOVA", precision = 0.1, sim = "sim2")
# makeRunsim(nSamps = 5, n = 5, G = 3, priorType = "default", sim = "sim2")

#----

# function 10: create shell files----
makeShSnellius <- function(n, G, priorType, precision = NULL, sim, wallTime) {
  shell <- paste0('#!/bin/bash

#SBATCH -J ', paste0(priorType, ifelse(!is.null(precision), paste0("_", precision), ""), "_n", n, "_G", G, "_", sim),'
#SBATCH -e .', paste0(priorType, ifelse(!is.null(precision), paste0("_", precision), ""), "_n", n, "_G", G, "_", sim),'.SERR
#SBATCH -o .', paste0(priorType, ifelse(!is.null(precision), paste0("_", precision), ""), "_n", n, "_G", G, "_", sim),'.SOUT
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t ', wallTime,'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aditibhangale@gmail.com

cd "$TMPDIR"

module load 2023
module load R/4.3.2-gfbf-2023a
export MKL_NUM_THREADS=1

export R_LIBS=$HOME/rpackages:$R_LIBS

cp $HOME/SR-SEM/stage2sem/functions_SR_SEM_covariate.R "$TMPDIR"
cp $HOME/SR-SEM/stage2sem/', paste0("runsim_s1covariate_", priorType, ifelse(!is.null(precision), paste0("_", precision), ""), 
                                    "_n", n, "_G", G, "_", sim, ".R"),' "$TMPDIR"

Rscript --vanilla ', paste0("runsim_s1covariate_", priorType, ifelse(!is.null(precision), paste0("_", precision), ""), 
                            "_n", n, "_G", G, "_", sim, ".R"),'

cp "$TMPDIR"/*.rds $HOME/SR-SEM/stage2sem/' 
                  
  )
  cat(shell, file = paste0("shell_s1covariate_", priorType, ifelse(!is.null(precision), paste0("_", precision), ""), 
                           "_n", n, "_G", G, "_", sim, ".sh"))
  invisible(NULL)
}

# makeShSnellius(n = 5, G = 3, priorType = "prophetic", precision = 0.1, 
#                sim = "sim2", wallTime = "5-00:00:00")

#----

