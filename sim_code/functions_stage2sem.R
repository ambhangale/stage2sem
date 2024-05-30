## Aditi M. Bhangale
## Last updated: 30 May 2024

# " Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem")
# getwd()

## README: only run stage1 for software-default, prophetic, and EB-FIML priors
## (based on our findings from the stage1paper)

#TODO do we even want to run the n = 6 conditions for this paper? we already know they'll do poorly

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = 1, n = 5, G = 3)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision <- 0.1
# priorType = "prophetic"
# iter = 100

#################################
# FUNCTIONS FOR SIMULATIONS ----
#################################

# function 0: generate level-specific (co)variance matrices----

getSigma <- function(return_mats = TRUE){
  
  ## person effects
  
  Inames_c <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # person-level indicator names
  Fnames_c <- c("f1@A", "f1@P") # person-level factor names
  
  # person-level lambda
  LAM_c <- matrix(0, 6, 2)
  LAM_c[1,1] <- 1
  LAM_c[2,1] <- 1.2
  LAM_c[3,1] <- 0.7
  LAM_c[4,2] <- 1
  LAM_c[5,2] <- 0.6
  LAM_c[6,2] <- 0.6
  
  dimnames(LAM_c) <- list(Inames_c, Fnames_c)
  
  # person-level factor cov matrix
  PHI_c <- matrix(c(0.40, 0.05, 0.05, 0.20), 2, 2)
  
  dimnames(PHI_c) <- list(Fnames_c, Fnames_c)
  
  # person-level indicator residual cov matrix
  PSI_c <- diag(6)
  PSI_c[1,1] <- 0.2
  PSI_c[2,2] <- 0.2
  PSI_c[3,3] <- 0.2
  PSI_c[4,4] <- 0.1
  PSI_c[1,4] <- PSI_c[4,1] <- 0.05
  PSI_c[5,5] <- 0.1
  PSI_c[6,6] <- 0.1
  PSI_c[3,6] <- PSI_c[6,3] <- 0.03
  
  dimnames(PSI_c) <- list(Inames_c, Inames_c)
  
  ## structure dyad effects
  
  Inames_d <- c("V1@AP", "V1@PA", "V2@AP", "V2@PA", "V3@AP", "V3@PA") # dyad-level indicator names
  Fnames_d <- c("f1@AP", "f1@PA") # dyad-level factor names
  
  # dyad-level lambda
  LAM_d <- matrix(0, 6, 2)
  LAM_d[1,1] <- LAM_d[2,2] <- 1
  LAM_d[3,1] <- LAM_d[4,2] <- 0.8
  LAM_d[5,1] <- LAM_d[6,2] <- 1.4
  
  dimnames(LAM_d) <- list(Inames_d, Fnames_d)
  
  # dyad-level factor cov matrix
  PHI_d <- matrix(c(0.60, 0.15, 0.15, 0.60), 2, 2)
  
  dimnames(PHI_d) <- list(Fnames_d, Fnames_d)
  
  # dyad-level indicator residual cov matrix
  PSI_d <- diag(6)
  PSI_d[1,1] <- PSI_d[2,2] <- 0.3
  PSI_d[3,3] <- PSI_d[4,4] <- 0.5
  PSI_d[3,4] <- PSI_d[4,3] <- 0.1
  PSI_d[5,5] <- PSI_d[6,6] <- 0.4
  PSI_d[5,6] <- PSI_d[6,5] <- -0.2
  
  dimnames(PSI_d) <- list(Inames_d, Inames_d)
  
  # person-level sigma
  SIGMA_c <- LAM_c %*% PHI_c %*% t(LAM_c) + PSI_c
  # dyad-level sigma
  SIGMA_d <- LAM_d %*% PHI_d %*% t(LAM_d) + PSI_d
  
  ## pop values for satmod
  FsatNames_c <- gsub("V", "f", Inames_c) # person-level satmod factor names
  FsatNames_d <- gsub("V", "f", Inames_d) # dyad-level satmod factor names
  
  dimnames(SIGMA_c) <- list(FsatNames_c, FsatNames_c)
  dimnames(SIGMA_d) <- list(FsatNames_d, FsatNames_d)
  
  # person-level corr matrix
  R_c <- cov2cor(SIGMA_c)
  # dyad-level corr matrix
  R_d <- cov2cor(SIGMA_d)
  
  mat_list <- list(LAM_c=LAM_c, PHI_c=PHI_c, PSI_c=PSI_c,
                   LAM_d=LAM_d, PHI_d=PHI_d, PSI_d=PSI_d,
                   SIGMA_c = SIGMA_c, SIGMA_d = SIGMA_d, R_c = R_c, R_d = R_d)
  
  if (return_mats)  return(mat_list)
  
  # creating a dataframe of person-level population values
  pop_c.cov <- as.data.frame(as.table(SIGMA_c))
  colnames(pop_c.cov) <- c("row", "col", "pop.cov")
  pop_c.cov$par_names <- paste0(pop_c.cov$row, "~~", pop_c.cov$col)
  pop_c.cov <- pop_c.cov[, c("par_names", "pop.cov")]
  
  # person-level SDs
  mat_c.SD <- diag(sqrt(diag(SIGMA_c)))
  dimnames(mat_c.SD) <- list(FsatNames_c, FsatNames_c)
  pop_c.SD <- as.data.frame(as.table(mat_c.SD))
  pop_c.SD <- pop_c.SD[pop_c.SD$Freq != 0, ]
  colnames(pop_c.SD) <- c("row", "col", "pop.SD")
  pop_c.SD$par_names <- paste0(pop_c.SD$row, "~~", pop_c.SD$col)
  pop_c.SD <- pop_c.SD[, c("par_names", "pop.SD")]
  rownames(pop_c.SD) <- NULL
  
  # person-level correlations
  # mat_c.cor <- solve(mat_c.SD) %*% SIGMA_c %*% solve(mat_c.SD) # use R_c instead
  pop_c.cor <- as.data.frame(as.table(R_c))
  colnames(pop_c.cor) <- c("row", "col", "pop.cor")
  pop_c.cor$par_names <- paste0(pop_c.cor$row, "~~", pop_c.cor$col)
  pop_c.cor <- pop_c.cor[, c("par_names", "pop.cor")]
  
  # creating a dataframe of dyad-level population values
  pop_d.cov <- as.data.frame(as.table(SIGMA_d))
  colnames(pop_d.cov) <- c("row", "col", "pop.cov")
  pop_d.cov$par_names <-  paste0(pop_d.cov$row, "~~", pop_d.cov$col)
  pop_d.cov <- pop_d.cov[, c("par_names", "pop.cov")]
  
  # dyad-level SDs
  mat_d.SD <- diag(sqrt(diag(SIGMA_d)))
  dimnames(mat_d.SD) <- list(FsatNames_d, FsatNames_d)
  pop_d.SD <- as.data.frame(as.table(mat_d.SD))
  pop_d.SD <- pop_d.SD[pop_d.SD$Freq != 0, ]
  pop_d.SD <- pop_d.SD[!duplicated(pop_d.SD$Freq), ]
  colnames(pop_d.SD) <- c("row", "col", "pop.SD")
  pop_d.SD$par_names <- paste0(pop_d.SD$row, "~~", pop_d.SD$col)
  pop_d.SD <- pop_d.SD[, c("par_names", "pop.SD")]
  rownames(pop_d.SD) <- NULL
  
  # dyad-level correlations
  # mat_d.cor <- solve(mat_d.SD) %*% SIGMA_d %*% solve(mat_d.SD) # use R_d instead
  pop_d.cor <- as.data.frame(as.table(R_d))
  colnames(pop_d.cor) <- c("row", "col", "pop.cor")
  pop_d.cor$par_names <- paste0(pop_d.cor$row, "~~", pop_d.cor$col)
  pop_d.cor <- pop_d.cor[, c("par_names", "pop.cor")]
  
  return(list(pop.cov = rbind(pop_c.cov, pop_d.cov),
              pop.cor = rbind(pop_c.cor, pop_d.cor),
              pop.SD = rbind(pop_c.SD, pop_d.SD)))
}

# getSigma()
# getSigma(return_mats = F)

#----

# function 1: simulate data for a single group----

genData <- function(n) {
  
  SIGMA_mats <- getSigma()
  SIGMA_c <- SIGMA_mats$SIGMA_c # save person- and dyad-level pop.cov matrices
  SIGMA_d <- SIGMA_mats$SIGMA_d
  
  Inames_c <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # person-level indicator names
  Inames_d <- c("V1@AP", "V1@PA", "V2@AP", "V2@PA", "V3@AP", "V3@PA") # dyad-level indicator names
  
  # changing names to generate data per SRM component
  dimnames(SIGMA_c) <- list(Inames_c, Inames_c)
  dimnames(SIGMA_d) <- list(Inames_d, Inames_d)
  
  ##MEAN VECTOR-----
  mu <- rep(0, 6) # group-mean centered SRM component variables
  
  ##DATA GENERATION-----
  
  library(mnormt) # for rmnorm()
  
  dat_c <- rmnorm(n = n, mean = mu, varcov = SIGMA_c)
  dat_d <- rmnorm(n = (n*(n - 1))/2, mean = mu, varcov = SIGMA_d) 
  # Ndyads = n * (n - 1) -- and we already have separate AP and PA columns, 
  # so we need only (n*(n - 1))/2 rows in total
  
  Ydat <- NULL
  Vnames <- c("V1", "V2", "V3")
  
  for(v in Vnames) {
    
    # actor effects
    Aname <- paste0(v, "@A")
    Amat <- matrix(dat_c[, Aname], nrow = n, ncol = n, byrow = FALSE)
    diag(Amat) <- NA
    
    # partner effects
    Pname <- paste0(v, "@P")
    Pmat <- matrix(dat_c[, Pname], nrow = n, ncol = n, byrow = TRUE)
    diag(Pmat) <- NA
    
    # dyad-level effects
    Rmat <- matrix(NA, n, n)
    
    # _ij & _ji effects
    APname <- paste0(v, "@AP")
    PAname <- paste0(v, "@PA")
    
    # build the dyad-level matrix
    foo <- which(lower.tri(Rmat, diag = FALSE), arr.ind = TRUE)
    
    Rmat[foo] <- dat_d[, APname]
    Rmat[foo[, 2:1]] <- dat_d[, PAname]
    
    Y_adj <- Amat + Pmat + Rmat # adjacency matrix for Y
    
    # converting the data to long format
    
    tempY_low <- cbind(foo, Y_adj[foo])
    colnames(tempY_low) <- c("Actor", "Partner", v)
    
    tempY_up <- cbind(foo[, 2:1], Y_adj[foo[, 2:1]])
    colnames(tempY_up) <- c("Actor", "Partner", v)
    
    tempY <- rbind(tempY_low, tempY_up)
    
    # merge with previous variable's long-Y
    if (is.null(Ydat)) {
      Ydat <- tempY
    } else {
      Ydat <- merge(Ydat, tempY) # necessary to generate multigroup data (next function)
    }
  } # end for loop for data generation
  return(Ydat)
  
}

# genData(n = 5)

#----

# function 2: simulate data for multiple groups----

genGroups <- function(MCSampID, n, G) {
  
  # set the seed
  library(parallel)
  library(portableParallelSeeds)
  mySeeds <- seedCreator(nReps = 5000, streamsPerRep = 1, seed = 20505)
  
  setSeeds(projSeeds = mySeeds, run = MCSampID)
  
  gDat <- lapply(1:G, function(g){
    
    gen <- genData(n = n)
    
    gen$Actor <- gen$Actor + g*100
    gen$Partner <- gen$Partner + g*100
    
    cbind(Group = g, gen)
  })
  
  gDat <- do.call("rbind", gDat)
  
  return(gDat)
}

# genGroups(MCSampID = 1, n = 5, G = 3)

# ----

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

# function 4: prophetic priors----

prophetic_priors <- function(pop_corMat, pop_SDvec, precision, default_prior) {
  
  priors <- default_prior
  
  # for SDs --- t-priors
  popSD_c <- pop_SDvec$popSD_c
  popSD_d <- pop_SDvec$popSD_d
  
  priors$rr_in_t$m <- popSD_c[4:6] # incoming 
  priors$rr_out_t$m <- popSD_c[1:3] # outgoing
  priors$rr_rel_t$m <- popSD_d[c(1, 3, 5)] # relationship
  
  priors$rr_in_t$sd <- priors$rr_out_t$sd <- 
    priors$rr_rel_t$sd <- rep(precision, times = 3) # set SD of t-priors to 0.1
  
  # for correlations --- beta priors
  ## begin: case-level hyperpars----
  popCorr_c <- pop_corMat$popCorr_c
  Fnames_c <- c("f1@A", "f1@P", "f2@A", "f2@P", "f3@A", "f3@P") # reorder names
  popCorr_c <- popCorr_c[Fnames_c, Fnames_c]
  
  # save  population correlation values
  targetVals_c <- list(gen1 = popCorr_c[2,1], ee21 = popCorr_c[3,1], 
                       ae21 = popCorr_c[4,1], ee31 = popCorr_c[5,1], 
                       ae31 = popCorr_c[6,1], ea21 = popCorr_c[3,2], 
                       aa21 = popCorr_c[4,2], ea31 = popCorr_c[5,2], 
                       aa31 = popCorr_c[6,2], gen2 = popCorr_c[4,3], 
                       ee32 = popCorr_c[5,3], ae32 = popCorr_c[6,3],
                       ea32 = popCorr_c[5,4], aa32 = popCorr_c[6,4],
                       gen3 = popCorr_c[6,5])
  
  newHyperpars_c <- mapply(function(tc, ...) optim(par = c(1.5, 1.5), fn = minLoss, ...)$par, 
                           targetCorr = targetVals_c, accept.dev = precision,
                           method = "L-BFGS-B", lower = 0, SIMPLIFY = FALSE)
  
  # save the new (strongly informative) alpha and beta parameters
  alpha_c <- sapply(newHyperpars_c, "[[", 1)
  beta_c <- sapply(newHyperpars_c, "[[", 2)
  
  # populate alpha hyperpars
  priors$case_beta[lower.tri(priors$case_beta, diag = FALSE)] <- alpha_c
  
  # populate beta parameters
  betamat_c <- matrix(NA, 6, 6)
  betamat_c[lower.tri(betamat_c, diag = FALSE)] <- beta_c
  betamat_c <- t(betamat_c)
  priors$case_beta[upper.tri(priors$case_beta, diag = FALSE)] <- 
    betamat_c[upper.tri(betamat_c, diag = FALSE)]
  ## end: case-level hyperpars----
  
  ## begin: dyad-level hyperpars----
  popCorr_d <- pop_corMat$popCorr_d
  
  # saving the population correlation values
  ## DIAG == dyadic reciprocity -- [2,1]; [4,3]; [6,5]
  ## ABOVE = inter -- [3,2]; [5,2]; [5,4]
  ## BELOW = intra -- [4,2]; [6,2]; [6,4]
  targetVals_d <- list(dyad1 = popCorr_d[2,1], 
                       dyad2 = popCorr_d[4,3], 
                       dyad3 = popCorr_d[6,5],
                       inter21 = popCorr_d[3,2], 
                       inter31 = popCorr_d[5,2], 
                       inter32 = popCorr_d[5,4],
                       intra21 = popCorr_d[4,2], 
                       intra31 = popCorr_d[6,2], 
                       intra32 = popCorr_d[6,4])
  
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

#----

# function 5: FIML-based priors----

FIML_priors <- function(data, rr.vars, IDout, IDin, IDgroup, precision = NULL, 
                        default_prior) {
  
  rr.data <- data
  priors <- default_prior # default priors
  
  library(srm)
  library(car)
    
    for (rr in 1:length(rr.vars)) {
      # uvSRM
      uv.syn <- paste0("%Person\n",
                       "f_", rr.vars[rr], "@A=~1*", rr.vars[rr], "@A\n", 
                       "f_", rr.vars[rr], "@P=~1*", rr.vars[rr], "@P\n", # F loadings
                       
                       rr.vars[rr], "@A~~0*", rr.vars[rr], "@A+0*", rr.vars[rr], "@P\n",
                       rr.vars[rr], "@P~~0*", rr.vars[rr], "@P\n", # residual (co)v
                       
                       "f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@A+f_", rr.vars[rr], "@P\n",
                       "f_", rr.vars[rr], "@P~~f_", rr.vars[rr], "@P\n", # A/P var + gen cov
                       
                       "%Dyad\n",
                       "f_", rr.vars[rr], "@AP=~1*", rr.vars[rr], "@AP\n",
                       "f_", rr.vars[rr], "@PA=~1*", rr.vars[rr], "@PA\n", # F loadings
                       
                       rr.vars[rr], "@AP~~0*", rr.vars[rr], "@AP+0*", rr.vars[rr], "@PA\n",
                       rr.vars[rr], "@PA~~0*", rr.vars[rr], "@PA\n", # residual (co)v
                       
                       "f_", rr.vars[rr], "@AP~~var", rr.vars[rr], "*f_", rr.vars[rr], "@AP+f_", rr.vars[rr], "@PA\n",
                       "f_", rr.vars[rr], "@PA~~var", rr.vars[rr], "*f_", rr.vars[rr], "@PA") # AP/PA var + dyadic cov
      
      uv.fit <- srm(model.syntax = uv.syn, data = rr.data, rrgroup_name = IDgroup,
                    person_names = c(IDout, IDin), 
                    fixed.groups = FALSE,
                    verbose = FALSE)
      
      MIuvcov_c <- uv.fit$sigma$U[[1]]
      MIuvdimNames_c <- paste0(rep(rr.vars[rr], each = 2), c("@A", "@P")) 
      dimnames(MIuvcov_c) <- list(MIuvdimNames_c, MIuvdimNames_c)
      if (MIuvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] >= 0) { # if actor var > 0
        AvarDelta_c <- deltaMethod(uv.fit, g. = paste0("sqrt(`f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`)"))
        priors$rr_out_t[rr.vars[rr], "m"] <- AvarDelta_c$Estimate
        priors$rr_out_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, AvarDelta_c$SE)
      } else {
        priors$rr_out_t[rr.vars[rr], "m"] <- 0
        priors$rr_out_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision,
                                                     uv.fit$parm.table[uv.fit$parm.table$par_names == 
                                                                         paste0("f_", rr.vars[rr],
                                                                                "@A~~f_", rr.vars[rr], "@A"), "se"])
      }
      if (MIuvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] >= 0) { # if partner var > 0
        PvarDelta_c <- deltaMethod(uv.fit, g. = paste0("sqrt(`f_", rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`)"))
        priors$rr_in_t[rr.vars[rr], "m"] <- PvarDelta_c$Estimate
        priors$rr_in_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, PvarDelta_c$SE)
      } else {
        priors$rr_in_t[rr.vars[rr], "m"] <- 0
        priors$rr_in_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, 
                                                    uv.fit$parm.table[uv.fit$parm.table$par_names == 
                                                                        paste0("f_", rr.vars[rr],
                                                                               "@P~~f_", rr.vars[rr], "@P"), "se"])
      }
      if (MIuvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] > 0 && 
          MIuvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] > 0) { # if both actor & partner var >0
        
        ## generalised cor
        uvcorDelta_c <- deltaMethod(uv.fit, g. = paste0("`f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@P`/(sqrt(`f_", 
                                                        rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`*`f_", 
                                                        rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`))"))
        MIuvcor_c <- uvcorDelta_c$Estimate
        MIuvcorSE_c <- uvcorDelta_c$SE
        MIuvcor_c <- ifelse(MIuvcor_c > 0.9, 0.9, ifelse(MIuvcor_c < -0.9, -0.9, MIuvcor_c))
        genR_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = MIuvcor_c, 
                                accept.dev = ifelse(!is.null(precision), precision, MIuvcorSE_c),
                                method = "L-BFGS-B", lower = 0)$par
        priors$case_beta[paste0(rr.vars[rr],"_in"),paste0(rr.vars[rr],"_out")] <- genR_hyperpars[1] # lower triangle
        priors$case_beta[paste0(rr.vars[rr],"_out"),paste0(rr.vars[rr],"_in")] <- genR_hyperpars[2] # upper triangle
        
      } else {
        priors$case_beta[paste0(rr.vars[rr],"_in"),
                         paste0(rr.vars[rr],"_out")] <- priors$case_beta[paste0(rr.vars[rr],"_out"),
                                                                         paste0(rr.vars[rr],"_in")] <- 1.5
      }
      
      MIuvcov_d <- uv.fit$sigma$D[[1]]
      MIuvdimNames_d <- paste0(rep(rr.vars[rr], each = 2), c("@AP", "@PA"))
      dimnames(MIuvcov_d) <- list(MIuvdimNames_d, MIuvdimNames_d)
      if (MIuvcov_d[paste0(rr.vars[rr], "@AP"), paste0(rr.vars[rr], "@AP")] > 0) { # if rel var > 0
        uvcorDelta_d <- deltaMethod(uv.fit, g. = paste0("`f_", rr.vars[rr], "@AP~~f_", rr.vars[rr], "@PA`/`f_", 
                                                        rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`"))
        MIuvcor_d <- uvcorDelta_d$Estimate
        MIuvcorSE_d <- uvcorDelta_d$SE
        MIuvcor_d <- ifelse(MIuvcor_d > 0.9, 0.9, ifelse(MIuvcor_d < -0.9, -0.9, MIuvcor_d))
        dyadicR_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                   targetCorr = MIuvcor_d, 
                                   accept.dev = ifelse(!is.null(precision), precision, MIuvcorSE_d),
                                   method = "L-BFGS-B", lower = 0)$par 
        priors$rr_beta_a[rr,rr] <- dyadicR_hyperpars[1]
        priors$rr_beta_b[rr,rr] <- dyadicR_hyperpars[2]
        
        relvarDelta_c <- deltaMethod(uv.fit, g. = paste0("sqrt(`f_", rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`)"))
        priors$rr_rel_t[rr.vars[rr], "m"] <- relvarDelta_c$Estimate
        priors$rr_rel_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision, relvarDelta_c$SE)
      } else {
        priors$rr_beta_a[rr,rr] <- priors$rr_beta_b[rr,rr] <- 1.5
        
        priors$rr_rel_t[rr.vars[rr], "m"] <- 0
        priors$rr_rel_t[rr.vars[rr], "sd"] <- ifelse(!is.null(precision), precision,
                                                     uv.fit$parm.table[uv.fit$parm.table$par_names == 
                                                                         paste0("f_", rr.vars[rr],
                                                                                "@AP~~f_", rr.vars[rr], "@AP"), "se"])
      }
      
      if (rr > 1L) for (kk in 1:(rr-1)) {
        # bvSRM
        bv.syn <- paste0("%Person\n",
                         "f_", rr.vars[rr], "@A=~1*", rr.vars[rr], "@A\n",
                         "f_", rr.vars[kk], "@A=~1*", rr.vars[kk], "@A\n",
                         
                         "f_", rr.vars[rr], "@P=~1*", rr.vars[rr], "@P\n",
                         "f_", rr.vars[kk], "@P=~1*", rr.vars[kk], "@P\n",
                         
                         rr.vars[rr], "@A~~0*", rr.vars[rr], "@A+0*", rr.vars[rr], "@P\n",
                         rr.vars[rr], "@P~~0*", rr.vars[rr], "@P\n",
                         rr.vars[kk], "@A~~0*", rr.vars[kk], "@A+0*", rr.vars[kk], "@P\n",
                         rr.vars[kk], "@P~~0*", rr.vars[kk], "@P\n",
                         
                         "f_", rr.vars[rr], "@A~~f_", rr.vars[rr], "@A+f_", rr.vars[rr], "@P+f_", rr.vars[kk], "@A+f_", rr.vars[kk], "@P\n",
                         "f_", rr.vars[kk], "@A~~f_", rr.vars[kk], "@A+f_", rr.vars[kk], "@P+f_", rr.vars[rr], "@P\n",
                         "f_", rr.vars[rr], "@P~~f_", rr.vars[rr], "@P+f_", rr.vars[kk], "@P\n",
                         "f_", rr.vars[kk], "@P~~f_", rr.vars[kk], "@P\n",
                         
                         "%Dyad\n",
                         "f_", rr.vars[rr], "@AP=~1*", rr.vars[rr], "@AP\n",
                         "f_", rr.vars[kk], "@AP=~1*", rr.vars[kk], "@AP\n",
                         
                         "f_", rr.vars[rr], "@PA=~1*", rr.vars[rr], "@PA\n",
                         "f_", rr.vars[kk], "@PA=~1*", rr.vars[kk], "@PA\n",
                         
                         rr.vars[rr], "@AP~~0*", rr.vars[rr], "@AP+0*", rr.vars[rr], "@PA\n",
                         rr.vars[rr], "@PA~~0*", rr.vars[rr], "@PA\n",
                         rr.vars[kk], "@AP~~0*", rr.vars[kk], "@AP+0*", rr.vars[kk], "@PA\n",
                         rr.vars[kk], "@PA~~0*", rr.vars[kk], "@PA\n",
                         
                         "f_", rr.vars[rr], "@AP~~var", rr.vars[rr], "*f_", rr.vars[rr], "@AP+dyad", 
                         rr.vars[rr], "*f_", rr.vars[rr], "@PA+intra", rr.vars[rr], 
                         rr.vars[kk], "*f_", rr.vars[kk], "@AP+inter", rr.vars[rr], 
                         rr.vars[kk], "*f_", rr.vars[kk], "@PA\n",
                         "f_", rr.vars[kk], "@AP~~var", rr.vars[kk],"*f_", rr.vars[kk], 
                         "@AP+dyad", rr.vars[kk],"*f_", rr.vars[kk], "@PA+inter", 
                         rr.vars[rr], rr.vars[kk], "*f_", rr.vars[rr], "@PA\n",
                         "f_", rr.vars[rr], "@PA~~var", rr.vars[rr], "*f_", rr.vars[rr], 
                         "@PA+intra", rr.vars[rr], rr.vars[kk], "*f_", rr.vars[kk], "@PA\n",
                         "f_", rr.vars[kk], "@PA~~var", rr.vars[kk],"*f_", rr.vars[kk], "@PA")
        
        bv.fit <- srm(model.syntax = bv.syn, data = rr.data, rrgroup_name = IDgroup,
                      person_names = c(IDout, IDin), 
                      fixed.groups = FALSE,
                      verbose = FALSE)
        
        # case-level hyperparameters
        MIbvcov_c <- bv.fit$sigma$U[[1]]
        MIbvdimNames_c <- c(paste0(c(rr.vars[kk], rr.vars[rr]), "@A"),
                            paste0(c(rr.vars[kk], rr.vars[rr]), "@P"))
        dimnames(MIbvcov_c) <- list(MIbvdimNames_c, MIbvdimNames_c)
        ## AA
        if (MIbvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@A"), paste0(rr.vars[kk], "@A")] > 0) {
          bvAADelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@A~~f_", rr.vars[kk], "@A`/(sqrt(`f_", 
                                                         rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`*`f_", 
                                                         rr.vars[kk], "@A~~f_", rr.vars[kk], "@A`))"))
          AA_c <- bvAADelta_c$Estimate
          AASE_c <- bvAADelta_c$SE
          AA_c <- ifelse(AA_c > 0.9, 0.9, ifelse(AA_c < -0.9, -0.9, AA_c))
          AA_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = AA_c, 
                                accept.dev = ifelse(!is.null(precision), precision, AASE_c),
                                method = "L-BFGS-B", lower = 0)$par ## AA corr
          priors$case_beta[paste0(rr.vars[rr], "_out"), paste0(rr.vars[kk], "_out")] <- AA_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_out"), paste0(rr.vars[rr], "_out")] <- AA_hyperpars[2] # upper triangle
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_out"), 
                           paste0(rr.vars[kk], "_out")] <- priors$case_beta[paste0(rr.vars[kk], "_out"), 
                                                                            paste0(rr.vars[rr], "_out")] <- 1.5
        }
        ## AP
        if (MIbvcov_c[paste0(rr.vars[rr], "@A"), paste0(rr.vars[rr], "@A")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@P"), paste0(rr.vars[kk], "@P")] > 0) {
          bvAPDelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@A~~f_", rr.vars[kk], "@P`/(sqrt(`f_", 
                                                         rr.vars[rr], "@A~~f_", rr.vars[rr], "@A`*`f_", 
                                                         rr.vars[kk], "@P~~f_", rr.vars[kk], "@P`))"))
          AP_c <- bvAPDelta_c$Estimate
          APSE_c <- bvAPDelta_c$SE
          AP_c <- ifelse(AP_c > 0.9, 0.9, ifelse(AP_c < -0.9, -0.9, AP_c))
          AP_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = AP_c, 
                                accept.dev = ifelse(!is.null(precision), precision, APSE_c),
                                method = "L-BFGS-B", lower = 0)$par ## AP corr
          priors$case_beta[paste0(rr.vars[rr], "_out"), paste0(rr.vars[kk], "_in")] <- AP_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_in"), paste0(rr.vars[rr], "_out")] <- AP_hyperpars[2] # upper triangle
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_out"), 
                           paste0(rr.vars[kk], "_in")] <- priors$case_beta[paste0(rr.vars[kk], "_in"), 
                                                                           paste0(rr.vars[rr], "_out")] <- 1.5
        }
        ## PA
        if (MIbvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@A"), paste0(rr.vars[kk], "@A")] > 0) {
          bvPADelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[kk], "@A~~f_", rr.vars[rr], "@P`/(sqrt(`f_", 
                                                         rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`*`f_", 
                                                         rr.vars[kk], "@A~~f_", rr.vars[kk], "@A`))"))
          PA_c <- bvPADelta_c$Estimate
          PASE_c <- bvPADelta_c$SE
          PA_c <- ifelse(PA_c > 0.9, 0.9, ifelse(PA_c < -0.9, -0.9, PA_c))
          PA_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = PA_c, 
                                accept.dev = ifelse(!is.null(precision), precision, PASE_c),
                                method = "L-BFGS-B", lower = 0)$par ## PA corr
          priors$case_beta[paste0(rr.vars[rr], "_in"), paste0(rr.vars[kk], "_out")] <- PA_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_out"), paste0(rr.vars[rr], "_in")] <- PA_hyperpars[2] # upper triangle
          
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_in"), 
                           paste0(rr.vars[kk], "_out")] <- priors$case_beta[paste0(rr.vars[kk], "_out"), 
                                                                            paste0(rr.vars[rr], "_in")] <- 1.5
        }
        ## PP
        if (MIbvcov_c[paste0(rr.vars[rr], "@P"), paste0(rr.vars[rr], "@P")] > 0 &&
            MIbvcov_c[paste0(rr.vars[kk], "@P"), paste0(rr.vars[kk], "@P")] > 0) {
          
          bvPPDelta_c <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@P~~f_", rr.vars[kk], "@P`/(sqrt(`f_", 
                                                         rr.vars[rr], "@P~~f_", rr.vars[rr], "@P`*`f_", 
                                                         rr.vars[kk], "@P~~f_", rr.vars[kk], "@P`))"))
          PP_c <- bvPPDelta_c$Estimate
          PPSE_c <- bvPPDelta_c$SE
          PP_c <- ifelse(PP_c > 0.9, 0.9, ifelse(PP_c < -0.9, -0.9, PP_c))
          PP_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                targetCorr = PP_c, 
                                accept.dev = ifelse(!is.null(precision), precision, PPSE_c),
                                method = "L-BFGS-B", lower = 0)$par ## PA corr
          priors$case_beta[paste0(rr.vars[rr], "_in"), paste0(rr.vars[kk], "_in")] <- PP_hyperpars[1] # lower triangle
          priors$case_beta[paste0(rr.vars[kk], "_in"), paste0(rr.vars[rr], "_in")] <- PP_hyperpars[2] # upper triangle
        } else {
          priors$case_beta[paste0(rr.vars[rr], "_in"), 
                           paste0(rr.vars[kk], "_in")] <- priors$case_beta[paste0(rr.vars[kk], "_in"), 
                                                                           paste0(rr.vars[rr], "_in")] <- 1.5
        }
        
        # dyad-level hyperparameters
        MIbvcov_d <- bv.fit$sigma$D[[1]]
        MIbvdimNames_d <- c(paste0(rep(rr.vars[kk], times = 2), c("@AP", "@PA")), 
                            paste0(rep(rr.vars[rr], times = 2), c("@AP", "@PA"))) #FIXME maybe just directly specify _ij/ji here?
        dimnames(MIbvcov_d) <- list(MIbvdimNames_d, MIbvdimNames_d)
        if (MIbvcov_d[paste0(rr.vars[rr], "@AP"), paste0(rr.vars[rr], "@AP")] > 0 && 
            MIbvcov_d[paste0(rr.vars[kk], "@AP"), paste0(rr.vars[kk], "@AP")] > 0) {
          ### intra-- below diagonal
          bvintraDelta_d <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[rr], "@AP~~f_", rr.vars[kk], "@AP`/(sqrt(`f_", 
                                                            rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`*`f_", 
                                                            rr.vars[kk], "@AP~~f_", rr.vars[kk], "@AP`))"))
          intra_d <- bvintraDelta_d$Estimate
          intraSE_d <- bvintraDelta_d$SE
          intra_d <- ifelse(intra_d > 0.9, 0.9, ifelse(intra_d < -0.9, -0.9, intra_d))
          intra_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                   targetCorr = intra_d,
                                   accept.dev = ifelse(!is.null(precision), precision, intraSE_d), 
                                   method = "L-BFGS-B", lower = 0)$par
          priors$rr_beta_a[rr,kk] <- intra_hyperpars[1]
          priors$rr_beta_b[rr,kk] <- intra_hyperpars[2]
          
          ### inter-- above diagonal
          bvinterDelta_d <- deltaMethod(bv.fit, g. = paste0("`f_", rr.vars[kk], "@AP~~f_", rr.vars[rr], "@PA`/(sqrt(`f_", 
                                                            rr.vars[rr], "@AP~~f_", rr.vars[rr], "@AP`*`f_", 
                                                            rr.vars[kk], "@AP~~f_", rr.vars[kk], "@AP`))"))
          inter_d <- bvinterDelta_d$Estimate
          interSE_d <- bvinterDelta_d$SE
          inter_d <- ifelse(inter_d > 0.9, 0.9, ifelse(inter_d < -0.9, -0.9, inter_d))
          inter_hyperpars <- optim(par = c(1.5, 1.5), fn = minLoss,
                                   targetCorr = inter_d,
                                   accept.dev = ifelse(!is.null(precision), precision, interSE_d), 
                                   method = "L-BFGS-B", lower = 0)$par
          priors$rr_beta_a[kk,rr] <- inter_hyperpars[1]
          priors$rr_beta_b[kk,rr] <- inter_hyperpars[2]
        } else {
          priors$rr_beta_a[rr,kk] <- priors$rr_beta_b[rr,kk] <- 
            priors$rr_beta_a[kk,rr] <- priors$rr_beta_b[kk,rr] <- 1.5
        }
      }
    }
  priors
}

#----

# function 6: set customised priors for MCMC stage----

set_priors <- function(data, rr.vars, IDout, IDin, IDgroup, priorType, targetCorr,
                       precision) {
  prior_env <- new.env()
  prior_env$default_prior <- srm_priors(data = data[rr.vars]) # default MCMC priors (diffuse priors)
  prior_env$pop_corMat <- list(popCorr_c = getSigma()$R_c,
                               popCorr_d = getSigma()$R_d)
  prior_env$pop_SDvec <- list(popSD_c = sqrt(diag(getSigma()$SIGMA_c)),
                              popSD_d = sqrt(diag(getSigma()$SIGMA_d)))
  
  if (priorType == "default") { # default (diffuse) priors
    srmPriors <- get("default_prior", envir = prior_env)
  } else if (priorType == "prophetic") { # prophetic priors
    srmPriors <- prophetic_priors(pop_corMat = get("pop_corMat", envir = prior_env), 
                                  pop_SDvec = get("pop_SDvec", envir = prior_env), 
                                  precision = precision, 
                                  default_prior = get("default_prior", envir = prior_env))
    
  } else if (priorType == "FIML") { # FIML-based priors (`srm`)
    srmPriors <- FIML_priors(data = data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                             IDgroup = IDgroup, precision = precision,
                             default_prior = get("default_prior", envir = prior_env))
  }
  return(srmPriors)
}

#----

# function 7: stage1 in `lavaan.srm`----

lavs1 <- function(MCSampID, n, G, rr.vars = c("V1", "V2", "V3"), IDout = "Actor", 
                  IDin = "Partner", IDgroup = "Group", priorType, 
                  precision = 0.1, iter = 2000, savefile = FALSE) {
  library(lavaan.srm)
  library(coda) # for gelman.diag()
  library(rstan) # for As.mcmc.list()
  
  s1_env <- new.env()
  s1_env$dat <- genGroups(MCSampID = MCSampID, n = n, G = G)
  s1_env$MCMC_pars <- c(paste0("pSigma[", outer(1:6, 1:6, FUN = paste, sep = ",")[upper.tri(diag(6), diag = TRUE)], "]"), # case level
                                     paste0("dSigma[", c(paste0(rep(1, each = 6), ",", 1:6), # relvar1; dyadcov11; intra/inter12; intra/inter13
                                                         paste0(rep(3, each = 4), ",", 3:6), #relvar2; dyadcov22; intra/inter23
                                                         paste0(rep(5, each = 2), ",", 5:6)), "]")) # relvar3; dyadcov33
  t0 <- Sys.time()
  
  if (priorType == "default") { # default
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, priorType = priorType)
    
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                 probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
    
    if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
      iter <- iter*2
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    }
  } else if (priorType == "prophetic") { # prophetic
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, priorType = priorType,
                            precision = precision)
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
    s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                 probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
    
    if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
      iter <- iter*2
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    }
  } else if (priorType == "FIML") { # FIML
    rr.data <- get("dat", envir = s1_env)
    
    s1_priors <- set_priors(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                            IDgroup = IDgroup, priorType = priorType, precision = precision)
    s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                    IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                    iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    
   s1long <- data.frame(summary(s1ests, as.stanfit = TRUE, 
                                 probs = NULL)$summary)[get("MCMC_pars", envir = s1_env),]
    
    if (any(s1long$n_eff < 100, na.rm = TRUE) | any(s1long$Rhat > 1.05, na.rm = TRUE)) {
      iter <- iter*2
      s1ests <- mvsrm(data = rr.data, rr.vars = rr.vars, IDout = IDout, IDin = IDin,
                      IDgroup = IDgroup, fixed.groups = T, init_r = 0.5,
                      iter = iter, priors = s1_priors, seed = 1512, verbose = F)
    }
  }
  t1 <- Sys.time()
  
  # assign simulation conditions as attributes to s1ests to use for stage2
  attr(s1ests, "MCSampID") <- MCSampID
  attr(s1ests, "n") <- n
  attr(s1ests, "G") <- G
  attr(s1ests, "priorType") <- priorType
  attr(s1ests, "precision") <- ifelse(priorType == "default", "NA", precision)
  attr(s1ests, "iter") <- iter
  attr(s1ests, "RunTime") <- difftime(t1, t0, units = "mins")
  
  attr(s1ests, "mPSRF") <- gelman.diag(As.mcmc.list(s1ests, pars = get("MCMC_pars", envir = s1_env)), 
                          autoburnin = FALSE)$mpsrf
  #FIXME below
  if (savefile) saveRDS(s1ests, paste0("s1_ID", MCSampID, ".nG", G, ".n", n, "_", 
                                    priorType,  "_", ifelse(priorType == "default", "NA", precision), ".rds"))
  return(s1ests)
  
}

lavs1(MCSampID = 1, n = 5, G = 3, priorType = "prophetic", precision = 0.1, iter = 100, savefile = F) -> foo
# lavs1(MCSampID = 1, n = 5, G = 3, priorType = "FIML", precision = 0.1, iter = 100, savefile = F)
# lavs1(MCSampID = 1, n = 5, G = 3, priorType = "default", iter = 100, savefile = F)

#----

# function 8: stage2 SR_SEM in lavaan.srm----

lavs2 <- function(MCSampID, n, G, priorType, precision, s1result) {
  library(lavaan.srm)
  
  s1ests <- s1result[[paste0("ID", MCSampID, ".nG", G, ".n", n, "_", 
                             priorType, "_", ifelse(priorType == "default", "NA", precision), "_s1ests")]]
  s1details <- s1result[[paste0("ID", MCSampID, ".nG", G, ".n", n, "_", 
                                priorType, "_", ifelse(priorType == "default", "NA", precision), "_analDetails")]]
  
  # stage2
  mod <- ' group: 1
  Factor_out =~ 1*V1_out + V2_out + V3_out
  Factor_in =~ 1*V1_in + V2_in + V3_in
  
  Factor_out ~~ Factor_out + Factor_in
  Factor_in ~~ Factor_in
  
  V1_out ~~ V1_out + V1_in
  V1_in ~~ V1_in
  V2_out ~~ V2_out
  V2_in ~~ V2_in
  V3_out ~~ V3_out + V3_in
  V3_in ~~ V3_in
  
  group: 2
  Factor_ij =~ 1*V1_ij + FL2*V2_ij + FL3*V3_ij
  Factor_ji =~ 1*V1_ji + FL2*V2_ji + FL3*V3_ji
  
  Factor_ij ~~ Fcov*Factor_ij + Factor_ji
  Factor_ji ~~ Fcov*Factor_ji
  
  V1_ij ~~ Ivar1*V1_ij
  V1_ji ~~ Ivar1*V1_ji
  V2_ij ~~ Ivar2*V2_ij + V2_ji
  V2_ji ~~ Ivar2*V2_ji
  V3_ij ~~ Ivar3*V3_ij + V3_ji
  V3_ji ~~ Ivar3*V3_ji
  '
  
  fit <- lavaan.srm(model = mod, data = s1ests, component = c("case", "dyad"), posterior.est = "mean")
  
  # fit@test$browne.residual.adf$stat # overall test statistic?
  # fit@test$browne.residual.adf$stat.group # level-specific test statistic?
  
  #TODO save overall fit measure (chi-sq)
  
} #TODO figure out how to cbind() MCSampID, n, G, analType, etc to the result

#----

# function 9: FIML for SR-SEM effects in`srm`----

ogsrm <- function(MCSampID, n, G, rr.vars = c("V1", "V2", "V3"), IDout = "Actor", 
                  IDin = "Partner", IDgroup = "Group", savefile = F) {
  library(srm)
  
  rr.data <- genGroups(MCSampID = MCSampID, n = n, G = G, rr.vars = rr.vars)
  
  # model
  mod_srm <- '
  %Person 
  Factor@A =~ 1*V1@A + V2@A + V3@A
  Factor@P =~ 1*V1@P + V2@P + V3@P
  
  Factor@A ~~ Factor@A + Factor@P
  Factor@P ~~ Factor@P
  
  V1@A ~~ V1@A + V1@P + 0*V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V1@P ~~ V1@P + 0*V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V2@A ~~ V2@A + V2@P + 0*V3@A + 0*V3@P
  V2@P ~~ V2@P + 0*V3@A + 0*V3@P
  V3@A ~~ V3@A + V3@P
  V3@P ~~ V3@P
  
  %Dyad 
  Factor@AP =~ 1*V1@AP + FL2*V2@AP + FL3*V3@AP
  Factor@PA =~ 1*V1@PA + FL2*V2@PA + FL3*V3@PA
  
  Factor@AP ~~ Fvar*Factor@AP + Factor@PA
  Factor@PA ~~ Fvar*Factor@PA
  
  V1@AP ~~ Ind1var*V1@AP + 0*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V1@PA ~~ Ind1var*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V2@AP ~~ Ind2var*V2@AP + V2@PA + 0*V3@AP + 0*V3@PA
  V2@PA ~~ Ind2var*V2@PA + 0*V3@AP + 0*V3@PA
  V3@AP ~~ Ind3var*V3@AP + V3@PA
  V3@PA ~~ Ind3var*V3@PA
  ' #TODO check the model code --- equality constraints, constraints to 0 and 1---are they all correct?
  
  fit_srm <- srm(mod_srm, rr.data, 
                 person_names = c(IDout, IDin), 
                 rrgroup_name = IDgroup, verbose = FALSE)
  
  if (!fit_srm$res_opt$converged) return(NULL)
  
  srm.parm <- fit_srm$parm.table
  
  popVals <- rbind(getSigma(return_mats = F)$popVals_c, getSigma(return_mats = F)$popVals_d)
  
  results_srm <- merge(srm.parm, popVals, by = "par_names")
  results_srm$level <- gsub("U", "case", gsub("D", "dyad", results_srm$level))
  results_srm <- results_srm[, c("par_names", "pop", "level", "est", "se")]
  
  # coverage
  results_srm$ci.lower <- results_srm$est - 1.96*results_srm$se
  results_srm$ci.upper <- results_srm$est + 1.96*results_srm$se
  results_srm$coverage <- results_srm$ci.lower < results_srm$pop & results_srm$pop < results_srm$ci.upper
  
  results_srm$bias <- results_srm$est - results_srm$pop
  results_srm$RB <- results_srm$bias / results_srm$pop
  
  results_srm <- cbind(MCSampID, n, G, priorType = "ogsrm", precision = "NA", 
                       iter = "NA", results_srm)
  
  if (savefile) saveRDS(results_srm, file = paste0("ID", MCSampID, ".nG", G, ".n", 
                                                   n, "-srmML-og-og", ".rds"))
  
  return(results_srm)
}

#----


