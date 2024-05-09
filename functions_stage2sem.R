## Aditi M. Bhangale
## Last updated: 9 May 2024

# " Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem")
# getwd()

#################################
# FUNCTIONS FOR SIMULATIONS ----
#################################

# function 0: generate level-specific (co)variance matrices----

getSigma <- function(return_mats = TRUE){
  
  ## case effects
  Inames_c <- c("V1@A", "V2@A", "V3@A", "V1@P", "V2@P", "V3@P") # case-level indicator names
  Fnames_c <- c("f1@A", "f1@P") # case-level factor names
  
  # case-level lambda
  LAM_c <- matrix(0, 6, 2)
  LAM_c[1,1] <- 1
  LAM_c[2,1] <- 1.2
  LAM_c[3,1] <- 0.7
  LAM_c[4,2] <- 1
  LAM_c[5,2] <- 0.6
  LAM_c[6,2] <- 0.6
  
  dimnames(LAM_c) <- list(Inames_c, Fnames_c)
  
  # case-level factor cov matrix
  PHI_c <- matrix(c(0.40, 0.05, 0.05, 0.20), 2, 2)
  
  dimnames(PHI_c) <- list(Fnames_c, Fnames_c)
  
  # case-level indicator residual cov matrix
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
  
  ## dyad effects
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
  
  # case-level sigma
  SIGMA_c <- LAM_c %*% PHI_c %*% t(LAM_c) + PSI_c
  # dyad-level sigma
  SIGMA_d <- LAM_d %*% PHI_d %*% t(LAM_d) + PSI_d
  
  ## pop values for satmod
  FsatNames_c <- c("f1@A", "f2@A", "f3@A", "f1@P", "f2@P", "f3@P") # case-level satmod factor names
  FsatNames_d <- c("f1@AP", "f1@PA", "f2@AP", "f2@PA", "f3@AP", "f3@PA") # dyad-level satmod factor names
  
  dimnames(SIGMA_c) <- list(FsatNames_c, FsatNames_c)
  dimnames(SIGMA_d) <- list(FsatNames_d, FsatNames_d)
  
  # case-level corr matrix
  R_c <- cov2cor(SIGMA_c)
  # dyad-level corr matrix
  R_d <- cov2cor(SIGMA_d)
  
  mat_list <- list(LAM_c=LAM_c, PHI_c=PHI_c, PSI_c=PSI_c,
                   LAM_d=LAM_d, PHI_d=PHI_d, PSI_d=PSI_d,
                   SIGMA_c = SIGMA_c, SIGMA_d = SIGMA_d, R_c = R_c, R_d = R_d)
  
  if (return_mats)  return(mat_list)
  
  # creating a dataframe of case-level population values
  pop_c.cov <- as.data.frame(as.table(SIGMA_c))
  colnames(pop_c.cov) <- c("row", "col", "pop.cov")
  pop_c.cov$par_names <- paste0(pop_c.cov$row, "~~", pop_c.cov$col)
  pop_c.cov <- pop_c.cov[, c("par_names", "pop.cov")]
  
  # case-level SDs
  mat_c.SD <- diag(sqrt(diag(SIGMA_c)))
  dimnames(mat_c.SD) <- list(FsatNames_c, FsatNames_c)
  pop_c.SD <- as.data.frame(as.table(mat_c.SD))
  pop_c.SD <- pop_c.SD[pop_c.SD$Freq != 0, ]
  colnames(pop_c.SD) <- c("row", "col", "pop.SD")
  pop_c.SD$par_names <- paste0(pop_c.SD$row, "~~", pop_c.SD$col)
  pop_c.SD <- pop_c.SD[, c("par_names", "pop.SD")]
  rownames(pop_c.SD) <- NULL
  
  # case-level correlations
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

getSigma()
getSigma(return_mats = F)

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

# function 4: thoughtful priors----

#----

# function 5: prophetic priors----

#----

# function 6: ANOVA-based methods-of-moments priors----

#----

# function 7: FIML-based priors----

#----

# function 8: set customised priors for MCMC stage----

#----

# function 9: stage1 in `lavaan.srm`----

#----

# function 10: stage2 SR_SEM in lavaan.srm----

#----

# function 11: FIML for SR-SEM effects in`srm`----

#----


