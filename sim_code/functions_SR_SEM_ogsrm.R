## Aditi M. Bhangale
## Last updated: 24 June 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

## FIML1S IN `srm`

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem/sim_code")
# getwd()

# MCSampID = 1; n = 6; G = 10
# rr.data <- genGroups(MCSampID = 1, n = 6, G = 10)
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"

source("functions_SR_SEM_lavs1.R") # for data generation functions

# function 1: FIML for SR-SEM effects in`srm`----

ogsrm <- function(MCSampID, n, G, IDout = "Actor", 
                  IDin = "Partner", IDgroup = "Group", savefile = F) {
  library(srm)
  t0 <- Sys.time()
  rr.data <- genGroups(MCSampID = MCSampID, n = n, G = G)
  
  # model
  mod_combi <- '
  %Person 
  Factor@A =~ 1*V1@A + V2@A + V3@A
  Factor@P =~ 1*V1@P + V2@P + V3@P
  
  Factor@A ~~ Factor@A + Factor@P
  Factor@P ~~ Factor@P
  
  V1@A ~~ V1@A + V1@P + 0*V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V1@P ~~ V1@P + 0*V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V2@A ~~ V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V2@P ~~ V2@P + 0*V3@A + 0*V3@P
  V3@A ~~ V3@A + V3@P
  V3@P ~~ V3@P
  
  %Dyad 
  Factor@AP =~ 1*V1@AP + FL2*V2@AP + FL3*V3@AP
  Factor@PA =~ 1*V1@PA + FL2*V2@PA + FL3*V3@PA
  
  Factor@AP ~~ Fvar*Factor@AP + Factor@PA
  Factor@PA ~~ Fvar*Factor@PA
  
  V1@AP ~~ Res1var*V1@AP + 0*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V1@PA ~~ Res1var*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V2@AP ~~ Res2var*V2@AP + V2@PA + 0*V3@AP + 0*V3@PA
  V2@PA ~~ Res2var*V2@PA + 0*V3@AP + 0*V3@PA
  V3@AP ~~ Res3var*V3@AP + V3@PA
  V3@PA ~~ Res3var*V3@PA
  '
  
  fit_combi <- srm(mod_combi, data = rr.data, 
                 person_names = c(IDout, IDin), 
                 rrgroup_name = IDgroup, verbose = FALSE)
  
  
  mod_sat <- '
  %Person
  f1@A =~ 1*V1@A
  f2@A =~ 1*V2@A
  f3@A =~ 1*V3@A
  
  f1@P =~ 1*V1@P
  f2@P =~ 1*V2@P
  f3@P =~ 1*V3@P
  
  V1@A ~~ 0*V1@A + 0*V1@P
  V1@P ~~ 0*V1@P
  V2@A ~~ 0*V2@A + 0*V2@P
  V2@P ~~ 0*V2@P
  V3@A ~~ 0*V3@A + 0*V3@P
  V3@P ~~ 0*V3@P
  
  f1@A ~~ f1@A + f2@A + f3@A + f1@P + f2@P + f3@P
  f2@A ~~        f2@A + f3@A + f1@P + f2@P + f3@P
  f3@A ~~               f3@A + f1@P + f2@P + f3@P
  f1@P ~~                      f1@P + f2@P + f3@P
  f2@P ~~                             f2@P + f3@P
  f3@P ~~                                    f3@P
  
  
  %Dyad
  f1@AP =~ 1*V1@AP
  f2@AP =~ 1*V2@AP
  f3@AP =~ 1*V3@AP
  
  # f1@PA =~ 1*V1@PA
  # f2@PA =~ 1*V2@PA
  # f3@PA =~ 1*V3@PA
  
  V1@AP ~~ 0*V1@AP + 0*V1@PA
  V1@PA ~~ 0*V1@PA
  V2@AP ~~ 0*V2@AP + 0*V2@PA
  V2@PA ~~ 0*V2@PA
  V3@AP ~~ 0*V3@AP + 0*V3@PA
  V3@PA ~~ 0*V3@PA
  
  f1@AP ~~ relvar1*f1@AP + intra12*f2@AP + intra13*f3@AP + dyad11*f1@PA  + inter12*f2@PA + inter13*f3@PA
  f2@AP ~~                 relvar2*f2@AP + intra23*f3@AP + inter12*f1@PA +  dyad22*f2@PA + inter23*f3@PA
  f3@AP ~~                                 relvar3*f3@AP + inter13*f1@PA + inter23*f2@PA +  dyad33*f3@PA
  f1@PA ~~                                                 relvar1*f1@PA + intra12*f2@PA + intra13*f3@PA
  f2@PA ~~                                                                 relvar2*f2@PA + intra23*f3@PA
  f3@PA ~~                                                                                 relvar3*f3@PA
  '
  
  fit_sat <- srm(mod_sat, data = rr.data, 
                 person_names = c(IDout, IDin), 
                 rrgroup_name = IDgroup, verbose = FALSE)
  
  mod_case <- ' 
  %Person 
  Factor@A =~ 1*V1@A + V2@A + V3@A
  Factor@P =~ 1*V1@P + V2@P + V3@P
  
  Factor@A ~~ Factor@A + Factor@P
  Factor@P ~~ Factor@P
  
  V1@A ~~ V1@A + V1@P + 0*V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V1@P ~~ V1@P + 0*V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V2@A ~~ V2@A + 0*V2@P + 0*V3@A + 0*V3@P
  V2@P ~~ V2@P + 0*V3@A + 0*V3@P
  V3@A ~~ V3@A + V3@P
  V3@P ~~ V3@P
  
  %Dyad
  f1@AP =~ 1*V1@AP
  f2@AP =~ 1*V2@AP
  f3@AP =~ 1*V3@AP
  
  # f1@PA =~ 1*V1@PA
  # f2@PA =~ 1*V2@PA
  # f3@PA =~ 1*V3@PA
  
  V1@AP ~~ 0*V1@AP + 0*V1@PA
  V1@PA ~~ 0*V1@PA
  V2@AP ~~ 0*V2@AP + 0*V2@PA
  V2@PA ~~ 0*V2@PA
  V3@AP ~~ 0*V3@AP + 0*V3@PA
  V3@PA ~~ 0*V3@PA
  
  f1@AP ~~ relvar1*f1@AP + intra12*f2@AP + intra13*f3@AP + dyad11*f1@PA  + inter12*f2@PA + inter13*f3@PA
  f2@AP ~~                 relvar2*f2@AP + intra23*f3@AP + inter12*f1@PA +  dyad22*f2@PA + inter23*f3@PA
  f3@AP ~~                                 relvar3*f3@AP + inter13*f1@PA + inter23*f2@PA +  dyad33*f3@PA
  f1@PA ~~                                                 relvar1*f1@PA + intra12*f2@PA + intra13*f3@PA
  f2@PA ~~                                                                 relvar2*f2@PA + intra23*f3@PA
  f3@PA ~~                                                                                 relvar3*f3@PA
  '
  
  fit_case <- srm(mod_case, data = rr.data, 
                 person_names = c(IDout, IDin), 
                 rrgroup_name = IDgroup, verbose = FALSE)
  
  mod_dyad <- '
  %Person
  f1@A =~ 1*V1@A
  f2@A =~ 1*V2@A
  f3@A =~ 1*V3@A
  
  f1@P =~ 1*V1@P
  f2@P =~ 1*V2@P
  f3@P =~ 1*V3@P
  
  V1@A ~~ 0*V1@A + 0*V1@P
  V1@P ~~ 0*V1@P
  V2@A ~~ 0*V2@A + 0*V2@P
  V2@P ~~ 0*V2@P
  V3@A ~~ 0*V3@A + 0*V3@P
  V3@P ~~ 0*V3@P
  
  f1@A ~~ f1@A + f2@A + f3@A + f1@P + f2@P + f3@P
  f2@A ~~        f2@A + f3@A + f1@P + f2@P + f3@P
  f3@A ~~               f3@A + f1@P + f2@P + f3@P
  f1@P ~~                      f1@P + f2@P + f3@P
  f2@P ~~                             f2@P + f3@P
  f3@P ~~                                    f3@P
  
  %Dyad 
  Factor@AP =~ 1*V1@AP + FL2*V2@AP + FL3*V3@AP
  Factor@PA =~ 1*V1@PA + FL2*V2@PA + FL3*V3@PA
  
  Factor@AP ~~ Fvar*Factor@AP + Factor@PA
  Factor@PA ~~ Fvar*Factor@PA
  
  V1@AP ~~ Res1var*V1@AP + 0*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V1@PA ~~ Res1var*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V2@AP ~~ Res2var*V2@AP + V2@PA + 0*V3@AP + 0*V3@PA
  V2@PA ~~ Res2var*V2@PA + 0*V3@AP + 0*V3@PA
  V3@AP ~~ Res3var*V3@AP + V3@PA
  V3@PA ~~ Res3var*V3@PA
  '
  
  fit_dyad <- srm(mod_dyad, data = rr.data, 
                  person_names = c(IDout, IDin), 
                  rrgroup_name = IDgroup, verbose = FALSE)
  
  # PARAMETER ESTIMATES
  if (fit_combi$res_opt$converged) {
  srm.parm <- fit_combi$parm.table
  
  popVals <- getSigma(return_mats = F)
  
  srm_ests <- merge(srm.parm, popVals, by = "par_names")
  srm_ests$level <- gsub("U", "case", gsub("D", "dyad", srm_ests$level))
  srm_ests <- srm_ests[, c("par_names", "level", "pop", "est", "se")]
  
  # coverage
  srm_ests$ci.lower <- srm_ests$est - qnorm(.975)*srm_ests$se 
  srm_ests$ci.upper <- srm_ests$est + qnorm(.975)*srm_ests$se
  srm_ests$coverage <- srm_ests$ci.lower < srm_ests$pop & srm_ests$pop < srm_ests$ci.upper
  
  srm_ests$bias <- srm_ests$est - srm_ests$pop
  srm_ests$RB <- srm_ests$bias / srm_ests$pop
  
  srm_ests <- cbind(MCSampID = MCSampID, n = n, G = G, condition = paste0(n, "-", G), 
                    analType = "FIML1S", s1priorType = "NA", 
                       s1iter = "NA", s1mPSRF = "NA", srm_ests) 
  } else {
    if (fit_case$res_opt$converged) {
      srm.parm_case <- fit_case$parm.table
      popVals <- getSigma(return_mats = F)
      
      srm_ests_case <- merge(srm.parm_case, popVals, by = "par_names")
      srm_ests_case$level <- "case"
      srm_ests_case <- srm_ests_case[, c("par_names", "level", "pop", "est", "se")]
      
      # coverage
      srm_ests_case$ci.lower <- srm_ests_case$est - qnorm(.975)*srm_ests_case$se 
      srm_ests_case$ci.upper <- srm_ests_case$est + qnorm(.975)*srm_ests_case$se
      srm_ests_case$coverage <- srm_ests_case$ci.lower < srm_ests_case$pop & srm_ests_case$pop < srm_ests_case$ci.upper
      
      srm_ests_case$bias <- srm_ests_case$est - srm_ests_case$pop
      srm_ests_case$RB <- srm_ests_case$bias / srm_ests_case$pop
      
      srm_ests_case <- cbind(MCSampID = MCSampID, n = n, G = G, condition = paste0(n, "-", G), 
                        analType = "FIML1S", s1priorType = "NA", 
                        s1iter = "NA", s1mPSRF = "NA", srm_ests_case) 
    } else {
      srm_ests_case <- NULL
    }
    
    if (fit_dyad$res_opt$converged) {
      srm.parm_dyad <- fit_dyad$parm.table
      popVals <- getSigma(return_mats = F)
      
      srm_ests_dyad <- merge(srm.parm_dyad, popVals, by = "par_names")
      srm_ests_dyad$level <- "dyad"
      srm_ests_dyad <- srm_ests_dyad[, c("par_names", "level", "pop", "est", "se")]
      
      # coverage
      srm_ests_dyad$ci.lower <- srm_ests_dyad$est - qnorm(.975)*srm_ests_dyad$se 
      srm_ests_dyad$ci.upper <- srm_ests_dyad$est + qnorm(.975)*srm_ests_dyad$se
      srm_ests_dyad$coverage <- srm_ests_dyad$ci.lower < srm_ests_dyad$pop & srm_ests_dyad$pop < srm_ests_dyad$ci.upper
      
      srm_ests_dyad$bias <- srm_ests_dyad$est - srm_ests_dyad$pop
      srm_ests_dyad$RB <- srm_ests_dyad$bias / srm_ests_dyad$pop
      
      srm_ests_dyad <- cbind(MCSampID = MCSampID, n = n, G = G, condition = paste0(n, "-", G), 
                             analType = "FIML1S", s1priorType = "NA", 
                             s1iter = "NA", s1mPSRF = "NA", s1Reff_outlier = "NA",
                             s1mPSRF_outlier = "NA", srm_ests_dyad) 
    } else {
      srm_ests_dyad <- NULL
    }
    srm_ests <- rbind(srm_ests_case, srm_ests_dyad)
  }
  
  # MODEL FIT
  dev_sat <- fit_sat$dev
  dev_combi <- fit_combi$dev
  dev_case <- fit_case$dev
  dev_dyad <- fit_dyad$dev
  
  df_combi <- 9 
  df_case <- 6
  df_dyad <- 3
  
  if (fit_sat$res_opt$converged) {
    if (fit_combi$res_opt$converged) {
      stat_combi <- dev_combi - dev_sat
      p_combi <- pchisq(stat_combi, df = df_combi, lower.tail = F)
    } else {
      stat_combi <- p_combi <- NULL
    }
    if (fit_case$res_opt$converged) {
      stat_case <- dev_case - dev_sat
      p_case <- pchisq(stat_case, df = df_case, lower.tail = F)
    } else {
      stat_case <- p_case <- NULL
    }
    if (fit_dyad$res_opt$converged) {
      stat_dyad <- dev_dyad - dev_sat
      p_dyad <- pchisq(stat_dyad, df = df_dyad, lower.tail = F)
    } else {
      stat_dyad <- p_dyad <- NULL
    }
  } else {
    stat_combi <- stat_case <- stat_dyad <- p_combi <- p_case <- p_dyad <- NULL
  }
  
  srm_LRT <- data.frame(cbind(MCSampID = MCSampID, n = n, G = G, condition = paste0(n, "-", G), 
                              analType = "FIML1S", s1priorType = "NA", s1iter = "NA", s1mPSRF = "NA", 
                              fitStat.type = "srm.LRT",
                              combi.stat = ifelse(is.null(stat_combi), "NA", stat_combi),
                              combi.df = df_combi,
                              combi.p = ifelse(is.null(p_combi), "NA", p_combi),
                              case.stat = ifelse(is.null(stat_case), "NA", stat_case),
                              case.df = df_case, 
                              case.p = ifelse(is.null(p_case), "NA", p_case),
                              dyad.stat = ifelse(is.null(stat_dyad), "NA", stat_dyad),
                              dyad.df = df_dyad, 
                              dyad.p = ifelse(is.null(p_dyad), "NA", p_dyad)))
  
  t1 <- Sys.time()
  
  srm_ests$RunTime <- difftime(t1, t0, units = "mins")
  
  srm_result <- list(srm_ests = srm_ests, srm_mod = srm_LRT)
  
  if (savefile) saveRDS(srm_result, file = paste0("ID", MCSampID, ".nG", G, ".n", 
                                                   n, "-srmML-og-og", ".rds"))
  
  return(srm_result)
}

# ogsrm(MCSampID = 1, n = 6, G = 10, savefile = T)

#----

# function 2: create runsim files for ogsrm----

makeRunsim <- function(nSamps, n, G, sim) {
  runsimfile <- paste0('## Aditi M. Bhangale
## Last updated:', Sys.Date(), 
                       
'\n# Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data

# runsim_ogsrm_n', n, '_G', G, '_', sim,'

source("functions_SR_SEM_ogsrm.R")

# specify conditions\n
                       
FIML1S_grid <- expand.grid(MCSampID = 1:', nSamps, ', n = ', n, ', G = ', G, ')\n
FIML1S_grid$row_num <- 1:nrow(FIML1S_grid)

# prepare parallel processing\n
library(doSNOW)

nClus <- 124
cl <- makeCluster(nClus)
registerDoSNOW(cl)

# run simulation\n',
paste0('ogResult <- foreach(row_num = 1:nrow(FIML1S_grid),
                    .packages = c("mnormt", "parallel", "portableParallelSeeds", 
                                  "srm")) %dopar% {
                                    
                                    out <- try(ogsrm(MCSampID = FIML1S_grid[row_num, ]$MCSampID, 
                                    n = FIML1S_grid[row_num, ]$n, 
                                                    G = FIML1S_grid[row_num, ]$G), silent = T)
                                    if(inherits(out, "try-error")) out <- NULL
                                    
                                    return(out)
                                  }
         
   # close cluster\n
   stopCluster(cl)
   
   saveRDS(ogResult, paste0("results_ogsrm_n', n, '_G', G, '_',sim, '_",Sys.Date(),".rds"))
         ')
  )
  
  cat(runsimfile, file = paste0("runsim_ogsrm_n", n, "_G", G, "_", sim, ".R"))
  invisible(NULL)
}

# makeRunsim(nSamps = 500, n = 6, G = 10, sim = "sim1")
# makeRunsim(nSamps = 500, n = 6, G = 25, sim = "sim1")
# makeRunsim(nSamps = 500, n = 8, G = 10, sim = "sim1")
# makeRunsim(nSamps = 500, n = 8, G = 25, sim = "sim1")
# makeRunsim(nSamps = 500, n = 10, G = 10, sim = "sim1")
# makeRunsim(nSamps = 500, n = 10, G = 25, sim = "sim1")
# makeRunsim(nSamps = 500, n = 20, G = 10, sim = "sim1")
# makeRunsim(nSamps = 500, n = 20, G = 25, sim = "sim1")

#----

# function 3: create shell files for ogsrm----

makeShSnellius <- function(n, G, sim, wallTime) {
  shell <- paste0('#!/bin/bash

#SBATCH -J ', paste0("ogsrm_n", n, "_G", G, "_", sim),'
#SBATCH -e .', paste0("ogsrm_n", n, "_G", G, "_", sim),'.SERR
#SBATCH -o .', paste0("ogsrm_n", n, "_G", G, "_", sim),'.SOUT
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

cp $HOME/SR-SEM/stage2sem/functions_SR_SEM_lavs1.R "$TMPDIR" 
cp $HOME/SR-SEM/stage2sem/functions_SR_SEM_ogsrm.R "$TMPDIR" 
cp $HOME/SR-SEM/stage2sem/', paste0("runsim_ogsrm_n", n, "_G", G, "_", sim, ".R"),' "$TMPDIR"

Rscript --vanilla ', paste0("runsim_ogsrm_n", n, "_G", G, "_", sim, ".R"),'

cp "$TMPDIR"/*.rds $HOME/SR-SEM/stage2sem/' 

)
  cat(shell, file = paste0("shell_ogsrm_n", n, "_G", G, "_", sim, ".sh"))
  invisible(NULL)
}

# makeShSnellius(n = 6, G = 10, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 6, G = 25, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 8, G = 10, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 8, G = 25, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 10, G = 10, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 10, G = 25, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 20, G = 10, sim = "sim1", wallTime = "5-00:00:00")
# makeShSnellius(n = 20, G = 25, sim = "sim1", wallTime = "5-00:00:00")

#----

