## Aditi M. Bhangale
## Last updated: 15 June 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

## FIML1S IN `srm`

# rm(list = ls())

setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem/sim_code")
# getwd()

# MCSampID = 1; n = 6; G = 10
# rr.data <- genGroups(MCSampID = 1, n = 6, G = 10)
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"

source("functions_SR_SEM_lavs1.R") # for data generation functions

# function 12: FIML for SR-SEM effects in`srm`----

#TODO check if any bugs & fix

ogsrm <- function(MCSampID, n, G, IDout = "Actor", 
                  IDin = "Partner", IDgroup = "Group", savefile = F) {
  library(srm)
  
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
  
  #### logLik(fit_combi)
  #### logLik(fit_sat)
  ### why the F is the df the same as the number of estimated parameters wtf
  
  #TODO fit a saturated model and compute the chi-square statistic using the -2*log likelihoods
  ## what if one model converges but the other doesn't
  
  if (!fit_combi$res_opt$converged) return(NULL)
  
  srm.parm <- fit_combi$parm.table
  
  popVals <- rbind(getSigma(return_mats = F)$popVals_c, getSigma(return_mats = F)$popVals_d)
  
  results_srm <- merge(srm.parm, popVals, by = "par_names")
  results_srm$level <- gsub("U", "case", gsub("D", "dyad", results_srm$level))
  results_srm <- results_srm[, c("par_names", "pop", "level", "est", "se")]
  
  # coverage
  results_srm$ci.lower <- results_srm$est - 1.96*results_srm$se # FIXME qnorm instead
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

#TODO don't forget to add the redundant columns so you can easily rbind() things

#----

