## Aditi M. Bhangale
## Last updated: 6 June 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

## FIML1S IN `srm`

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem/sim_code")
# getwd()

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = 1, n = 5, G = 3)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision <- 0.1
# priorType = "prophetic"
# iter = 100

# function 12: FIML for SR-SEM effects in`srm`----

#TODO check if any bugs & fix

ogsrm <- function(MCSampID, n, G, rr.vars = c("V1", "V2", "V3"), IDout = "Actor", 
                  IDin = "Partner", IDgroup = "Group", savefile = F) {
  library(srm)
  
  rr.data <- genGroups(MCSampID = MCSampID, n = n, G = G# , rr.vars = rr.vars
  )
  
  # model
  mod_srm <- '
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
  
  V1@AP ~~ Ind1var*V1@AP + 0*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V1@PA ~~ Ind1var*V1@PA + 0*V2@AP + 0*V2@PA + 0*V3@AP + 0*V3@PA
  V2@AP ~~ Ind2var*V2@AP + V2@PA + 0*V3@AP + 0*V3@PA
  V2@PA ~~ Ind2var*V2@PA + 0*V3@AP + 0*V3@PA
  V3@AP ~~ Ind3var*V3@AP + V3@PA
  V3@PA ~~ Ind3var*V3@PA
  ' #TODO check the model code --- equality constraints, constraints to 0 and 1---are they all correct? --- acc. to Nestler code
  
  fit_srm <- srm(mod_srm, rr.data, 
                 person_names = c(IDout, IDin), 
                 rrgroup_name = IDgroup, verbose = FALSE)
  
  #TODO fit a saturated model and compute the chi-square statistic using the -2*log likelihoods
  
  if (!fit_srm$res_opt$converged) return(NULL)
  
  srm.parm <- fit_srm$parm.table
  
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

#----

