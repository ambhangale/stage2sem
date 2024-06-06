## Aditi M. Bhangale
## Last updated: 6 June 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

## PREPARE STAGE1 `lavaan.srm` RESULTS FOR STAGE 2 + STAGE 2 `lavaan.srm`

# rm(list = ls())

# setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem")
# getwd()

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = 1, n = 5, G = 3)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision <- 0.1
# priorType = "prophetic"
# iter = 100

# function 10: function to flag outliers in lavs1 output----
#TODO
#----

# function 11: stage2 SR_SEM in lavaan.srm----

#TODO figure out how to cbind() MCSampID, n, G, analType, etc to the result

### THE MODEL IS CORRECT, THIS IS JUST HOW NESTLER ET AL., FIT IT

lavs2 <- function(s1ests) {
  library(lavaan.srm)
  
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
  
  Factor_ij ~~ Fvar*Factor_ij + Factor_ji
  Factor_ji ~~ Fvar*Factor_ji
  
  V1_ij ~~ Ivar1*V1_ij
  V1_ji ~~ Ivar1*V1_ji
  V2_ij ~~ Ivar2*V2_ij + V2_ji
  V2_ji ~~ Ivar2*V2_ji
  V3_ij ~~ Ivar3*V3_ij + V3_ji
  V3_ji ~~ Ivar3*V3_ji
  '
  
  fit <- lavaan.srm(model = mod, data = s1ests, component = c("case", "dyad"), posterior.est = "mean")
  # lavaan.srm(model = mod, data = s1ests, component = c("case", "dyad"), posterior.est = "mean", test = c("satorra.bentler","scaled.shifted","browne.residual.adf"))
  
  if(lavInspect(fit, "converged")) {
    s2ests <- parameterEstimates(fit)
    s2ests$lhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests$lhs))))
    s2ests$rhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests$rhs))))
    s2ests$par_names <- paste0(s2ests$lhs, s2ests$op, s2ests$rhs)
    # remove factor loadings constrained for identification and redundant rows
    s2ests <- s2ests[!(s2ests$par_names %in% c("Factor@A=~V1@A", "Factor@P=~V1@P", 
                                               "Factor@AP=~V1@AP", "Factor@PA=~V1@PA",
                                               "Factor@PA=~V2@PA", "Factor@PA=~V3@PA",
                                               "Factor@PA~~Factor@PA", "V1@PA~~V1@PA",
                                               "V2@PA~~V2@PA", "V3@PA~~V3@PA")), ]
    popVals <- getSigma(return_mats = F)
    s2ests <- merge(s2ests, popVals, by = "par_names")
    s2ests$group <- ifelse(s2ests$group == 1, "case", "dyad")
    s2ests$MCSampID <- attr(s1ests, "MCSampID")
    s2ests$n <- attr(s1ests, "n"); s2ests$G <- attr(s1ests, "G")
    s2ests$condition <- paste0(attr(s1ests, "n"), "-", attr(s1ests, "G"))
    s2ests$analType <- "2SMLE"
    s2ests$s1priorType <- paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision"))
    s2ests$s1iter <- attr(s1ests, "iter")
    s2ests$s1mPSRF <- attr(s1ests, "mPSRF")
    s2ests$coverage <- s2ests$ci.lower < s2ests$pop & s2ests$pop < s2ests$ci.upper
    
    s2ests <- s2ests[, c("MCSampID", "n", "G", "condition", "analType", "s1priorType", 
                         "s1iter", "s1mPSRF", "par_names", "group", 
                         "pop", "est", "se", "ci.lower", "ci.upper", "coverage")] # remove non-redundant rows and reorder
    
    s2mod <- c(MCSampID = attr(s1ests, "MCSampID"), n = attr(s1ests, "n"),
               G = attr(s1ests, "G"), 
               condition = paste0(attr(s1ests, "n"), "-", attr(s1ests, "G")), 
               analType = "2SMLE", 
               s1priorType  = paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision")), 
               s1iter = attr(s1ests, "iter"), s1mPSRF = attr(s1ests, "mPSRF"),
               combi.ADF = fit@test$browne.residual.adf$stat, 
               combi.df = fit@test$browne.residual.adf$df, 
               combi.p = fit@test$browne.residual.adf$pvalue,
               case.ADF = fit@test$browne.residual.adf$stat.group[1],
               dyad.ADF = fit@test$browne.residual.adf$stat.group[2]) # TODO also add columns in s2ests?
    
  } else {
    s2result <- NULL
  }
  
  #TODO saveRDS and if result is NULL, then what?
  
  #TODO save overall fit measure (chi-sq)
  
}

#----

