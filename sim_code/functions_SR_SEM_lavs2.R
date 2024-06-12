## Aditi M. Bhangale
## Last updated: 12 June 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

## PREPARE STAGE1 `lavaan.srm` RESULTS FOR STAGE 2 + STAGE 2 `lavaan.srm`

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

# function 10: function to flag outliers in lavs1 output----
#TODO
#----

# function 11: stage2 SR_SEM in lavaan.srm----
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
  
  fit <- lavaan.srm(model = mod, data = s1ests, component = c("case", "dyad"), posterior.est = "mean", 
  test = c("satorra.bentler","scaled.shifted","browne.residual.adf"))
  
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
    names(s2ests)[names(s2ests) == "group"] <- "level"
    s2ests$MCSampID <- attr(s1ests, "MCSampID")
    s2ests$n <- attr(s1ests, "n"); s2ests$G <- attr(s1ests, "G")
    s2ests$condition <- paste0(attr(s1ests, "n"), "-", attr(s1ests, "G"))
    s2ests$analType <- "2SMLE"
    s2ests$s1priorType <- paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision"))
    s2ests$s1iter <- attr(s1ests, "iter")
    s2ests$s1mPSRF <- attr(s1ests, "mPSRF")
    s2ests$coverage <- s2ests$ci.lower < s2ests$pop & s2ests$pop < s2ests$ci.upper
    
    s2ests <- s2ests[, c("MCSampID", "n", "G", "condition", "analType", "s1priorType", 
                         "s1iter", "s1mPSRF", "par_names", "level", 
                         "pop", "est", "se", "ci.lower", "ci.upper", "coverage")] # remove non-redundant rows and reorder
    
    fitStats <- fit@test # save all test statistics in a single object
    
    standard <- c(fitStat.type = fitStats$standard$test, 
                  combi.stat = fitStats$standard$stat, 
                  combi.df = fitStats$standard$df, 
                  combi.p = fitStats$standard$pvalue,
                  case.stat = fitStats$standard$stat.group[1], 
                  case.df = 6, case.p = pchisq(fitStats$standard$stat.group[1], 
                                               df = 6, lower = F),
                  dyad.stat = fitStats$standard$stat.group[2], 
                  dyad.df = 12, dyad.p = pchisq(fitStats$standard$stat.group[2], 
                                                df = 12, lower = F)) # Standard chi-square test
    
    adf <- c(fitStat.type = fitStats$browne.residual.adf$test, 
             combi.stat = fitStats$browne.residual.adf$stat, 
             combi.df = fitStats$browne.residual.adf$df, 
             combi.p = fitStats$browne.residual.adf$pvalue,
             case.stat = fitStats$browne.residual.adf$stat.group[1], 
             case.df = 6, case.p = pchisq(fitStats$browne.residual.adf$stat.group[1], 
                                          df = 6, lower = F),
             dyad.stat = fitStats$browne.residual.adf$stat.group[2], 
             dyad.df = 12, dyad.p = pchisq(fitStats$browne.residual.adf$stat.group[2], 
                                           df = 12, lower = F)) # Browne's residual-based ADF
    
    sb <- c(fitStat.type = fitStats$satorra.bentler$test, 
            combi.stat = fitStats$satorra.bentler$stat, 
            combi.df = fitStats$satorra.bentler$df, 
            combi.p = fitStats$satorra.bentler$pvalue,
            case.stat = fitStats$satorra.bentler$stat.group[1], 
            case.df = 6, case.p = pchisq(fitStats$satorra.bentler$stat.group[1], 
                                         df = 6, lower = F),
            dyad.stat = fitStats$satorra.bentler$stat.group[2], 
            dyad.df = 12, dyad.p = pchisq(fitStats$satorra.bentler$stat.group[2], 
                                          df = 12, lower = F)) # Satorra and Bentler's corrected statistic
    ss <- c(fitStat.type = fitStats$scaled.shifted$test, 
            combi.stat = fitStats$scaled.shifted$stat, 
            combi.df = fitStats$scaled.shifted$df, 
            combi.p = fitStats$scaled.shifted$pvalue,
            case.stat = fitStats$scaled.shifted$stat.group[1], 
            case.df = 6, case.p = pchisq(fitStats$scaled.shifted$stat.group[1], 
                                         df = 6, lower = F),
            dyad.stat = fitStats$scaled.shifted$stat.group[2], 
            dyad.df = 12, dyad.p = pchisq(fitStats$scaled.shifted$stat.group[2], 
                                          df = 12, lower = F)) # Scaled-shifted statistic
    
    Ncase <- attr(s1ests, "nobs")["case"]; Ndyad <- attr(s1ests, "nobs")["dyad"]
    Ntotal <- Ncase + Ndyad
    
    yb_corrected.combistat <- fitStats$browne.residual.adf$stat / (1 + (fitStats$browne.residual.adf$stat/Ntotal))
    yb_corrected.casestat <- fitStats$browne.residual.adf$stat.group[1] / (1 + (fitStats$browne.residual.adf$stat.group[1]/(Ncase)))
    yb_corrected.dyadstat <- fitStats$browne.residual.adf$stat.group[2] / (1 + (fitStats$browne.residual.adf$stat.group[2]/(Ndyad)))
    yb_corrected.combip <- pchisq(yb_corrected.combistat, df = 18, lower = F)
    yb_corrected.casep <- pchisq(yb_corrected.casestat, df = 6, lower = F)
    yb_corrected.dyadp <- pchisq(yb_corrected.dyadstat, df = 12, lower = F)
    yb_corrected <- c(fitStat.type = "yuan.bentler.corrected.adf", 
                  combi.stat = yb_corrected.combistat, 
                  combi.df = 18, 
                  combi.p = yb_corrected.combip,
                  case.stat = yb_corrected.casestat, 
                  case.df = 6, case.p = yb_corrected.casep,
                  dyad.stat = yb_corrected.dyadstat, 
                  dyad.df = 12, dyad.p = yb_corrected.dyadp) # Yuan and Bentler's small-sample correction for ADF
    
    
    yb_F.combistat <- ((Ntotal - 18)/(Ntotal*18))*fitStats$browne.residual.adf$stat
    yb_F.casestat <- ((Ncase - 6)/(Ntotal*6))*fitStats$browne.residual.adf$stat.group[1]
    yb_F.dyadstat <- ((Ncase - 12)/(Ntotal*12))*fitStats$browne.residual.adf$stat.group[2]
    yb_F.combip <- pchisq(yb_F.combistat, df = 18, lower = F)
    yb_F.casep <- pchisq(yb_F.casestat, df = 6, lower = F)
    yb_F.dyadp <- pchisq(yb_F.dyadstat, df = 12, lower = F)
    yb_F <- c(fitStat.type = "yuan.bentler.F", 
              combi.stat = yb_F.combistat, 
              combi.df = paste0(18, ",", (Ntotal-18)), 
              combi.p = yb_F.combip,
              case.stat = yb_F.casestat, 
              case.df = paste0(6, ",", (Ncase-6)), case.p = yb_F.casep,
              dyad.stat = yb_F.dyadstat, 
              dyad.df = paste0(12, ",", (Ndyad-12)), dyad.p = yb_F.dyadp) # Yuan and Bentler's F statistic based on ADF
    
    s2mod <- data.frame(cbind(MCSampID = attr(s1ests, "MCSampID"),
                              n = attr(s1ests, "n"),
                              G = attr(s1ests, "G"),
                              condition = paste0(attr(s1ests, "n"), "-", attr(s1ests, "G")),
                              analType = "2SMLE",
                              s1priorType = paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision")),
                              s1iter = attr(s1ests, "iter"),
                              s1mPSRF = attr(s1ests, "mPSRF"), rbind(standard, adf, sb, ss, yb_corrected, yb_F)))
    
    
    #TODO don't forget to add redundant columns to make sure you can rbind() with the ogsrm() output later on
    
  } else {
    s2result <- NULL
  }
  
  #TODO saveRDS and if result is NULL, then what?
  
  #TODO save overall fit measure (chi-sq)
  
}

#----

