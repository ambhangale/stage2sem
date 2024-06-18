## Aditi M. Bhangale
## Last updated: 18 June 2024

# "Comparing maximum likelihood to two-stage estimation for structural equation 
# models of social-network data"
## Using Nestler et al., 2020 population values

## PREPARE STAGE1 `lavaan.srm` RESULTS FOR STAGE 2 + STAGE 2 `lavaan.srm`

# rm(list = ls())

setwd("/Users/Aditi_2/Desktop/UvA/SR-SEM_job/stage2sem/sim_code")
# getwd()

# MCSampID = 1; n = 5; G = 3
# rr.data <- genGroups(MCSampID = 1, n = 5, G = 3)
# rr.vars <- c("V1", "V2", "V3")
# IDout <- "Actor"; IDin <- "Partner"; IDgroup <- "Group"
# precision <- 0.1
# priorType = "prophetic"
# iter = 100

# function 10: function to flag outliers in lavs1 output----
#TODO create a logical (TRUE/FALSE) attribute so that lavs2 runs only if it is outliers == FALSE
## otherwise lavs2 returns NULL
#----

# function 11: stage2 SR_SEM in lavaan.srm----

lavs2 <- function(s1ests, savefile = FALSE) {
  library(lavaan.srm)
  
  # stage2
  mod_combi <- ' group: 1
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
  
  fit_combi <- lavaan.srm(model = mod_combi, data = s1ests, component = c("case", "dyad"), posterior.est = "mean", 
  test = c("satorra.bentler","scaled.shifted","browne.residual.adf"))
  
  mod_case <- '
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
  '
  
  fit_case <- lavaan.srm(model = mod_case, data = s1ests, component = "case", posterior.est = "mean", 
                         test = c("satorra.bentler","scaled.shifted","browne.residual.adf"))
  
  mod_dyad <- '
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
  
  fit_dyad <- lavaan.srm(model = mod_dyad, data = s1ests, component = "dyad", posterior.est = "mean", 
                         test = c("satorra.bentler","scaled.shifted","browne.residual.adf"))
  
  # PARAMETER ESTIMATES
  if (lavInspect(fit_combi, "converged")) {
    s2ests <- parameterEstimates(fit_combi)
    s2ests$lhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests$lhs))))
    s2ests$rhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests$rhs))))
    s2ests$par_names <- paste0(s2ests$lhs, s2ests$op, s2ests$rhs)
    ## remove factor loadings constrained for identification and redundant rows
    s2ests <- s2ests[!(s2ests$par_names %in% c("Factor@A=~V1@A", "Factor@P=~V1@P", 
                                               "Factor@AP=~V1@AP", "Factor@PA=~V1@PA",
                                               "Factor@PA=~V2@PA", "Factor@PA=~V3@PA",
                                               "Factor@PA~~Factor@PA", "V1@PA~~V1@PA",
                                               "V2@PA~~V2@PA", "V3@PA~~V3@PA")), ]
    popVals <- getSigma(return_mats = F)
    s2ests <- merge(s2ests, popVals, by = "par_names")
    s2ests$group <- ifelse(s2ests$group == 1, "case", "dyad")
    names(s2ests)[names(s2ests) == "group"] <- "level" # rename 'group' column to 'level'
    s2ests$MCSampID <- attr(s1ests, "MCSampID")
    s2ests$n <- attr(s1ests, "n"); s2ests$G <- attr(s1ests, "G")
    s2ests$condition <- paste0(attr(s1ests, "n"), "-", attr(s1ests, "G"))
    s2ests$analType <- "2SMLE"
    s2ests$s1priorType <- paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision"))
    s2ests$s1iter <- attr(s1ests, "iter")
    s2ests$s1mPSRF <- attr(s1ests, "mPSRF")
    s2ests$coverage <- s2ests$ci.lower < s2ests$pop & s2ests$pop < s2ests$ci.upper
    s2ests$bias <- s2ests$est - s2ests$pop
    s2ests$RB <- s2ests$bias / s2ests$pop
    
    s2ests <- s2ests[, c("MCSampID", "n", "G", "condition", "analType", "s1priorType", 
                         "s1iter", "s1mPSRF", "par_names", "level", 
                         "pop", "est", "se", "ci.lower", "ci.upper", "coverage", "bias", "RB")] # remove non-redundant rows and reorder
  
  } else {
    if (lavInspect(fit_case, "converged")) {
    s2ests_case <- parameterEstimates(fit_case)
    s2ests_case$lhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests_case$lhs))))
    s2ests_case$rhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests_case$rhs))))
    s2ests_case$par_names <- paste0(s2ests_case$lhs, s2ests_case$op, s2ests_case$rhs)
    ## remove factor loadings constrained for identification and redundant rows
    s2ests_case <- s2ests_case[!(s2ests_case$par_names %in% c("Factor@A=~V1@A", "Factor@P=~V1@P")), ]
    popVals <- getSigma(return_mats = F)
    s2ests_case <- merge(s2ests_case, popVals, by = "par_names")
    s2ests_case$level <- "case"
    s2ests_case$MCSampID <- attr(s1ests, "MCSampID")
    s2ests_case$n <- attr(s1ests, "n"); s2ests_case$G <- attr(s1ests, "G")
    s2ests_case$condition <- paste0(attr(s1ests, "n"), "-", attr(s1ests, "G"))
    s2ests_case$analType <- "2SMLE"
    s2ests_case$s1priorType <- paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision"))
    s2ests_case$s1iter <- attr(s1ests, "iter")
    s2ests_case$s1mPSRF <- attr(s1ests, "mPSRF")
    s2ests_case$coverage <- s2ests_case$ci.lower < s2ests_case$pop & s2ests_case$pop < s2ests_case$ci.upper
    s2ests_case$bias <- s2ests_case$est - s2ests_case$pop
    s2ests_case$RB <- s2ests_case$bias / s2ests_case$pop
    
    s2ests_case <- s2ests_case[, c("MCSampID", "n", "G", "condition", "analType", "s1priorType", 
                         "s1iter", "s1mPSRF", "par_names", "level", 
                         "pop", "est", "se", "ci.lower", "ci.upper", "coverage", "bias", "RB")] # remove non-redundant rows and reorder
    } else {
      s2ests_case <- NULL
    }
    if (lavInspect(fit_dyad, "converged")) {
      s2ests_dyad <- parameterEstimates(fit_dyad)
      s2ests_dyad$lhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests_dyad$lhs))))
      s2ests_dyad$rhs <- gsub("_out", "@A", gsub("_in", "@P", gsub("_ij", "@AP", gsub("_ji", "@PA", s2ests_dyad$rhs))))
      s2ests_dyad$par_names <- paste0(s2ests_dyad$lhs, s2ests_dyad$op, s2ests_dyad$rhs)
      ## remove factor loadings constrained for identification and redundant rows
      s2ests_dyad <- s2ests_dyad[!(s2ests_dyad$par_names %in% c("Factor@AP=~V1@AP", "Factor@PA=~V1@PA",
                                                                "Factor@PA=~V2@PA", "Factor@PA=~V3@PA",
                                                                "Factor@PA~~Factor@PA", "V1@PA~~V1@PA",
                                                                "V2@PA~~V2@PA", "V3@PA~~V3@PA")), ]
      popVals <- getSigma(return_mats = F)
      s2ests_dyad <- merge(s2ests_dyad, popVals, by = "par_names")
      s2ests_dyad$level <- "dyad"
      s2ests_dyad$MCSampID <- attr(s1ests, "MCSampID")
      s2ests_dyad$n <- attr(s1ests, "n"); s2ests_dyad$G <- attr(s1ests, "G")
      s2ests_dyad$condition <- paste0(attr(s1ests, "n"), "-", attr(s1ests, "G"))
      s2ests_dyad$analType <- "2SMLE"
      s2ests_dyad$s1priorType <- paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision"))
      s2ests_dyad$s1iter <- attr(s1ests, "iter")
      s2ests_dyad$s1mPSRF <- attr(s1ests, "mPSRF")
      s2ests_dyad$coverage <- s2ests_dyad$ci.lower < s2ests_dyad$pop & s2ests_dyad$pop < s2ests_dyad$ci.upper
      s2ests_dyad$bias <- s2ests_dyad$est - s2ests_dyad$pop
      s2ests_dyad$RB <- s2ests_dyad$bias / s2ests_dyad$pop
      
      s2ests_dyad <- s2ests_dyad[, c("MCSampID", "n", "G", "condition", "analType", "s1priorType", 
                                     "s1iter", "s1mPSRF", "par_names", "level", 
                                     "pop", "est", "se", "ci.lower", "ci.upper", "coverage", "bias", "RB")] # remove non-redundant rows and reorder
    } else {
      s2ests_dyad <- NULL
      }
    s2ests <- rbind(s2ests_case, s2ests_dyad)
  }
  
  # MODEL FIT
  if (lavInspect(fit_combi, "converged")) {
    # model fit--combi
    standard.LRT_combi <- lavTestLRT(fit_combi, method = "standard")
    standard.stat_combi <- standard.LRT_combi$`Chisq diff`[2]
    df_combi <- standard.LRT_combi$`Df diff`[2]
    standard.p_combi <- standard.LRT_combi$`Pr(>Chisq)`[2]
    
    adf.LRT_combi <- lavTestLRT(fit_combi, type = "browne.residual.adf")
    adf.stat_combi <- adf.LRT_combi$`Chisq diff`[2]
    adf.p_combi <- adf.LRT_combi$`Pr(>Chisq)`[2]
    
    sb.LRT_combi <- lavTestLRT(fit_combi, method = "satorra.bentler.2001")
    sb.stat_combi <- sb.LRT_combi$`Chisq diff`[2]
    sb.p_combi <- sb.LRT_combi$`Pr(>Chisq)`[2]
    
    ss.LRT_combi <- lavTestLRT(fit_combi, test = "scaled.shifted")
    ss.stat_combi <- ss.LRT_combi$`Chisq diff`[2]
    ss.p_combi <- ss.LRT_combi$`Pr(>Chisq)`[2]
    
    N_combi <- attr(s1ests, "nobs")["case"] + attr(s1ests, "nobs")["dyad"]
    yb.corrected.stat_combi <- adf.stat_combi / (1 + (adf.stat_combi/N_combi))
    yb.corrected.p_combi <- pchisq(yb.corrected.stat_combi, df = df_combi, lower = F)
    
    yb.F.stat_combi <- ((N_combi - df_combi)/(N_combi * df_combi))*adf.stat_combi
    yb.F.df2_combi <- N_combi - df_combi
    yb.F.p_combi <- pf(yb.F.stat_combi, df1 = df_combi, df2 = yb.F.df2_combi, lower.tail = F)
  } else {
    standard.stat_combi <- standard.p_combi <- adf.stat_combi <- adf.p_combi <-
      sb.stat_combi <- sb.p_combi <- ss.stat_combi <- ss.p_combi <- 
      yb.corrected.stat_combi <- yb.corrected.p_combi <- yb.F.stat_combi <- 
      yb.F.p_combi <- NULL
    
    N_combi <- attr(s1ests, "nobs")["combi"]
    yb.F.df2_combi <- N_combi - df_combi
  }
  
  if (lavInspect(fit_case, "converged")) {
    # model fit--case
    standard.LRT_case <- lavTestLRT(fit_case, method = "standard")
    standard.stat_case <- standard.LRT_case$`Chisq diff`[2]
    df_case <- standard.LRT_case$`Df diff`[2]
    standard.p_case <- standard.LRT_case$`Pr(>Chisq)`[2]
    
    adf.LRT_case <- lavTestLRT(fit_case, type = "browne.residual.adf")
    adf.stat_case <- adf.LRT_case$`Chisq diff`[2]
    adf.p_case <- adf.LRT_case$`Pr(>Chisq)`[2]
    
    sb.LRT_case <- lavTestLRT(fit_case, method = "satorra.bentler.2001")
    sb.stat_case <- sb.LRT_case$`Chisq diff`[2]
    sb.p_case <- sb.LRT_case$`Pr(>Chisq)`[2]
    
    ss.LRT_case <- lavTestLRT(fit_case, test = "scaled.shifted")
    ss.stat_case <- ss.LRT_case$`Chisq diff`[2]
    ss.p_case <- ss.LRT_case$`Pr(>Chisq)`[2]
    
    N_case <- attr(s1ests, "nobs")["case"]
    yb.corrected.stat_case <- adf.stat_case / (1 + (adf.stat_case/N_case))
    yb.corrected.p_case <- pchisq(yb.corrected.stat_case, df = df_case, lower = F)
    
    yb.F.stat_case <- ((N_case - df_case)/(N_case * df_case))*adf.stat_case
    yb.F.df2_case <- N_case - df_case
    yb.F.p_case <- pf(yb.F.stat_case, df1 = df_case, df2 = yb.F.df2_case, lower.tail = F)
  } else {
    standard.stat_case <- standard.p_case <- adf.stat_case <- adf.p_case <-
      sb.stat_case <- sb.p_case <- ss.stat_case <- ss.p_case <- 
      yb.corrected.stat_case <- yb.corrected.p_case <- yb.F.stat_case <- 
      yb.F.p_case <- NULL
    
    N_case <- attr(s1ests, "nobs")["case"]
    yb.F.df2_case <- N_case - df_case
  }
  
  if (lavInspect(fit_dyad, "converged")) {
    # model fit--dyad
    standard.LRT_dyad <- lavTestLRT(fit_dyad, method = "standard")
    standard.stat_dyad <- standard.LRT_dyad$`Chisq diff`[2]
    df_dyad <- standard.LRT_dyad$`Df diff`[2]
    standard.p_dyad <- standard.LRT_dyad$`Pr(>Chisq)`[2]
    
    adf.LRT_dyad <- lavTestLRT(fit_dyad, type = "browne.residual.adf")
    adf.stat_dyad <- adf.LRT_dyad$`Chisq diff`[2]
    adf.p_dyad <- adf.LRT_dyad$`Pr(>Chisq)`[2]
    
    sb.LRT_dyad <- lavTestLRT(fit_dyad, method = "satorra.bentler.2001")
    sb.stat_dyad <- sb.LRT_dyad$`Chisq diff`[2]
    sb.p_dyad <- sb.LRT_dyad$`Pr(>Chisq)`[2]
    
    ss.LRT_dyad <- lavTestLRT(fit_dyad, test = "scaled.shifted")
    ss.stat_dyad <- ss.LRT_dyad$`Chisq diff`[2]
    ss.p_dyad <- ss.LRT_dyad$`Pr(>Chisq)`[2]
    
    N_dyad <- attr(s1ests, "nobs")["dyad"]
    yb.corrected.stat_dyad <- adf.stat_dyad / (1 + (adf.stat_dyad/N_dyad))
    yb.corrected.p_dyad <- pchisq(yb.corrected.stat_dyad, df = df_dyad, lower = F)
    
    yb.F.stat_dyad <- ((N_dyad - df_dyad)/(N_dyad * df_dyad))*adf.stat_dyad
    yb.F.df2_dyad <- N_dyad - df_dyad
    yb.F.p_dyad <- pf(yb.F.stat_dyad, df1 = df_dyad, df2 = yb.F.df2_dyad, lower.tail = F)
  } else {
    standard.stat_dyad <- standard.p_dyad <- adf.stat_dyad <- adf.p_dyad <-
      sb.stat_dyad <- sb.p_dyad <- ss.stat_dyad <- ss.p_dyad <- 
      yb.corrected.stat_dyad <- yb.corrected.p_dyad <- yb.F.stat_dyad <- 
      yb.F.p_dyad <- NULL
    
    N_dyad <- attr(s1ests, "nobs")["dyad"]
    yb.F.df2_dyad <- N_dyad - df_dyad
    
  }
  
  standard <- c(fitStat.type = "standard", 
                combi.stat = ifelse(!is.null(standard.stat_combi), standard.stat_combi, "NA"), 
                combi.df = df_combi, 
                combi.p = ifelse(!is.null(standard.p_combi), standard.p_combi, "NA"),
                case.stat = ifelse(!is.null(standard.stat_case), standard.stat_case, "NA"), 
                case.df = df_case, case.p = ifelse(!is.null(standard.p_case), standard.p_case, "NA"),
                dyad.stat = ifelse(!is.null(standard.stat_dyad), standard.stat_dyad, "NA"), 
                dyad.df = df_dyad, dyad.p = ifelse(!is.null(standard.p_dyad), standard.p_dyad, "NA")) # Standard chi-square test
  
  adf <- c(fitStat.type = "browne.residual.adf", 
                combi.stat = ifelse(!is.null(adf.stat_combi), adf.stat_combi, "NA"), 
                combi.df = df_combi, 
                combi.p = ifelse(!is.null(adf.p_combi), adf.p_combi, "NA"),
                case.stat = ifelse(!is.null(adf.stat_case), adf.stat_case, "NA"), 
                case.df = df_case, case.p = ifelse(!is.null(adf.p_case), adf.p_case, "NA"),
                dyad.stat = ifelse(!is.null(adf.stat_dyad), adf.stat_dyad, "NA"), 
                dyad.df = df_dyad, dyad.p = ifelse(!is.null(adf.p_dyad), adf.p_dyad, "NA")) # Browne's residual-based ADF
  
  sb <- c(fitStat.type = "satorra.bentler.2001", 
           combi.stat = ifelse(!is.null(sb.stat_combi), sb.stat_combi, "NA"), 
           combi.df = df_combi, 
           combi.p = ifelse(!is.null(sb.p_combi), sb.p_combi, "NA"),
           case.stat = ifelse(!is.null(sb.stat_case), sb.stat_case, "NA"), 
           case.df = df_case, case.p = ifelse(!is.null(sb.p_case), sb.p_case, "NA"),
           dyad.stat = ifelse(!is.null(sb.stat_dyad), sb.stat_dyad, "NA"), 
           dyad.df = df_dyad, dyad.p = ifelse(!is.null(sb.p_dyad), sb.p_dyad, "NA")) # Satorra and Bentler's corrected statistic
  
  ss <- c(fitStat.type = "scaled.shifted", 
           combi.stat = ifelse(!is.null(ss.stat_combi), ss.stat_combi, "NA"), 
           combi.df = df_combi, 
           combi.p = ifelse(!is.null(ss.p_combi), ss.p_combi, "NA"),
           case.stat = ifelse(!is.null(ss.stat_case), ss.stat_case, "NA"), 
           case.df = df_case, case.p = ifelse(!is.null(ss.p_case), ss.p_case, "NA"),
           dyad.stat = ifelse(!is.null(ss.stat_dyad), ss.stat_dyad, "NA"), 
           dyad.df = df_dyad, dyad.p = ifelse(!is.null(ss.p_dyad), ss.p_dyad, "NA")) # Scaled-shifted statistic
  
  yb.corrected <- c(fitStat.type = "yuan.bentler.corrected.adf", 
           combi.stat = ifelse(!is.null(yb.corrected.stat_combi), yb.corrected.stat_combi, "NA"), 
           combi.df = df_combi, 
           combi.p = ifelse(!is.null(yb.corrected.p_combi), yb.corrected.p_combi, "NA"),
           case.stat = ifelse(!is.null(yb.corrected.stat_case), yb.corrected.stat_case, "NA"), 
           case.df = df_case, case.p = ifelse(!is.null(yb.corrected.p_case), yb.corrected.p_case, "NA"),
           dyad.stat = ifelse(!is.null(yb.corrected.stat_dyad), yb.corrected.stat_dyad, "NA"), 
           dyad.df = df_dyad, dyad.p = ifelse(!is.null(yb.corrected.p_dyad), yb.corrected.p_dyad, "NA")) # Yuan and Bentler's small-sample correction for ADF
  
  yb.F <- c(fitStat.type = "yuan.bentler.F",
           combi.stat = ifelse(!is.null(yb.F.stat_combi), yb.F.stat_combi, "NA"), 
           combi.df = paste0(df_combi, ",", yb.F.df2_combi), 
           combi.p = ifelse(!is.null(yb.F.p_combi), yb.F.p_combi, "NA"),
           case.stat = ifelse(!is.null(yb.F.stat_case), yb.F.stat_case, "NA"), 
           case.df = paste0(df_case, ",", yb.F.df2_case), 
           case.p = ifelse(!is.null(yb.F.p_case), yb.F.p_case, "NA"),
           dyad.stat = ifelse(!is.null(yb.F.stat_dyad), yb.F.stat_dyad, "NA"), 
           dyad.df = paste0(df_dyad, ",", yb.F.df2_dyad), 
           dyad.p = ifelse(!is.null(yb.F.p_dyad), yb.F.p_dyad, "NA"))
  
  s2mod <- data.frame(cbind(MCSampID = attr(s1ests, "MCSampID"),
                            n = attr(s1ests, "n"),
                            G = attr(s1ests, "G"),
                            condition = paste0(attr(s1ests, "n"), "-", attr(s1ests, "G")),
                            analType = "2SMLE",
                            s1priorType = paste0("MCMC-", attr(s1ests, "priorType"), "-", attr(s1ests, "precision")),
                            s1iter = attr(s1ests, "iter"),
                            s1mPSRF = attr(s1ests, "mPSRF"), rbind(standard, adf, sb, ss, yb.corrected, yb.F)))
  
  s2result <- list(s2ests = s2ests, s2mod = s2mod)
  
  if (savefile) saveRDS(s2result, file = paste0("s2_ID", attr(s1ests, "MCSampID"),
                                                ".nG", attr(s1ests, "G"), ".n", 
                                                attr(s1ests, "n"), "_", attr(s1ests, "priorType"),
                                                "_", attr(s1ests, "precision"), ".rds"))
   return(s2result)
}

# lavs2(s1ests = s1ests, savefile = F)

#----

