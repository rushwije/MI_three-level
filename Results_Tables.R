#######################################################################################
#   Filename    :	Results_Tables.R 
#                   
#   Project     :   This script will generate the tables in the supplementary section of
#                   BiomJ article "Evaluation of approaches for accommodating interactions and 
#                   non-linear terms in multiple imputation of incomplete three-level data"  
#                   [Tables S2-S14]
#   Author      :   Rushani Wijesuriya                                                                
#   Date        :   07.09.2020
########################################################################################

#************************************Notes for running this script**********************

# 1. Insert the path to the  Code_and_Data folder in line 65 
# 2. Run the script from begining to end, the tables will be saved  in the clus10 folder under 
#    each missing data mechanism under each analysis model
#    For example tables S2 and S3 in folder Interaction-> MAR-CATS -> Clus10 and tables S4 and S5 
#    in the folder Interaction -> MAR-inflated -> Clus10 etc.


rm(list = ls())

options(scipen=999) # To disable scientific notation

library(rsimsum)  #To summarize results from the simulation stduy
library(readxl)  # To read in the excel files
library(WriteXLS)  # To export results 
library(gt)  #To generate the gt table
library(tidyverse)  #To generate the gt table
library(glue)  #To generate the gt table
library(tibble)  #To generate the gt table

## Folder names for changing the directory
Analysis <- c("Interaction", "Interaction2", "Squared")
MD.mech <- c("MAR-CATS", "MAR-inflated")
Clus <- c("Clus40", "Clus10")

## Variance component parameters
VC1 <- 0.5
VC2 <- 0.45
VC3 <- 0.05

# MI method names under each analysis
Method.list <- list(c("JM-1L-DI-wide", "FCS-1L-DI-wide", "JM-2L-wide", "FCS-2L-wide", 
    "SMC-JM-2L-DI", "SMC-SM-2L-DI", "SMC-JM-3L"), c("JM-1L-DI-wide", "JM-1L-DI-wide-JAV", 
    "FCS-1L-DI-wide", "FCS-1L-DI-wide-Passive_c", "FCS-1L-DI-wide-Passive_all", "JM-2L-wide-JAV", 
    "FCS-2L-wide-Passive_c", "FCS-2L-wide-Passive_all", "SMC-JM-2L-DI", "SMC-SM-2L-DI", 
    "SMC-JM-3L"), c("JM-1L-DI-wide", "JM-1L-DI-wide-JAV", "FCS-1L-DI-wide", "FCS-1L-DI-wide-Passive", 
    "JM-2L-wide-JAV", "FCS-2L-wide-Passive", "SMC-JM-2L-DI", "SMC-SM-2L-DI", "SMC-JM-3L"))
# True parameters
beta1 <- c(-0.07, -0.024, -0.024)
beta3 <- c(0.013, 0.023, -0.009)

# To save results from the regression coefficient estimates
clus.list.est <- list()
# To save results from the variance component estimates
clus.list.vc <- list()

for (k in 1:length(Analysis)) {
    
    for (i in 1:length(MD.mech)) {
        
        for (j in 1:length(Clus)) {
            
            
            ## Set path to the working directory with results (See README file)
            setwd(file.path("...", "/Code_and_Data/Intermediate_results/", 
                Analysis[[k]], MD.mech[[i]], Clus[[j]]))
            
            ## Import all the results files in the folder
            temp <- list.files(pattern = glob2rx("*.xlsx"))
            
            ############################################### Regression coefficient estimates ###
            
            ## Extract the results for the regression coefficient estimates only
            temp1 <- temp[grep("est", temp)]
            
            if (k == 1) {
                dup1 <- c("MVNIslwide_results.est.xlsx", "FCSslwide_results.est.xlsx", 
                  "MVNImlJAV_results.est.xlsx", "FCSmlJAV_results.est.xlsx", "MVNImlDI_results.est.xlsx", 
                  "MDMB_results.est.xlsx", "BLIMP_results.est.xlsx")
            } else if (k == 2) {
                dup1 <- c("MVNIslwide_results.est.xlsx", "MVNIslwide_results_JAV.est.xlsx", 
                  "FCSslwide_results.est.xlsx", "FCSslwidepassive1_results.est.xlsx", 
                  "FCSslwidepassive2_results.est.xlsx", "MVNImlJAV_results_JAV.est.xlsx", 
                  "FCSmlJAV_resultsP1.est.xlsx", "FCSmlJAV_resultsP2.est.xlsx", "MVNImlDI_results.est.xlsx", 
                  "MDMB_results.est.xlsx", "BLIMP_results.est.xlsx")
            } else {
                dup1 <- c("MVNIslwide_results.est.xlsx", "MVNIslwide_results.sq.est.xlsx", 
                  "FCSslwide_results.est.xlsx", "FCSslwide_results.passive.est.xlsx", 
                  "MVNImlJAV_results.sq.est.xlsx", "FCSmlJAV_results.passive.est.xlsx", 
                  "MVNImlDI_results.est.xlsx", "MDMB_results.est.xlsx", "BLIMP_results.est.xlsx")
            }
            
            # Order the extracted results list
            temp1 <- temp1[order(match(temp1, dup1))]
            
            # Create empty matrix to hold values of the estimates and compute perfromance
            # measures
            MSIM_dep <- matrix(NA, nrow = length(temp1) * 1000, ncol = 4)
            MSIM_inter <- matrix(NA, nrow = length(temp1) * 1000, ncol = 4)
            
            # Extract the main effect and interaction effect estimates
            dep <- c()
            inter <- c()
            
            for (m in 1:length(temp1)) {
                data <- read_excel(temp1[m])
                data <- data[-1]
                dep <- c(dep, as.vector(as.matrix(data[1, ])))
                inter <- c(inter, as.vector(as.matrix(data[7, ])))
            }
            
            MSIM_dep[, 3] <- dep
            MSIM_inter[, 3] <- inter
            
            ## Extract the results for the standard error estimates only
            temp2 <- temp[grep("sd", temp)]
            
            if (k == 1) {
                dup2 <- c("MVNIslwide_results.sd.xlsx", "FCSslwide_results.sd.xlsx", 
                  "MVNImlJAV_results.sd.xlsx", "FCSmlJAV_results.sd.xlsx", "MVNImlDI_results.sd.xlsx", 
                  "MDMB_results.sd.xlsx", "BLIMP_results.sd.xlsx")
            } else if (k == 2) {
                dup2 <- c("MVNIslwide_results.sd.xlsx", "MVNIslwide_results_JAV.sd.xlsx", 
                  "FCSslwide_results.sd.xlsx", "FCSslwidepassive1_results.sd.xlsx", 
                  "FCSslwidepassive2_results.sd.xlsx", "MVNImlJAV_results_JAV.sd.xlsx", 
                  "FCSmlJAV_resultsP1.sd.xlsx", "FCSmlJAV_resultsP2.sd.xlsx", "MVNImlDI_results.sd.xlsx", 
                  "MDMB_results.sd.xlsx", "BLIMP_results.sd.xlsx")
            } else {
                dup2 <- c("MVNIslwide_results.sd.xlsx", "MVNIslwide_results.sq.sd.xlsx", 
                  "FCSslwide_results.sd.xlsx", "FCSslwide_results.passive.sd.xlsx", 
                  "MVNImlJAV_results.sq.sd.xlsx", "FCSmlJAV_results.passive.sd.xlsx", 
                  "MVNImlDI_results.sd.xlsx", "MDMB_results.sd.xlsx", "BLIMP_results.sd.xlsx")
            }
            
            # Order the extracted results list
            temp2 <- temp2[order(match(temp2, dup2))]
            
            # Extract the standard error estimates of the main effect and the interaction
            # effect
            sd.dep <- c()
            sd.inter <- c()
            
            for (m in 1:length(temp2)) {
                data <- read_excel(temp2[m])
                data <- data[-1]
                sd.dep <- c(sd.dep, as.vector(as.matrix(data[1, ])))
                sd.inter <- c(sd.inter, as.vector(as.matrix(data[7, ])))
            }
            
            MSIM_dep[, 4] <- sd.dep
            MSIM_inter[, 4] <- sd.inter
            
            MSIM_dep[, 2] <- MSIM_inter[, 2] <- rep(Method.list[[k]], each = 1000)
            
            MSIM_dep[, 1] <- MSIM_inter[, 1] <- rep(1:1000, length(temp1))
            
            colnames(MSIM_dep) <- colnames(MSIM_inter) <- c("dataset", "Method", 
                "b", "se")
            
            MSIM_dep <- as.data.frame(MSIM_dep)
            MSIM_inter <- as.data.frame(MSIM_inter)
            
            MSIM_dep$dataset <- as.numeric(levels(MSIM_dep$dataset))[MSIM_dep$dataset]
            MSIM_dep$b <- as.numeric(levels(MSIM_dep$b))[MSIM_dep$b]
            MSIM_dep$se <- as.numeric(levels(MSIM_dep$se))[MSIM_dep$se]
            
            MSIM_inter$dataset <- as.numeric(levels(MSIM_inter$dataset))[MSIM_inter$dataset]
            MSIM_inter$b <- as.numeric(levels(MSIM_inter$b))[MSIM_inter$b]
            MSIM_inter$se <- as.numeric(levels(MSIM_inter$se))[MSIM_inter$se]
            
            MSIM_dep <- MSIM_dep[order(MSIM_dep$dataset), ]
            MSIM_inter <- MSIM_inter[order(MSIM_inter$dataset), ]
            
            ## Compute the performance measures
            
            ## For the main effect estimate
            S1 <- simsum(data = MSIM_dep, estvarname = "b", true = beta1[k], se = "se", 
                methodvar = "Method", ref = "JM-1L-DI-wide")
            
            ## Specify the order in which methods will be presented in tables
            if (k == 1) {
                od <- c(1, 2, 4, 3, 5, 7, 6)
            } else if (k == 2) {
                od <- c(1, 3, 5, 4, 8, 7, 2, 6, 9, 11, 10)
            } else {
                od <- c(1, 3, 4, 6, 2, 5, 7, 9, 8)
            }
            
            # 1. Average estimate
            est <- cbind(get_data(S1, stats = c("thetamean")), order = od)
            
            # 2. Bias
            bias <- cbind(get_data(S1, stats = c("bias")), order = od)
            
            # 3. Relative bias
            RB <- cbind((get_data(S1, stats = c("bias"))[2]/beta1[k]) * 100, order = od)
            
            # 4. Empirical SE
            emp <- cbind(get_data(S1, stats = c("empse")), order = od)
            
            # 5. Model based SE
            mod <- cbind(get_data(S1, stats = c("modelse")), order = od)
            
            # 6. Coverage
            cov <- cbind(get_data(S1, stats = c("cover")), order = od)
            
            prevdeplist <- list(est, bias, RB, emp, mod, cov)
            sorted.est1 <- lapply(prevdeplist, function(d) {
                d[order(d$order), ]
            })
            
            
            ## For the interaction effect estimate
            S2 <- simsum(data = MSIM_inter, estvarname = "b", true = beta3[k], se = "se", 
                methodvar = "Method", ref = "JM-1L-DI-wide")
            
            # 1. Average estimate
            est1 <- cbind(get_data(S2, stats = c("thetamean")), order = od)
            
            # 2. Bias
            bias1 <- cbind(get_data(S2, stats = c("bias")), order = od)
            
            # 3. Relative bias
            RB1 <- cbind((get_data(S2, stats = c("bias"))[2]/beta3[k]) * 100, order = od)
            
            # 4. Empirical SE
            emp1 <- cbind(get_data(S2, stats = c("empse")), order = od)
            
            # 5. Model based SE
            mod1 <- cbind(get_data(S2, stats = c("modelse")), order = od)
            
            # 6. Coverage
            cov1 <- cbind(get_data(S2, stats = c("cover")), order = od)
            
            interlist <- list(est1, bias1, RB1, emp1, mod1, cov1)
            sorted.est2 <- lapply(interlist, function(d) {
                d[order(d$order), ]
            })
            
            
            results_est <- matrix(NA, nrow = length(temp1), ncol = 13)
            
            results_est[, 1] <- Method.list[[k]]
            
            results_est[, 2] <- round(sorted.est1[[1]][, 2], 3)
            results_est[, 3] <- round(sorted.est1[[2]][, 2], 4)
            results_est[, 4] <- round(sorted.est1[[3]][, 1], 1)
            results_est[, 5] <- round(sorted.est1[[4]][, 2], 4)
            results_est[, 6] <- round(sorted.est1[[5]][, 2], 4)
            results_est[, 7] <- round(sorted.est1[[6]][, 2], 6)
            results_est[, 8] <- round(sorted.est2[[1]][, 2], 3)
            results_est[, 9] <- round(sorted.est2[[2]][, 2], 4)
            results_est[, 10] <- round(sorted.est2[[3]][, 1], 1)
            results_est[, 11] <- round(sorted.est2[[4]][, 2], 4)
            results_est[, 12] <- round(sorted.est2[[5]][, 2], 4)
            results_est[, 13] <- round(sorted.est2[[6]][, 2], 3)
            
            
            colnames(results_est) <- c("Method", "Main effect.Average estimate", 
                "Main effect.Bias", "Main effect.Relative Bias(%)", "Main effect.Emp SE", 
                "Main effect.Model SE", "Main effect.coverage", "Interaction effect.Average estimate", 
                "Interaction effect.Bias", "Interaction effect.Relative Bias(%)", 
                "Interaction effect.Emp SE", "Interaction effect.Model SE", "Interaction effect.coverage")
            
            clus.list.est[[j]] <- results_est
            
            ############################################### Variance component (VC) estimates ###
            
            ## Extract the results for the estimates only
            temp3 <- temp[grep("RE", temp)]
            
            if (k == 1) {
                dup3 <- c("MVNIslwide_results.RE.xlsx", "FCSslwide_results.RE.xlsx", 
                  "MVNImlJAV_results.RE.xlsx", "FCSmlJAV_results.RE.xlsx", "MVNImlDI_results.RE.xlsx", 
                  "MDMB_results.RE.xlsx", "BLIMP_results.RE.xlsx")
            } else if (k == 2) {
                dup3 <- c("MVNIslwide_results.RE.xlsx", "MVNIslwide_results_JAV.RE.xlsx", 
                  "FCSslwide_results.RE.xlsx", "FCSslwidepassive1_results.RE.xlsx", 
                  "FCSslwidepassive2_results.RE.xlsx", "MVNImlJAV_results_JAV.RE.xlsx", 
                  "FCSmlJAV_resultsP1.RE.xlsx", "FCSmlJAV_resultsP2.RE.xlsx", "MVNImlDI_results.RE.xlsx", 
                  "MDMB_results.RE.xlsx", "BLIMP_results.RE.xlsx")
            } else {
                dup3 <- c("MVNIslwide_results.RE.xlsx", "MVNIslwide_results.sq.RE.xlsx", 
                  "FCSslwide_results.RE.xlsx", "FCSslwide_results.passive.RE.xlsx", 
                  "MVNImlJAV_results.sq.RE.xlsx", "FCSmlJAV_results.passive.RE.xlsx", 
                  "MVNImlDI_results.RE.xlsx", "MDMB_results.RE.xlsx", "BLIMP_results.RE.xlsx")
            }
            
            temp3 <- temp3[order(match(temp3, dup3))]
            
            MSIM_L1 <- MSIM_L2 <- MSIM_L3 <- matrix(NA, nrow = length(temp1) * 1000, 
                ncol = 4)
            
            L1 <- L2 <- L3 <- c()
            
            for (m in 1:length(temp3)) {
                data <- read_excel(temp3[m], col_types = "numeric")
                data <- data[-1]
                L1 <- c(L1, as.vector(as.matrix(data[3, ])))
                L2 <- c(L2, as.vector(as.matrix(data[2, ])))
                L3 <- c(L3, as.vector(as.matrix(data[1, ])))
                
            }
            
            MSIM_L1[, 3] <- (L1 * L1)
            MSIM_L2[, 3] <- (L2 * L2)
            MSIM_L3[, 3] <- (L3 * L3)
            
            MSIM_L1[, 2] <- MSIM_L2[, 2] <- MSIM_L3[, 2] <- rep(Method.list[[k]], 
                each = 1000)
            
            MSIM_L1[, 1] <- MSIM_L1[, 4] <- MSIM_L2[, 1] <- MSIM_L2[, 4] <- MSIM_L3[, 
                1] <- MSIM_L3[, 4] <- rep(1:1000, length(temp1))
            
            colnames(MSIM_L1) <- colnames(MSIM_L2) <- colnames(MSIM_L3) <- c("dataset", 
                "Method", "b", "se")
            
            MSIM_L1 <- as.data.frame(MSIM_L1)
            MSIM_L2 <- as.data.frame(MSIM_L2)
            MSIM_L3 <- as.data.frame(MSIM_L3)
            
            MSIM_L1$dataset <- as.numeric(levels(MSIM_L1$dataset))[MSIM_L1$dataset]
            MSIM_L1$b <- as.numeric(levels(MSIM_L1$b))[MSIM_L1$b]
            MSIM_L1$se <- as.numeric(levels(MSIM_L1$se))[MSIM_L1$se]
            
            MSIM_L2$dataset <- as.numeric(levels(MSIM_L2$dataset))[MSIM_L2$dataset]
            MSIM_L2$b <- as.numeric(levels(MSIM_L2$b))[MSIM_L2$b]
            MSIM_L2$se <- as.numeric(levels(MSIM_L2$se))[MSIM_L2$se]
            
            MSIM_L3$dataset <- as.numeric(levels(MSIM_L3$dataset))[MSIM_L3$dataset]
            MSIM_L3$b <- as.numeric(levels(MSIM_L3$b))[MSIM_L3$b]
            MSIM_L3$se <- as.numeric(levels(MSIM_L3$se))[MSIM_L3$se]
            
            MSIM_L1 <- MSIM_L1[order(MSIM_L1$dataset), ]
            MSIM_L2 <- MSIM_L2[order(MSIM_L2$dataset), ]
            MSIM_L3 <- MSIM_L3[order(MSIM_L3$dataset), ]
            
            ############################################### Compute the performance measures ###
            
            ## For the level 1 VC
            S1 <- simsum(data = MSIM_L1, estvarname = "b", true = VC1, se = "se", 
                methodvar = "Method", ref = "JM-1L-DI-wide")
            
            # 1. Bias
            bias <- cbind(get_data(S1, stats = c("bias")), order = od)
            
            # 2. Relative bias
            RB <- cbind((get_data(S1, stats = c("bias"))[2]/VC1) * 100, order = od)
            
            # 3. Empirical standard error
            emp <- cbind(get_data(S1, stats = c("empse")), order = od)
            
            l1list <- list(bias, RB, emp)
            
            sorted.VC1 <- lapply(l1list, function(d) {
                d[order(d$order), ]
            })
            
            ## For the level 2 VC
            S2 <- simsum(data = MSIM_L2, estvarname = "b", true = VC2, se = "se", 
                methodvar = "Method", ref = "JM-1L-DI-wide")
            
            # 1. Bias
            bias2 <- cbind(get_data(S2, stats = c("bias")), order = od)
            
            # 2. Relative bias
            RB2 <- cbind((get_data(S2, stats = c("bias"))[2]/VC2) * 100, order = od)
            
            # 3. Empirical standard error
            emp2 <- cbind(get_data(S2, stats = c("empse")), order = od)
            
            
            l2list <- list(bias2, RB2, emp2)
            
            sorted.VC2 <- lapply(l2list, function(d) {
                d[order(d$order), ]
            })
            
            ## For the level 3 VC
            S3 <- simsum(data = MSIM_L3, estvarname = "b", true = VC3, se = "se", 
                methodvar = "Method", ref = "JM-1L-DI-wide")
            
            # 1. Bias
            bias3 <- cbind(get_data(S3, stats = c("bias")), order = od)
            
            # 2. Relative bias
            RB3 <- cbind((get_data(S3, stats = c("bias"))[2]/VC3) * 100, order = od)
            
            # 3. Empirical standard error
            emp3 <- cbind(get_data(S3, stats = c("empse")), order = od)
            
            l3list <- list(bias3, RB3, emp3)
            
            sorted.VC3 <- lapply(l3list, function(d) {
                d[order(d$order), ]
            })
            
            results_VC <- matrix(NA, length(temp1), ncol = 10)
            
            results_VC[, 1] <- Method.list[[k]]
            results_VC[, 2] <- round(sorted.VC3[[1]][, 2], 4)
            results_VC[, 3] <- round(sorted.VC3[[2]][, 1], 2)
            results_VC[, 4] <- round(sorted.VC3[[3]][, 2], 3)
            results_VC[, 5] <- round(sorted.VC2[[1]][, 2], 4)
            results_VC[, 6] <- round(sorted.VC2[[2]][, 1], 2)
            results_VC[, 7] <- round(sorted.VC2[[3]][, 2], 3)
            results_VC[, 8] <- round(sorted.VC1[[1]][, 2], 4)
            results_VC[, 9] <- round(sorted.VC1[[2]][, 1], 2)
            results_VC[, 10] <- round(sorted.VC1[[3]][, 2], 3)
            
            
            results_VC <- as.data.frame(results_VC)
            colnames(results_VC) <- c("Method", "Level 3 VC.Bias", "Level 3 VC.Relative Bias(%)", 
                "Level 3 VC.Emp SE", "Level 2 VC.Bias", "Level 2 VC.Relative Bias(%)", 
                "Level 2 VC.Emp SE", "Level 1 VC.Bias", "Level 1 VC.Relative Bias(%)", 
                "Level 1 VC.Emp SE")
            
            clus.list.vc[[j]] <- results_VC
            
        }
        
        # Save results(will be saved in the clus10 folder under each missing data mechanism under
        #each analysis model
        
        # 1. Performance measures for estimating regression coefficients
        Table_est <- rbind(clus.list.est[[1]], clus.list.est[[2]])
        Table_est <- as.data.frame(Table_est)
        Table_est$Cluster <- rep(c("40 clusters", "10 clusters"), each = length(temp1))
        write.csv(Table_est, paste0("Results_est", MD.mech[[i]], ".csv"))
        
        # Or as a gt table
        Table_est <- as_tibble(Table_est)
        Table_est <- Table_est %>% gt(groupname_col = "Cluster") %>% tab_spanner_delim(delim = ".")
        gtsave(Table_est, paste0("Results_est", MD.mech[[i]], ".pdf"))
        
        
        # 2. Performance measures for estimating VCs
        
        Table_VC <- rbind(clus.list.vc[[1]], clus.list.vc[[2]])
        Table_VC <- as.data.frame(Table_VC)
        Table_VC$Cluster <- rep(c("40 clusters", "10 clusters"), each = length(temp1))
        write.csv(Table_VC, paste0("Results_VC", MD.mech[[i]], ".csv"))
        
        # Or as a gt table
        Table_VC <- as_tibble(Table_VC)
        Table_VC <- Table_VC %>% gt(groupname_col = "Cluster") %>% tab_spanner_delim(delim = ".")
        gtsave(Table_VC, paste0("Table_VC", MD.mech[[i]], ".pdf"))
        
    }
    
}

