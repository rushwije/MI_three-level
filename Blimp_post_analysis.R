#######################################################################################
#   Filename    :	Results_Tables.R 
#                   
#   Project     :   This script is to be used for post imputation analysis and pooling for 
#                   SMC-JM-3L in Blimp and accompanies the 
#                   BiomJ article "Evaluation of approaches for accommodating interactions and 
#                   non-linear terms in multiple imputation of incomplete three-level data"  
#                   [Tables S2-S13]
#   Author      :   Rushani Wijesuriya                                                                
#   Date        :   07.09.2020
########################################################################################

#************************************Notes for running this script**********************

# 1. Insert the path to the folder where the imputations are saved
# 2. Results will be saved in the respective folder under each simulation scenario


# Required packages
library(mitml)  #For fitting the lmer model
library(lme4)  # For pooling the results
library(writexl)  #For writing/exporting the results

## Clear the workspace
rm(list = ls())

## Folder names for changing the directory
Analysis <- c("Interaction", "Interaction2", "Squared")
MD.mech <- c("MAR-CATS", "MAR-inflated")
Clus <- c("Clus40", "Clus10")

for (j in 1:length(Analysis)) {
    
    for (k in 1:length(MD.mech)) {
        
        for (l in 1:length(Clus)) {
            
            # Set working directory to where the imputations are saved
            setwd(file.path("...", Analysis[[j]], MD.mech[[k]], Clus[[l]], "/blimpimps"))
            
            ## Load the simulated data
            temp <- list.files(pattern = glob2rx("*.csv"))
            data <- lapply(temp, read.csv, header = FALSE)
            
            ## Generate empty matrices to save the results
            blimp_results.est <- matrix(NA, nrow = 7, ncol = length(temp))
            blimp_results.sd <- matrix(NA, nrow = 7, ncol = length(temp))
            blimp_results.RE <- matrix(NA, nrow = 3, ncol = length(temp))
            blimp_results.ICC <- matrix(NA, nrow = 2, ncol = length(temp))
            blimp_results.CI <- c()
            
            
            ## Analysis
            for (i in 1:length(temp)) {
                
                impdata <- data[[i]]
                names(impdata) <- c("impno", "c_id", "time", "napscore_z", "school", 
                  "c_age", "c_gender", "c_ses", "c_nap1_z", "prev_dep", "prev_sdq")
                
                impdata$c_id <- as.factor(impdata$c_id)
                
                mylist <- list()
                
                for (m in 1:20) {
                  mylist[[m]] <- impdata[impdata$impno == m, ]
                }
                
                
                # Fit the analysis of interest on the imputed datasets
                if (j == 1) {
                  mods <- lapply(mylist, function(d) {
                    lmer(napscore_z ~ prev_dep + time + prev_dep * time + c_age + 
                      as.factor(c_gender) + c_nap1_z + c_ses + (1 | school/c_id), 
                      data = d)
                  })
                } else if (j == 2) {
                  mods <- lapply(mylist, function(d) {
                    lmer(napscore_z ~ prev_dep + time + prev_dep * c_ses + c_age + 
                      as.factor(c_gender) + c_nap1_z + c_ses + (1 | school/c_id), 
                      data = d)
                  })
                } else {
                  mods <- lapply(mylist, function(d) {
                    d$prev_dep2 <- d$prev_dep * d$prev_dep
                    lmer(napscore_z ~ prev_dep + time + c_age + as.factor(c_gender) + 
                      c_nap1_z + c_ses + prev_dep2 + (1 | school/c_id), data = d)
                  })
                }
                
                # Combine the estimates
                MI_est <- testEstimates(mods, var.comp = TRUE, df.com = NULL)
                
                CI <- matrix(NA, 7, 2)
                Conf <- confint(MI_est)
                CI[, 1] <- as.vector(Conf[2:8, 1])
                CI[, 2] <- as.vector(Conf[2:8, 2])
                blimp_results.CI <- cbind(blimp_results.CI, CI)
                
                # Store the estimates
                blimp_results.est[, i] <- MI_est$estimates[2:8, 1]
                blimp_results.sd[, i] <- MI_est$estimates[2:8, 2]
                blimp_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
                blimp_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
                blimp_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
                
                blimp_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 
                  1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
                blimp_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
                  1])/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 
                  1])
                
                
            }
            
            rownames(blimp_results.est) <- rownames(blimp_results.sd) <- rownames(blimp_results.CI) <- c("prev_dep", 
                "time", "c_age", "c_gender", "c_nap1_z", "c_ses", "inter")
            colnames(blimp_results.est) <- colnames(blimp_results.sd) <- colnames(blimp_results.ICC) <- colnames(blimp_results.RE) <- c(seq(1:length(temp)))
            rownames(blimp_results.RE) <- c("level 3", "level 2", "level 1")
            rownames(blimp_results.ICC) <- c("level 3", "level 2")
            colnames(blimp_results.CI) <- c(rep(1:length(temp), each = 2))
            
            ## Save the results
            write_xlsx(as.data.frame(cbind(rownames(blimp_results.est), blimp_results.est)), 
                "BLIMP_results.est.xlsx")
            write_xlsx(as.data.frame(cbind(rownames(blimp_results.sd), blimp_results.sd)), 
                "BLIMP_results.sd.xlsx")
            write_xlsx(as.data.frame(cbind(rownames(blimp_results.RE), blimp_results.RE)), 
                "BLIMP_results.RE.xlsx")
            write_xlsx(as.data.frame(cbind(rownames(blimp_results.ICC), blimp_results.ICC)), 
                "BLIMP_results.ICC.xlsx")
            write_xlsx(as.data.frame(cbind(rownames(blimp_results.CI), blimp_results.CI)), 
                "BLIMP_results.CI.xlsx")
        }
    }
}
