#######################################################################################
#   Filename    :	Blimp_prep.R 
#                   
#   Project     :   This script is to be used for recoding the data prior to running the 
#                   simulations on the Blimp application and accomponies the
#                   BiomJ article "Evaluation of approaches for accommodating interactions and 
#                   non-linear terms in multiple imputation of incomplete three-level data" 
#                   
#   Author      :   Rushani Wijesuriya                                                                
#   Date        :   07.09.2020
########################################################################################

#************************************Notes for running this script**********************

# 1. The input file path in this script should be changed to the folder where the orginal data under 
# each simulation  (line 36)
# 2. The output file path should be changed to a new folder to save the newly generated datasets 
# after transformation (line 43)

rm(list = ls())

library(DataCombine)  # To generate the lagged variables

## Folder names for changing the directory
Analysis <- c("Interaction", "Interaction2", "Squared")
MD.mech <- c("MAR-CATS", "MAR-inflated")
Clus <- c("Clus40", "Clus10")

for (j in 1:length(Analysis)) {
    
    for (k in 1:length(MD.mech)) {
        
        for (l in 1:length(Clus)) {
            
            ## Specify the input path
            setwd(file.path("...", Analysis[[j]], MD.mech[[k]], Clus[[l]], "/Data"))
            
            ## Load the orginal data
            data <- lapply(list.files(pattern = glob2rx("*.csv")), read.csv)

            
            ## Specify the output path
            setwd(file.path("...", Analysis[[j]], MD.mech[[k]], Clus[[l]], "/Blimpdata"))

            ## Transform the data
            for (i in 1:length(data)) {
                
                simdataL <- data[[i]]
                print(i)
                
                simdataL$c_id <- as.numeric(simdataL$c_id)
                
                # Create previous wave depression
                simdataL <- slide(simdataL, Var = "c_dep", GroupVar = "c_id", slideBy = -1)
                colnames(simdataL)[colnames(simdataL) == "c_dep-1"] <- "prev_dep"
                simdataL$c_dep <- NULL
                
                
                # Create previous wave SDQ variable
                simdataL <- slide(simdataL, Var = "p_sdq", GroupVar = "c_id", slideBy = -1)
                colnames(simdataL)[colnames(simdataL) == "p_sdq-1"] <- "prev_sdq"
                simdataL$p_sdq <- NULL
                
                # Remove unwanted waves
                simdataL <- subset(simdataL, wave != 2 & wave != 4 & wave != 6)
                
                # Create dummy indicators for the nominal/ordinal variables
                simdataL$c_gender <- ifelse(simdataL$c_gender == "male", 1, 0)
                simdataL$child <- NULL
                
                # Recode missing values
                simdataL$prev_dep[is.na(simdataL$prev_dep)] <- 999
                simdataL$c_ses[is.na(simdataL$c_ses)] <- 999
                
                x <- matrix(NA, ncol = ncol(simdataL), nrow = nrow(simdataL))
                
                for (m in 1:ncol(simdataL)) x[, m] <- simdataL[, m]
                
                write.table(x, paste0("data", i, ".csv"), sep = ",", col.names = FALSE,
                row.names = FALSE, quote = FALSE)

            }
        
        }
    }
}
