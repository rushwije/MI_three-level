#######################################################################################
#   Filename    :	MI_methods_analysis1.R 
#                   
#                   MI methods considered under the analysis model involving a quadratic effect
#                   of the exposure
#                                                                             
#   Project     :   BiomJ article "Evaluation of approaches for accommodating interactions and 
#                   non-linear terms in multiple imputation of incomplete three-level data"                                                             
#   Author      :   Rushani Wijesuriya                                                                
#   Date        :   07.09.2020
########################################################################################

#********Notes for running the MI approaches in this script**********************

#For the following MI approaches fractions of all the repititions were run in parallel for each 
#simulation scenario (number of higher-level clusters x missing data mechanism) for efficiency 

# 1.JM-1L-DI-wide  (5 sets of 1000/5 repetitions)
# 2.JM-1L-DI-wide-JAV (50 sets of 1000/50 repetitions)
# 3.JM-2L-wide-JAV (5 sets of 1000/5 repetitions)
# 4.SMC-JM-2L-DI  (50 sets of 1000/50 repetitions)
# 5.SMC-SM-2L-DI  (10 sets of 1000/10 repetitions)

#Individual seeds used for each simulation scenario and the paralley run simulations are given 
#below as comments under each MI approach 

#********************************************************************************
## clear the workspace
rm(list = ls())

## load the required packages
require(lme4)  #for fitting the lmer model
require(jomo)  # for single-level and two-level JM and two-level SMC JM
require(mitml)  # for pooling the results
require(mice)  #for single and two-level FCS
require(xlsx)  #for writing/exporting the results
require(DataCombine)  #for generating the lagged wave variables
require(mdmb)  #for two-level SMC-SM

# Change the working directory to the relevant folder where the data are stored
setwd("C:\\temp")

## Load the simulated data
data <- lapply(list.files(pattern = glob2rx("*.csv")), read.csv)

###############################################
###              FCS-1L-DI-wide             ###
###############################################

set.seed(99231) #MAR-CATS, 40 school clsuters
# set.seed(90271570) #MAR-CATS, 10 school clusters
# set.seed(8293216) #MAR-inflated, 40 school clusters
# set.seed(90277770)  #MAR-inflated, 10 school clusters

## Generate empty matrices to save the results
FCSslwide_results.est <- matrix(NA, nrow = 7, ncol = length(data))
FCSslwide_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
FCSslwide_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
FCSslwide_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
FCSslwide_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$child), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Create school dummy indicators
    simdataw$school <- as.factor(simdataw$school)
    
    ## Remove unwanted variables
    simdataw1 <- simdataw[, !names(simdataw) %in% c("c_id", "child", "napscore_z.2", 
        "napscore_z.4", "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", 
        "c_dep.5", "c_dep.7")]
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    
    # Specify the predictor matrix
    pred <- make.predictorMatrix(simdataw1)
    pred[!rownames(pred) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_ses"), ] <- 0
    diag(pred) <- 0
    
    ## Specify the imputation method for incomplete variables
    meth <- rep("", ncol(pred))
    meth[substr(row.names(pred), 1, 5) %in% c("c_dep")] <- "norm"
    meth[substr(row.names(pred), 1, 5) %in% c("c_ses")] <- "norm"
    
    ## Perform imputations
    imp2 <- mice(data = simdataw1, m = M, predictorMatrix = pred, method = meth, 
        maxit = 10)

    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        # 1. extract the mth imputed dataset
        datw <- complete(imp2, m)
        
        
        # 2. Remove unwanted variables
        datw <- datw[, names(datw) %in% c("napscore_z.3", "napscore_z.5", "napscore_z.7", 
            "c_age", "c_gender", "c_dep.2", "c_dep.4", "c_dep.6", "c_ses", "c_nap1_z", 
            "id", "school")]
        
        
        # 3. Rename depression variables for reshape
        colnames(datw)[colnames(datw) == "c_dep.2"] <- "prev_dep.3"
        colnames(datw)[colnames(datw) == "c_dep.4"] <- "prev_dep.5"
        colnames(datw)[colnames(datw) == "c_dep.6"] <- "prev_dep.7"
        
        # 4. Reshape to long
        datL <- reshape(datw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z"), times = c(3, 5, 7), direction = "long")
        
        datL$prev_dep2 <- datL$prev_dep * datL$prev_dep
        datL <- datL[order(datL$id), ]
        
        # 5. Save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_dep2 + (1 | school/id), data = d)
    })
    
    # Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    # Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    FCSslwide_results.CI <- cbind(FCSslwide_results.CI, CI)
    
    # Store the estimates
    FCSslwide_results.est[, i] <- MI_est$estimates[2:8, 1]
    FCSslwide_results.sd[, i] <- MI_est$estimates[2:8, 2]
    FCSslwide_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    FCSslwide_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    FCSslwide_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    FCSslwide_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 1] + 
        MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    FCSslwide_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1])/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
}

rownames(FCSslwide_results.est) <- rows
colnames(FCSslwide_results.est) <- c(seq(1:length(data)))

rownames(FCSslwide_results.sd) <- rows
colnames(FCSslwide_results.sd) <- c(seq(1:length(data)))

rownames(FCSslwide_results.RE) <- c("level 3", "level 2", "level 1")
colnames(FCSslwide_results.RE) <- c(seq(1:length(data)))

rownames(FCSslwide_results.ICC) <- c("level 3", "level 2")
colnames(FCSslwide_results.ICC) <- c(seq(1:length(data)))

rownames(FCSslwide_results.CI) <- rows
colnames(FCSslwide_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(FCSslwide_results.est, "FCSslwide_results.est.xlsx")
write.xlsx(FCSslwide_results.sd, "FCSslwide_results.sd.xlsx")
write.xlsx(FCSslwide_results.RE, "FCSslwide_results.RE.xlsx")
write.xlsx(FCSslwide_results.ICC, "FCSslwide_results.ICC.xlsx")
write.xlsx(FCSslwide_results.CI, "FCSslwide_results.CI.xlsx")

###############################################
###        FCS-1L-DI-wide-Passive           ###
###############################################

set.seed(23357819) #MAR-CATS, 40 school clsuters
# set.seed(54546178) #MAR-CATS, 10 school clusters
# set.seed(237819) #MAR-inflated, 40 school clusters
# set.seed(58566178)  #MAR-inflated, 10 school clusters

## Generate empty matrices to save the results
FCSslwidepassive_results.est <- matrix(NA, nrow = 7, ncol = length(data))
FCSslwidepassive_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
FCSslwidepassive_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
FCSslwidepassive_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
FCSslwidepassive_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Create school dummy indicators
    simdataw$school <- as.factor(simdataw$school)
    
    ## Remove unwanted variables
    simdataw1 <- simdataw[, !names(simdataw) %in% c("c_id", "child", "napscore_z.2", 
        "napscore_z.4", "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", 
        "c_dep.5", "c_dep.7")]
    
    ## Generate the squared terms
    simdataw1$c_depsq.2 <- simdataw1$c_dep.2 * simdataw1$c_dep.2
    simdataw1$c_depsq.4 <- simdataw1$c_dep.4 * simdataw1$c_dep.4
    simdataw1$c_depsq.6 <- simdataw1$c_dep.6 * simdataw1$c_dep.6
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    
    # create the predictor matrix 
    pred <- make.predictorMatrix(simdataw1)
    pred[!rownames(pred) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_ses", "c_depsq.2", 
        "c_depsq.4", "c_depsq.6"), ] <- 0
    pred["c_dep.2", "c_depsq.2"] <- 0
    pred["c_dep.4", "c_depsq.4"] <- 0
    pred["c_dep.6", "c_depsq.6"] <- 0
    diag(pred) <- 0
    
    ## Sepcify the imputation method
    meth <- make.method(simdataw1)
    meth[substr(row.names(pred), 1, 6) %in% c("c_dep.")] <- "norm"
    meth[substr(row.names(pred), 1, 5) %in% c("c_ses")] <- "norm"
    meth["c_depsq.2"] <- "~I(c_dep.2*c_dep.2)"
    meth["c_depsq.4"] <- "~I(c_dep.4*c_dep.4)"
    meth["c_depsq.6"] <- "~I(c_dep.6*c_dep.6)"
    
    ## Perform imputations
    imp2 <- mice(data = simdataw1, m = M, predictorMatrix = pred, method = meth, 
        maxit = 10)
    
    ## Analysis
    
    mylist <- list()
    for (m in 1:M) {
        
        # 1. Extract the mth imputed dataset
        datw <- complete(imp2, m)
        
        # 2. Remove unwanted variables
        datw <- datw[, names(datw) %in% c("napscore_z.3", "napscore_z.5", "napscore_z.7", 
            "c_age", "c_gender", "c_dep.2", "c_dep.4", "c_dep.6", "c_depsq.2", "c_depsq.4", 
            "c_depsq.6", "c_ses", "c_nap1_z", "id", "school")]
        
        # 3. Rename depression variables for reshape
        colnames(datw)[colnames(datw) == "c_dep.2"] <- "prev_dep.3"
        colnames(datw)[colnames(datw) == "c_dep.4"] <- "prev_dep.5"
        colnames(datw)[colnames(datw) == "c_dep.6"] <- "prev_dep.7"
        
        
        colnames(datw)[colnames(datw) == "c_depsq.2"] <- "prev_depsq.3"
        colnames(datw)[colnames(datw) == "c_depsq.4"] <- "prev_depsq.5"
        colnames(datw)[colnames(datw) == "c_depsq.6"] <- "prev_depsq.7"
        
        # 4. Reshape to long
        datL <- reshape(datw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7"), c("prev_depsq.3", 
                "prev_depsq.5", "prev_depsq.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z", "prev_depsq"), times = c(3, 5, 7), direction = "long")
        
        datL <- datL[order(datL$id), ]
        
        # 5. Save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_depsq + (1 | school/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    FCSslwidepassive_results.CI <- cbind(FCSslwidepassive_results.CI, CI)
    
    ## Store the estimates
    FCSslwidepassive_results.est[, i] <- MI_est$estimates[2:8, 1]
    FCSslwidepassive_results.sd[, i] <- MI_est$estimates[2:8, 2]
    FCSslwidepassive_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    FCSslwidepassive_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    FCSslwidepassive_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    FCSslwidepassive_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    FCSslwidepassive_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
        1])/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
    
}
rownames(FCSslwidepassive_results.est) <- rows
colnames(FCSslwidepassive_results.est) <- c(seq(1:length(data)))

rownames(FCSslwidepassive_results.sd) <- rows
colnames(FCSslwidepassive_results.sd) <- c(seq(1:length(data)))

rownames(FCSslwidepassive_results.RE) <- c("level 3", "level 2", "level 1")
colnames(FCSslwidepassive_results.RE) <- c(seq(1:length(data)))

rownames(FCSslwidepassive_results.ICC) <- c("level 3", "level 2")
colnames(FCSslwidepassive_results.ICC) <- c(seq(1:length(data)))

rownames(FCSslwidepassive_results.CI) <- rows
colnames(FCSslwidepassive_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(FCSslwidepassive_results.est, "FCSslwidepassive_results.est.xlsx")
write.xlsx(FCSslwidepassive_results.sd, "FCSslwidepassive_results.sd.xlsx")
write.xlsx(FCSslwidepassive_results.RE, "FCSslwidepassive_results.RE.xlsx")
write.xlsx(FCSslwidepassive_results.ICC, "FCSslwidepassive_results.ICC.xlsx")
write.xlsx(FCSslwidepassive_results.CI, "FCSslwidepassive_results.CI.xlsx")

###############################################
###         FCS-2L-wide-Passive             ###
###############################################

set.seed(79378020) #MAR-CATS, 40 school clsuters
# set.seed(26656781) #MAR-CATS, 10 school clusters
# set.seed(765237819) #MAR-inflated, 40 school clusters
# set.seed(26657781)  #MAR-inflated, 10 school clusters

## Generate empty matrices to save the results
FCSmlwidepassive_results.est <- matrix(NA, nrow = 7, ncol = length(data))
FCSmlwidepassive_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
FCSmlwidepassive_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
FCSmlwidepassive_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
FCSmlwidepassive_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$child), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Generate the squared terms
    simdataw$c_depsq.2 <- simdataw$c_dep.2 * simdataw$c_dep.2
    simdataw$c_depsq.4 <- simdataw$c_dep.4 * simdataw$c_dep.4
    simdataw$c_depsq.6 <- simdataw$c_dep.6 * simdataw$c_dep.6
    
    ## Remove unwanted variables
    simdataw1 <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7", 
        "child")]
    
    ## set the number of imputations and iterations
    M <- 20
    
    # create the predictor matrix
    pred <- make.predictorMatrix(simdataw1)
    pred[!rownames(pred) %in% c("school"), "school"] <- -2
    pred[!rownames(pred) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_ses", "c_depsq.2", 
        "c_depsq.4", "c_depsq.6"), ] <- 0
    pred["c_dep.2", "c_depsq.2"] <- 0
    pred["c_dep.4", "c_depsq.4"] <- 0
    pred["c_dep.6", "c_depsq.6"] <- 0
    pred[, "c_id"] <- 0
    diag(pred) <- 0
    
    ## Sepcify the imputation method
    meth <- mice::make.method(data = simdataw1)
    meth[substr(row.names(pred), 1, 5) %in% c("c_dep")] <- "2l.pan"
    meth[substr(row.names(pred), 1, 5) %in% c("c_ses")] <- "2l.pan"
    meth["c_depsq.2"] <- "~I(c_dep.2*c_dep.2)"
    meth["c_depsq.4"] <- "~I(c_dep.4*c_dep.4)"
    meth["c_depsq.6"] <- "~I(c_dep.6*c_dep.6)"
    
    # Perform the imputations
    imp4 <- mice(simdataw1, m = M, maxit = 10, predictorMatrix = pred, method = meth)
    
    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        # 1. Extract the ith imputed data set
        imputed <- complete(imp4, m)
        imputed <- imputed[, !names(imputed) %in% c("p_sdq.2", "p_sdq.4", "p_sdq.6")]
        
        # 2. Rename for reshape
        colnames(imputed)[colnames(imputed) == "c_dep.2"] <- "prev_dep.3"
        colnames(imputed)[colnames(imputed) == "c_dep.4"] <- "prev_dep.5"
        colnames(imputed)[colnames(imputed) == "c_dep.6"] <- "prev_dep.7"
        
        colnames(imputed)[colnames(imputed) == "c_depsq.2"] <- "prev_depsq.3"
        colnames(imputed)[colnames(imputed) == "c_depsq.4"] <- "prev_depsq.5"
        colnames(imputed)[colnames(imputed) == "c_depsq.6"] <- "prev_depsq.7"
        
        # 3. Reshape to long format
        datL <- reshape(imputed, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7"), c("prev_depsq.3", 
                "prev_depsq.5", "prev_depsq.7")), idvar = "c_id", v.names = c("prev_dep", 
            "napscore_z", "prev_depsq"), times = c(3, 5, 7), direction = "long")
        datL <- datL[order(datL$c_id), ]
        
        # 5. Save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_depsq + (1 | school/c_id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Store the estimates
    FCSmlwidepassive_results.est[, i] <- MI_est$estimates[2:8, 1]
    FCSmlwidepassive_results.sd[, i] <- MI_est$estimates[2:8, 2]
    FCSmlwidepassive_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    FCSmlwidepassive_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    FCSmlwidepassive_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    FCSmlwidepassive_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    FCSmlwidepassive_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
        1])/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
    # Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    FCSmlwidepassive_results.CI <- cbind(FCSmlwidepassive_results.CI, CI)
    
}

rownames(FCSmlwidepassive_results.est) <- rows
colnames(FCSmlwidepassive_results.est) <- c(seq(1:length(data)))

rownames(FCSmlwidepassive_results.sd) <- rows
colnames(FCSmlwidepassive_results.sd) <- c(seq(1:length(data)))

rownames(FCSmlwidepassive_results.RE) <- c("level 3", "level 2", "level 1")
colnames(FCSmlwidepassive_results.RE) <- c(seq(1:length(data)))

rownames(FCSmlwidepassive_results.ICC) <- c("level 3", "level 2")
colnames(FCSmlwidepassive_results.ICC) <- c(seq(1:length(data)))

rownames(FCSmlwidepassive_results.CI) <- rows
colnames(FCSmlwidepassive_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(FCSmlwidepassive_results.est, "FCSmlwidepassive_results.est.xlsx")
write.xlsx(FCSmlwidepassive_results.sd, "FCSmlwidepassive_results.sd.xlsx")
write.xlsx(FCSmlwidepassive_results.RE, "FCSmlwidepassive_results.RE.xlsx")
write.xlsx(FCSmlwidepassive_results.CI, "FCSmlwidepassive_results.CI.xlsx")
write.xlsx(FCSmlwidepassive_results.ICC, "FCSmlwidepassive_results.ICC.xlsx")

###############################################
###           JM-1L-DI-wide                 ###
###############################################

## Commands used for running the simulations in parallel
args <- commandArgs(trailingOnly = TRUE)  #These arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum <- as.numeric(args[1])

T1 <- list.files(pattern = "*.csv")

s2 <- seq(200, 1000, by = 200)
s1 <- s2 - 199

temp <- T1[s1[datnum]:s2[datnum]]

# Setting random seeds for the parallely run replications

#MAR-CATS, 40 school clsuters
set.seed(23901984)
seed=sample(1e8,5,replace=F)[datnum]

#MAR-CATS, 10 school clsuters
# set.seed(23922984)
# seed=sample.int(1e8,5,replace = F)[datnum] 

#MAR-inflated, 40 school clsuters
# set.seed(239984)
# seed=sample.int(1e8,5,replace = F)[datnum] 

#MAR-inflated, 10 school clsuters
# set.seed(23112984)
# seed=sample.int(1e8,5,replace = F)[datnum] 

set.seed(seed)

## Load the simulated data
data <- lapply(temp, read.csv)

## Generate empty matrices to save the results
MVNIslwide_results.est <- matrix(NA, nrow = 7, ncol = length(data))
MVNIslwide_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
MVNIslwide_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
MVNIslwide_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
MVNIslwide_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7")]
    
    ## Set number of imputations and burn in and thining intervals
    M <- 20
    nburn <- 1000
    NB <- 100
    
    ## Create a dataframe with variables to be imputed
    myvars <- names(simdataw) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_ses")
    dataimp <- simdataw[myvars]
    
    ## create school dummy indicators (edit accordingly for 10 school clusters)
    school_DI <- data.frame(model.matrix(simdataw$c_id ~ as.factor(simdataw$school) - 
        1, simdataw))
    names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x, y) paste(x, 
        seq(1, y), sep = "_"), "schoo_Ind", 40))
    school_DI <- school_DI[, 1:39]
    
    ## Create a dataframe with complete variables
    datacomp <- cbind(Intercept = rep(1, nrow(simdataw)), simdataw[, names(simdataw) %in% 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7", "c_age", "c_gender", "c_nap1_z", 
            "p_sdq.2", "p_sdq.4", "p_sdq.6")], school_DI)
    
    ## Perform imputations without the random effects (SL imputation)
    imp1 <- jomo1con(Y = dataimp, X = datacomp, nimp = M, nburn = nburn, nbetween = NB)
    
    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        datw <- imp1[imp1$Imputation == m, ]
        
        # 1. Attach the school variable and remove school indicator variables
        datw <- cbind(datw, simdataw$school)
        
        # 2. Remove unwanted variables
        datw <- datw[, names(datw) %in% c("napscore_z.3", "napscore_z.5", "napscore_z.7", 
            "c_age", "c_gender", "c_dep.2", "c_dep.4", "c_dep.6", "c_ses", "c_nap1_z", 
            "id", "simdataw$school")]
        names(datw)[names(datw) == "simdataw$school"] <- "school"
        
        # 3. Rename depression variables for reshape
        colnames(datw)[colnames(datw) == "c_dep.2"] <- "prev_dep.3"
        colnames(datw)[colnames(datw) == "c_dep.4"] <- "prev_dep.5"
        colnames(datw)[colnames(datw) == "c_dep.6"] <- "prev_dep.7"
        
        # 4. Reshape to long and generate the squared term
        datL <- reshape(datw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z"), times = c(3, 5, 7), direction = "long")
        
        datL$prev_dep2 <- datL$prev_dep * datL$prev_dep
        
        datL <- datL[order(datL$school, datL$id), ]
        
        # 5. Save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_dep2 + (1 | school/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    MVNIslwide_results.CI <- cbind(MVNIslwide_results.CI, CI)
    
    
    ## Store the estimates
    MVNIslwide_results.est[, i] <- MI_est$estimates[2:8, 1]
    MVNIslwide_results.sd[, i] <- MI_est$estimates[2:8, 2]
    MVNIslwide_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    MVNIslwide_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    MVNIslwide_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    MVNIslwide_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 1] + 
        MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    MVNIslwide_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1])/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
}

rows <- c("prev_dep", "time", "c_age", "c_gender", "c_nap1_z", "c_ses", "interaction")

rownames(MVNIslwide_results.est) <- rows
colnames(MVNIslwide_results.est) <- c(seq(1:length(data)))

rownames(MVNIslwide_results.sd) <- rows
colnames(MVNIslwide_results.sd) <- c(seq(1:length(data)))

rownames(MVNIslwide_results.RE) <- c("level 3", "level 2", "level 1")
colnames(MVNIslwide_results.RE) <- c(seq(1:length(data)))

rownames(MVNIslwide_results.ICC) <- c("level 3", "level 2")
colnames(MVNIslwide_results.ICC) <- c(seq(1:length(data)))

rownames(MVNIslwide_results.CI) <- rows
colnames(MVNIslwide_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(MVNIslwide_results.est, paste0("MVNIslwide_results.est", datnum, ".xlsx"))
write.xlsx(MVNIslwide_results.sd, paste0("MVNIslwide_results.sd", datnum, ".xlsx"))
write.xlsx(MVNIslwide_results.RE, paste0("MVNIslwide_results.RE", datnum, ".xlsx"))
write.xlsx(MVNIslwide_results.CI, paste0("MVNIslwide_results.CI", datnum, ".xlsx"))
write.xlsx(MVNIslwide_results.ICC, paste0("MVNIslwide_results.ICC", datnum, ".xlsx"))

###############################################
###          JM-1L-DI-wide-JAV              ###
###############################################

## Commands used for running the simulations in parallel
args <- commandArgs(trailingOnly = TRUE)  #These arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum <- as.numeric(args[1])

T1 <- list.files(pattern = "*.csv")

s2 <- seq(20, 1000, by = 20)
s1 <- s2 - 19

temp <- T1[s1[datnum]:s2[datnum]]

# Setting random seeds for the parallely run replications

#MAR-CATS, 40 school clsuters
set.seed(3444627)
seed=sample(1e8,50,replace=F)[datnum]

#MAR-CATS, 10 school clsuters
# set.seed(26915984)
# seed=sample.int(1e8,50,replace = F)[datnum] 

#MAR-inflated, 40 school clsuters
# set.seed(3554617)
# seed=sample.int(1e8,50,replace = F)[datnum] 

#MAR-inflated, 10 school clsuters
# set.seed(26315384)
# seed=sample.int(1e8,50,replace = F)[datnum] 

set.seed(seed)

## Load the simulated data
data <- lapply(temp, read.csv)

## Generate empty matrices to save the results
MVNIslwideJAV_results.est <- matrix(NA, nrow = 7, ncol = length(data))
MVNIslwideJAV_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
MVNIslwideJAV_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
MVNIslwideJAV_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
MVNIslwideJAV_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7")]
    
    ## Set number of imputations
    M <- 20
    nburn <- 1000
    NB <- 100
    
    ## Generate the quadratic terms
    simdataw$c_depsq.2 <- simdataw$c_dep.2 * simdataw$c_dep.2
    simdataw$c_depsq.4 <- simdataw$c_dep.4 * simdataw$c_dep.4
    simdataw$c_depsq.6 <- simdataw$c_dep.6 * simdataw$c_dep.6
    
    ## Create a dataframe with variables to be imputed
    myvars <- names(simdataw) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_depsq.2", 
        "c_depsq.4", "c_depsq.6", "c_ses")
    dataimp <- simdataw[myvars]
    
    ## Create school dummy indicators(edit accordingly for 10 school clusters)
    school_DI <- data.frame(model.matrix(simdataw$c_id ~ as.factor(simdataw$school) - 
        1, simdataw))
    names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x, y) paste(x, 
        seq(1, y), sep = "_"), "schoo_Ind", 40))
    school_DI <- school_DI[, 1:39]
    
    ## Create a dataframe with complete variables
    datacomp <- cbind(Intercept = rep(1, nrow(simdataw)), simdataw[, names(simdataw) %in% 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7", "c_age", "c_gender", "c_nap1_z", 
            "p_sdq.2", "p_sdq.4", "p_sdq.6")], school_DI)
    
    ## Perform imputations without the random effects (SL imputation)
    imp1 <- jomo1con(Y = dataimp, X = datacomp, nimp = M, nburn = nburn, nbetween = NB)
    
    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        # 1. Extract the mth imputed dataset
        datw <- imp1[imp1$Imputation == m, ]
        
        # 2. Attach the school variable and remove school indicator variables
        datw <- cbind(datw, simdataw$school)
        
        # 3. Remove unwanted variables
        datw <- datw[, names(datw) %in% c("napscore_z.3", "napscore_z.5", "napscore_z.7", 
            "c_age", "c_gender", "c_dep.2", "c_dep.4", "c_dep.6", "c_depsq.2", "c_depsq.4", 
            "c_depsq.6", "c_ses", "c_nap1_z", "id", "simdataw$school")]
        
        colnames(datw)[colnames(datw) == "simdataw$school"] <- "school"
        
        # 4. Rename depression variables for reshape
        colnames(datw)[colnames(datw) == "c_dep.2"] <- "prev_dep.3"
        colnames(datw)[colnames(datw) == "c_dep.4"] <- "prev_dep.5"
        colnames(datw)[colnames(datw) == "c_dep.6"] <- "prev_dep.7"
        
        colnames(datw)[colnames(datw) == "c_depsq.2"] <- "prev_depsq.3"
        colnames(datw)[colnames(datw) == "c_depsq.4"] <- "prev_depsq.5"
        colnames(datw)[colnames(datw) == "c_depsq.6"] <- "prev_depsq.7"
        
        # 5. Reshape to long
        datL <- reshape(datw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7"), c("prev_depsq.3", 
                "prev_depsq.5", "prev_depsq.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z", "prev_depsq"), times = c(3, 5, 7), direction = "long")
        
        datL <- datL[order(datL$school, datL$id), ]
        
        # 6. save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_depsq + (1 | school/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    MVNIslwideJAV_results.CI <- cbind(MVNIslwideJAV_results.CI, CI)
    
    ## Store the estimates
    MVNIslwideJAV_results.est[, i] <- MI_est$estimates[2:8, 1]
    MVNIslwideJAV_results.sd[, i] <- MI_est$estimates[2:8, 2]
    MVNIslwideJAV_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    MVNIslwideJAV_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    MVNIslwideJAV_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    MVNIslwideJAV_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    MVNIslwideJAV_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
        1])/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
}

rownames(MVNIslwideJAV_results.est) <- rows
colnames(MVNIslwideJAV_results.est) <- c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.sd) <- rows
colnames(MVNIslwideJAV_results.sd) <- c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.RE) <- c("level 3", "level 2", "level 1")
colnames(MVNIslwideJAV_results.RE) <- c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.ICC) <- c("level 3", "level 2")
colnames(MVNIslwideJAV_results.ICC) <- c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.CI) <- rows
colnames(MVNIslwideJAV_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(MVNIslwideJAV_results.est, paste0("MVNIslwideJAV_results.est", datnum, ".xlsx"))
write.xlsx(MVNIslwideJAV_results.sd, paste0("MVNIslwideJAV_results.sd", datnum, ".xlsx"))
write.xlsx(MVNIslwideJAV_results.RE, paste0("MVNIslwideJAV_results.RE", datnum, ".xlsx"))
write.xlsx(MVNIslwideJAV_results.CI, paste0("MVNIslwideJAV_results.CI", datnum, ".xlsx"))
write.xlsx(MVNIslwideJAV_results.ICC, paste0("MVNIslwideJAV_results.ICC", datnum, ".xlsx"))

###############################################
###             JM-2L-wide-JAV              ###
###############################################

## Commands used for running the simulations in parallel
args <- commandArgs(trailingOnly = TRUE)  #These arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum<-as.numeric(args[1])  

T1 <- list.files(pattern = "*.csv")

s2 <- seq(200, 1000, by = 200)
s1 <- s2 - 199

temp <- T1[s1[datnum]:s2[datnum]]

# setting random seeds for the parallely run replications

#MAR-CATS, 40 school clsuters
set.seed(367201)
seed=sample(1e8,50,replace=F)[datnum]

#MAR-CATS, 10 school clsuters
# set.seed(23422914)
# seed=sample.int(1e8,50,replace = F)[datnum] 

#MAR-inflated, 40 school clsuters
# set.seed(12434179)
# seed=sample.int(1e8,50,replace = F)[datnum] 

#MAR-inflated, 10 school clsuters
# set.seed(23411914)
# seed=sample.int(1e8,50,replace = F)[datnum] 

set.seed(seed)

## Load the simulated data
data <- lapply(temp, read.csv)

## Generate empty matrices to save the results
MVNImlwideJAV_results.est <- matrix(NA, nrow = 7, ncol = length(data))
MVNImlwideJAV_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
MVNImlwideJAV_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
MVNImlwideJAV_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
MVNImlwideJAV_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]

    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7")]
    
    simdataw$c_depsq.2 <- simdataw$c_dep.2 * simdataw$c_dep.2
    simdataw$c_depsq.4 <- simdataw$c_dep.4 * simdataw$c_dep.4
    simdataw$c_depsq.6 <- simdataw$c_dep.6 * simdataw$c_dep.6
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    nburn <- 1000
    NB <- 100
    
    
    ## Create a dataframe with variables to be imputed
    myvars <- names(simdataw) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_depsq.2", 
        "c_depsq.4", "c_depsq.6", "c_ses")
    dataimp <- simdataw[myvars]
    
    ## Create a dataframe with complete variables
    datacomp <- cbind(Intercept = rep(1, nrow(simdataw)), simdataw[, names(simdataw) %in% 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7", "c_age", "c_gender", "c_nap1_z", 
            "p_sdq.2", "p_sdq.4", "p_sdq.6")])
    
    # Create a data frame with column of 1's for random intercept
    datcompRE <- cbind(Intercept = rep(1, nrow(simdataw)))
    
    ## Perform imputations without the random effects
    imp5 <- jomo(Y = dataimp, X = datacomp, Z = datcompRE, clus = simdataw$school, 
        nimp = M, nburn = nburn, nbetween = NB)
    
    mylist <- list()
    
    ##Analysis
    for (m in 1:M) {
        
        # 1. Extract the ith imputed data set
        imputed <- imp5[imp5$Imputation == m, ]
        imputed <- imputed[, !names(imputed) %in% c("p_sdq.2", "p_sdq.4", "p_sdq.6", 
            "Intercept", "Imputation")]
        
        # 2. Rename depression for reshape
        colnames(imputed)[colnames(imputed) == "c_dep.2"] <- "prev_dep.3"
        colnames(imputed)[colnames(imputed) == "c_dep.4"] <- "prev_dep.5"
        colnames(imputed)[colnames(imputed) == "c_dep.6"] <- "prev_dep.7"
        
        colnames(imputed)[colnames(imputed) == "c_depsq.2"] <- "prev_depsq.3"
        colnames(imputed)[colnames(imputed) == "c_depsq.4"] <- "prev_depsq.5"
        colnames(imputed)[colnames(imputed) == "c_depsq.6"] <- "prev_depsq.7"
        
        # 3. Reshape to long format
        datL <- reshape(imputed, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7"), c("prev_depsq.3", 
                "prev_depsq.5", "prev_depsq.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z", "prev_depsq"), times = c(3, 5, 7), direction = "long")
        
        datL <- datL[order(datL$clus, datL$id), ]
        
        # 4. save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_depsq + (1 | clus/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Store the estimates
    MVNImlwideJAV_results.est[, i] <- MI_est$estimates[2:8, 1]
    MVNImlwideJAV_results.sd[, i] <- MI_est$estimates[2:8, 2]
    MVNImlwideJAV_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    MVNImlwideJAV_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    MVNImlwideJAV_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    MVNImlwideJAV_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    MVNImlwideJAV_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
        1])/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    MVNImlwideJAV_results.CI <- cbind(MVNImlwideJAV_results.CI, CI)
}

rownames(MVNImlwideJAV_results.est) <- rows
colnames(MVNImlwideJAV_results.est) <- c(seq(1:length(data)))

rownames(MVNImlwideJAV_results.sd) <- rows
colnames(MVNImlwideJAV_results.sd) <- c(seq(1:length(data)))

rownames(MVNImlwideJAV_results.RE) <- c("level 3", "level 2", "level 1")
colnames(MVNImlwideJAV_results.RE) <- c(seq(1:length(data)))

rownames(MVNImlwideJAV_results.ICC) <- c("level 3", "level 2")
colnames(MVNImlwideJAV_results.ICC) <- c(seq(1:length(data)))

rownames(MVNImlwideJAV_results.CI) <- rows
colnames(MVNImlwideJAV_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(MVNImlwideJAV_results.est, paste0("MVNImlwideJAV_results.est", datnum, ".xlsx"))
write.xlsx(MVNImlwideJAV_results.sd, paste0("MVNImlwideJAV_results.sd", datnum, ".xlsx"))
write.xlsx(MVNImlwideJAV_results.RE, paste0("MVNImlwideJAV_results.RE", datnum, ".xlsx"))
write.xlsx(MVNImlwideJAV_results.CI, paste0("MVNImlwideJAV_results.CI", datnum, ".xlsx"))
write.xlsx(MVNImlwideJAV_results.ICC, paste0("MVNImlwideJAV_results.ICC", datnum, ".xlsx"))

###############################################
###              SMC-JM-2L-DI               ###
###############################################

## Commands used for running the simulations in parallel
args <- commandArgs(trailingOnly = TRUE)  #These arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum<-as.numeric(args[1])   

T1 <- list.files(pattern = "*.csv")

s2 <- seq(20, 1000, by = 20)
s1 <- s2 - 19

temp <- T1[s1[datnum]:s2[datnum]]

#setting random seeds for the parallely run replications

#MAR-CATS, 40 school clsuters
set.seed(367201)
seed=sample(1e8,50,replace=F)[datnum]

#MAR-CATS, 10 school clsuters
# set.seed(23422914)
# seed=sample.int(1e8,50,replace = F)[datnum] 

#MAR-inflated, 40 school clsuters
# set.seed(12434179)
# seed=sample.int(1e8,50,replace = F)[datnum] 

#MAR-inflated, 10 school clsuters
# set.seed(23411914)
# seed=sample.int(1e8,50,replace = F)[datnum] 

set.seed(seed)

## Load the simulated data
data <- lapply(temp, read.csv)

## Generate empty matrices to save the results
SMCJMDI_results.est <- matrix(NA, nrow = 7, ncol = length(data))
SMCJMDI_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
SMCJMDI_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
SMCJMDI_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
SMCJMDI_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$child), ]
    
    ## Rearrange the data set for imputations
    # 1. SDQ variable
    simdataL <- slide(simdataL, Var = "p_sdq", GroupVar = "c_id", slideBy = -1)
    simdataL <- simdataL[, !names(simdataL) %in% c("p_sdq")]
    colnames(simdataL)[colnames(simdataL) == "p_sdq-1"] <- "prev_sdq"
    
    # 2. Depression values at waves 2, 4 and 6
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "prev_sdq", "c_dep"), 
        idvar = "c_id", timevar = "wave", direction = "wide")
    
    colnames(simdataw)[colnames(simdataw) == "c_dep.2"] <- "prev_dep.3"
    colnames(simdataw)[colnames(simdataw) == "c_dep.4"] <- "prev_dep.5"
    colnames(simdataw)[colnames(simdataw) == "c_dep.6"] <- "prev_dep.7"
    
    # 3. Remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "prev_sdq.2", 
        "napscore_z.4", "prev_sdq.4", "napscore_z.6", "prev_sdq.6", "c_dep.3", "c_dep.5", 
        "c_dep.7")]
    
    ## Reshape back to long
    simdataL <- reshape(simdataw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7"), c("prev_sdq.3", "prev_sdq.5", 
            "prev_sdq.7")), idvar = "c_id", v.names = c("prev_dep", "napscore_z", 
        "prev_sdq"), times = c(3, 5, 7), direction = "long")
    
    simdataL$c_id <- as.numeric(simdataL$c_id)
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Recode sex as a factor
    simdataL$c_gender <- as.factor(simdataL$c_gender)
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    nburn <- 500
    NB <- 10
    
    ## Create school dummy indicators(edit accordingly for 10 school clusters)
    school_DI <- data.frame(model.matrix(simdataL$c_id ~ as.factor(simdataL$school) - 
        1, simdataL))
    names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x, y) paste(x, 
        seq(1, y), sep = "_"), "schoo_Ind", 40))
    school_DI <- school_DI[, 1:39]
    
    ## Define a dataframe with all variables required
    data_jomo <- cbind(simdataL[!names(simdataL) %in% c("child", "school")], school_DI)
    
    # Specify the levels of each variable in the imputation model
    mylevel <- c(2, 2, 2, 2, 2, 1, 1, 1, 1, rep(2, times = 39))
    
    ## Formula of the substantive lmer model (change accordingly for 10 school
    # clusters)
    formula <- as.formula(napscore_z ~ prev_dep + time + I(prev_dep^2) + prev_sdq + 
        c_age + c_gender + c_ses + c_nap1_z + schoo_Ind_1 + schoo_Ind_2 + schoo_Ind_3 + 
        schoo_Ind_4 + schoo_Ind_5 + schoo_Ind_6 + schoo_Ind_7 + schoo_Ind_8 + schoo_Ind_9 + 
        schoo_Ind_10 + schoo_Ind_11 + schoo_Ind_12 + schoo_Ind_13 + schoo_Ind_14 + 
        schoo_Ind_15 + schoo_Ind_16 + schoo_Ind_17 + schoo_Ind_18 + schoo_Ind_19 + 
        schoo_Ind_20 + schoo_Ind_21 + schoo_Ind_22 + schoo_Ind_23 + schoo_Ind_24 + 
        schoo_Ind_25 + schoo_Ind_26 + schoo_Ind_27 + schoo_Ind_28 + schoo_Ind_29 + 
        schoo_Ind_30 + schoo_Ind_31 + schoo_Ind_32 + schoo_Ind_33 + schoo_Ind_34 + 
        schoo_Ind_35 + schoo_Ind_36 + schoo_Ind_37 + schoo_Ind_38 + schoo_Ind_39 + 
        (1 | c_id))
    
    ## Running the imputations
    imp5 <- jomo.lmer(formula, data_jomo, level = mylevel, nimp = M, nburn = nburn, 
        nbetween = NB)
    
    ##Analysis
    mylist <- list()
   
    for (m in 1:M) {
        
        # 1.Extract the ith imputed data set and attach the school variable
        imputed <- imp5[imp5$Imputation == m, ]
        imputed <- imputed[, names(imputed) %in% c("prev_dep", "c_age", "c_gender", 
            "c_ses", "c_nap1_z", "time", "napscore_z", "clus")]
        imputed$school <- simdataL$school

        imputed$prev_dep2 <- imputed$prev_dep * imputed$prev_dep
        # 2. save the dataset in a list
        mylist[[m]] <- imputed
    }
    
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + c_age + c_gender + c_nap1_z + c_ses + 
            prev_dep2 + (1 | school/clus), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    SMCJMDI_results.CI <- cbind(SMCJMDI_results.CI, CI)
    
    # Store the estimates
    SMCJMDI_results.est[, i] <- MI_est$estimates[2:8, 1]
    SMCJMDI_results.sd[, i] <- MI_est$estimates[2:8, 2]
    SMCJMDI_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    SMCJMDI_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    SMCJMDI_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    SMCJMDI_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
        1] + MI_est$var.comp[3, 1])
    SMCJMDI_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1])/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
}

rownames(SMCJMDI_results.est) <- rows
colnames(SMCJMDI_results.est) <- c(seq(1:length(data)))

rownames(SMCJMDI_results.sd) <- rows
colnames(SMCJMDI_results.sd) <- c(seq(1:length(data)))

rownames(SMCJMDI_results.RE) <- c("level 3", "level 2", "level 1")
colnames(SMCJMDI_results.RE) <- c(seq(1:length(data)))

rownames(SMCJMDI_results.ICC) <- c("level 3", "level 2")
colnames(SMCJMDI_results.ICC) <- c(seq(1:length(data)))

rownames(SMCJMDI_results.CI) <- rows
colnames(SMCJMDI_results.CI) <- c(rep(1:length(data), each = 2))

## Save the results
write.xlsx(SMCJMDI_results.est, paste0("SMCJMDI_results.est", datnum, ".xlsx"))
write.xlsx(SMCJMDI_results.sd, paste0("SMCJMDI_results.sd", datnum, ".xlsx"))
write.xlsx(SMCJMDI_results.RE, paste0("SMCJMDI_results.RE", datnum, ".xlsx"))
write.xlsx(SMCJMDI_results.CI, paste0("SMCJMDI_results.CI", datnum, ".xlsx"))
write.xlsx(SMCJMDI_results.ICC, paste0("SMCJMDI_results.ICC", datnum, ".xlsx"))

###############################################
###               SMC-SM-2L-DI              ###
###############################################

args <- commandArgs(trailingOnly = TRUE)  #These arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum <- as.numeric(args[1])

T1 <- list.files(pattern = "*.csv")

s2 <- seq(100, 1000, by = 100)
s1 <- s2 - 99

temp <- T1[s1[datnum]:s2[datnum]]

# setting random seeds for the parallely run replications

#MAR-CATS, 40 school clsuters
set.seed(790138093)
seed=sample(1e8,10,replace=F)[datnum]

#MAR-CATS, 10 school clsuters
# set.seed(23422914)
# seed=sample.int(1e8,10,replace = F)[datnum] 

#MAR-inflated, 40 school clsuters
# set.seed(434156412)
# seed=sample.int(1e8,10,replace = F)[datnum] 

#MAR-inflated, 10 school clsuters
# set.seed(1381112)
# seed=sample.int(1e8,10,replace = F)[datnum] 

set.seed(seed)

## Load the simulated data
data <- lapply(temp, read.csv)

## Generate empty matrices to save the results
SMCSMDI_results.est <- matrix(NA, nrow = 7, ncol = length(data))
SMCSMDI_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
SMCSMDI_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
SMCSMDI_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
SMCSMDI_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$child), ]
    
    ## Create the previous wave deperssion symptoms
    simdataL <- slide(simdataL, Var = "c_dep", GroupVar = "c_id", slideBy = -1)
    colnames(simdataL)[colnames(simdataL) == "c_dep-1"] <- "prev_dep"
    simdataL$c_dep <- NULL
    
    ## create previous wave SDQ variable
    simdataL <- slide(simdataL, Var = "p_sdq", GroupVar = "c_id", slideBy = -1)
    colnames(simdataL)[colnames(simdataL) == "p_sdq-1"] <- "prev_sdq"
    simdataL$p_sdq <- NULL
    
    ## Remove unwanted waves
    simdataL <- subset(simdataL, wave != 2 & wave != 4 & wave != 6)
    simdataL$c_gender <- ifelse(simdataL$c_gender == "male", 1, 0)
    
    ## Create DI for schools
    simdataL$school <- simdataL$school - 1
    simdataL$c_id <- as.integer(simdataL$c_id)
    
    ## Model specification
    iter <- 2000
    burnin <- 1000  ##burn in for 1000 iteration and one imputed dataset 
    Nimp <- 20  #is saved every 100th iteration=20imp
    
    ## Substantive model formula
    Y_formula <- napscore_z ~ prev_dep + wave + I(prev_dep^2) + c_age + c_gender + 
        c_ses + c_nap1_z + school + (1 | c_id)
    
    ## Model for the outcome/Y variable
    dep <- list(model = "mlreg", formula = Y_formula)
    
    ## Covariate models 
    
    # 1.prev_dep
    X1_formula <- prev_dep ~ wave + c_age + c_ses + c_gender + c_nap1_z + prev_sdq + 
        school + (1 | c_id)
    
    ind_x <- list(model = "mlreg", formula = X1_formula, sampling_level = "c_id")
    
    # 2. SES
    X2_formula <- c_ses ~ c_age + c_gender + c_nap1_z + school
    
    ind_x2 <- list(model = "linreg", formula = X2_formula, variable_level = "c_id")
    
    ## Ordering of the conditional covariate models
    ind <- list(c_ses = ind_x2, prev_dep = ind_x)
    
    ## Estimate model
    mod1 <- mdmb::frm_fb(simdataL, dep, ind, burnin = burnin, iter = iter, Nimp = 20, 
        aggregation = TRUE)
    
    # Convert output into list of imputed datasets
    datlist <- mdmb::frm2datlist(mod1)
    
    ## Analysis
    mylist <- list()
    
    for (j in 1:Nimp) {
        
        # 1. Extract the jth imputed dataset
        datL <- datlist[[j]]
        
        #2. Remove unwanted variables
        datL <- datL[, names(datL) %in% c("c_id", "wave", "napscore_z", "c_age", 
            "c_gender", "c_nap1_z", "prev_dep", "prev_sdq", "c_ses", "school")]
        
        # 3. Generate the squared term
        datL$prevdep2 <- datL$prev_dep * datL$prev_dep
        
        # 4. Save the dataset in a list
        mylist[[j]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + wave + c_age + as.factor(c_gender) + c_nap1_z + 
            c_ses + prevdep2 + (1 | school/c_id), data = d)
    })
    
    MI_est <- testEstimates(mods, var.comp = TRUE, df.com = NULL)
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    SMCSMDI_results.CI <- cbind(SMCSMDI_results.CI, CI)
    
    
    ## Store the estimates
    SMCSMDI_results.est[, i] <- MI_est$estimates[2:8, 1]
    SMCSMDI_results.sd[, i] <- MI_est$estimates[2:8, 2]
    SMCSMDI_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    SMCSMDI_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    SMCSMDI_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    
    SMCSMDI_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 1] + MI_est$var.comp[1, 
        1] + MI_est$var.comp[3, 1])
    SMCSMDI_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1])/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
}

rownames(SMCSMDI_results.est) <- rows
colnames(SMCSMDI_results.est) <- c(seq(1:length(temp)))

rownames(SMCSMDI_results.sd) <- rows
colnames(SMCSMDI_results.sd) <- c(seq(1:length(temp)))

rownames(SMCSMDI_results.RE) <- c("level 3", "level 2", "level 1")
colnames(SMCSMDI_results.RE) <- c(seq(1:length(temp)))

rownames(SMCSMDI_results.ICC) <- c("level 3", "level 2")
colnames(SMCSMDI_results.ICC) <- c(seq(1:length(temp)))

rownames(SMCSMDI_results.CI) <- rows
colnames(SMCSMDI_results.CI) <- c(rep(1:length(temp), each = 2))

## Save the results
write.xlsx(SMCSMDI_results.est, paste0("SMCSMDI_results.est", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.sd, paste0("SMCSMDI_results.sd", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.RE, paste0("SMCSMDI_results.RE", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.CI, paste0("SMCSMDI_results.CI", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.ICC, paste0("SMCSMDI_results.ICC", datnum, ".xlsx"))
