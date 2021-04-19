#######################################################################################
#   Filename    :	MI_methods_analysis1.R 
#                   
#                   MI methods considered under the analysis model involving an interaction 
#                   between the time-varying exposure and time 
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
# 2.SMC-JM-2L-DI  (50 sets of 1000/50 repetitions)
# 3.SMC-SM-2L-DI  (10 sets of 1000/10 repetitions)

#Individual seeds used for each simulation scenario and the paralley run simulations are given 
#below as comments under each MI approach 

#********************************************************************************
# clear the workspace
rm(list = ls())

# Change the working directory to the relevant folder where the data are stored
setwd("C:\\temp")

## load the required packages
require(lme4)  #for fitting the linear mixed model
require(jomo)  # for single-level JM, two-level JM and two-level SMC-JM
require(mitml)  #for pooling the results from MI
require(mice)  #for single and two-level FCS
require(xlsx)  #for writing/exporting the results
require(DataCombine)  #for generating the lagged wave variables
require(mdmb)  #for two-level SMC-SM


# Import the simulated datasets
data <- lapply(list.files(pattern = glob2rx("*.csv")), read.csv)

###############################################
###              FCS-1L-DI-wide             ###
###############################################

set.seed(45271) #MAR-CATS, 40 school clsuters
# set.seed(6774257) #MAR-CATS, 10 school clusters
# set.seed(81112684) #MAR-Inflated, 40 school clusters
# set.seed(9757010)  #MAR-Inflated, 10 school clusters

## Generate empty matrices to save the results
FCSslwide_results.est <- matrix(NA, nrow = 7, ncol = length(data))
FCSslwide_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
FCSslwide_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
FCSslwide_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
FCSslwide_results.CI <- c()

for (i in 1:length(data)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7")]
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    
    ## Create school dummy indicators
    simdataw$school <- as.factor(simdataw$school)
    
    ## Perform imputations
    imp2 <- mice(data = simdataw[, !names(simdataw) %in% c("c_id", "child")], m = M, 
        method = "norm", maxit = 10)
    
    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        # 1. Extract the mth imputed dataset
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
        datL <- datL[order(datL$id), ]
        
        # 5. Save the dataset in a list
        mylist[[m]] <- datL
        
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + prev_dep * time + c_age + c_gender + 
            c_nap1_z + c_ses + (1 | school/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    FCSslwide_results.CI <- cbind(FCSslwide_results.CI, CI)
    
    ## Store the estimates
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

rows <- c("prev_dep", "time", "c_age", "c_gender", "c_nap1_z", "c_ses", "interaction")

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

# Save the results
write.xlsx(FCSslwide_results.est, "FCSslwide_results.est.xlsx")
write.xlsx(FCSslwide_results.sd, "FCSslwide_results.sd.xlsx")
write.xlsx(FCSslwide_results.RE, "FCSslwide_results.RE.xlsx")
write.xlsx(FCSslwide_results.ICC, "FCSslwide_results.ICC.xlsx")
write.xlsx(FCSslwide_results.CI, "FCSslwide_results.CI.xlsx")

###############################################
###                JM-2L-wide               ###
###############################################

set.seed(67809245)#MAR-CATS, 40 school clsuters
# set.seed(87930271) #MAR-CATS, 10 school clusters
# set.seed(80012684) #MAR-inflated, 40 school clusters
# set.seed(4558010)  #MAR-inflated, 10 school clusters

## Generate empty matrices to save the results
MVNImlwide_results.est <- matrix(NA, nrow = 7, ncol = length(data))
MVNImlwide_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
MVNImlwide_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
MVNImlwide_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
MVNImlwide_results.CI <- c()

for (i in 1:length(temp)) {
    
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7")]
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    nburn <- 1000
    NB <- 100
    
    ## Create a dataframe with variables to be imputed
    myvars <- names(simdataw) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_ses")
    dataimp <- simdataw[myvars]
    
    ## Create a dataframe with complete variables
    datacomp <- cbind(Intercept = rep(1, nrow(simdataw)), simdataw[, names(simdataw) %in% 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7", "c_age", "c_gender", "c_nap1_z", 
            "p_sdq.2", "p_sdq.4", "p_sdq.6")])
    
    ## Create a data frame with column of 1's for random intercept
    datcompRE <- cbind(Intercept = rep(1, nrow(simdataw)))
    
    ## Perform imputations without the random effects
    imp5 <- jomo1rancon(Y = dataimp, X = datacomp, Z = datcompRE, clus = simdataw$school, 
        nimp = M, nburn = nburn, nbetween = NB)
    
    ## Analysis
    
    mylist <- list()
    for (m in 1:M) {
        
        # 1. Extract the ith imputed data set
        imputed <- imp5[imp5$Imputation == m, ]
        
        # 2. Rmove unwanted variables
        imputed <- imputed[, !names(imputed) %in% c("p_sdq.2", "p_sdq.4", "p_sdq.6", 
            "Intercept", "Imputation")]
        
        # 3. Rename depression for reshape
        colnames(imputed)[colnames(imputed) == "c_dep.2"] <- "prev_dep.3"
        colnames(imputed)[colnames(imputed) == "c_dep.4"] <- "prev_dep.5"
        colnames(imputed)[colnames(imputed) == "c_dep.6"] <- "prev_dep.7"
        
        # 4.Reshape to long format
        datL <- reshape(imputed, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z"), times = c(3, 5, 7), direction = "long")
        
        datL <- datL[order(datL$clus, datL$id), ]
        
        # 5. save the dataset in a list
        mylist[[m]] <- datL
    }
    
    # Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + prev_dep * time + c_age + c_gender + 
            c_nap1_z + c_ses + (1 | clus/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    # Store the estimates
    MVNImlwide_results.est[, i] <- MI_est$estimates[2:8, 1]
    MVNImlwide_results.sd[, i] <- MI_est$estimates[2:8, 2]
    MVNImlwide_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    MVNImlwide_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    MVNImlwide_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    MVNImlwide_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 1] + 
        MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    MVNImlwide_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1])/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    MVNImlwide_results.CI <- cbind(MVNImlwide_results.CI, CI)
}

rownames(MVNImlwide_results.est) <- rows
colnames(MVNImlwide_results.est) <- c(seq(1:length(temp)))

rownames(MVNImlwide_results.sd) <- rows
colnames(MVNImlwide_results.sd) <- c(seq(1:length(temp)))

rownames(MVNImlwide_results.RE) <- c("level 3", "level 2", "level 1")
colnames(MVNImlwide_results.RE) <- c(seq(1:length(temp)))

rownames(MVNImlwide_results.ICC) <- c("level 3", "level 2")
colnames(MVNImlwide_results.ICC) <- c(seq(1:length(temp)))

rownames(MVNImlwide_results.CI) <- rows
colnames(MVNImlwide_results.CI) <- c(rep(1:length(temp), each = 2))

# Save the results
write.xlsx(MVNImlwide_results.est, "MVNImlwide_results.est.xlsx")
write.xlsx(MVNImlwide_results.sd, "MVNImlwide_results.sd.xlsx")
write.xlsx(MVNImlwide_results.RE, "MVNImlwide_results.RE.xlsx")
write.xlsx(MVNImlwide_results.CI, "MVNImlwide_results.CI.xlsx")
write.xlsx(MVNImlwide_results.ICC, "MVNImlwide_results.ICC.xlsx")

###############################################
###                FCS-2L-wide              ###
###############################################

set.seed(7889134)#MAR-CATS, 40 school clsuters
# set.seed(9017829) #MAR-CATS, 10 school clusters
# set.seed(9278901) #MAR-inflated, 40 school clusters
# set.seed(87572801)  #MAR-inflated, 10 school clusters

## Generate empty matrices to save the results
FCSmlwide_results.est <- matrix(NA, nrow = 7, ncol = length(data))
FCSmlwide_results.sd <- matrix(NA, nrow = 7, ncol = length(data))
FCSmlwide_results.RE <- matrix(NA, nrow = 3, ncol = length(data))
FCSmlwide_results.ICC <- matrix(NA, nrow = 2, ncol = length(data))
FCSmlwide_results.CI <- c()

for (i in 1:length(data)) {
    simdataL <- data[[i]]
    simdataL <- simdataL[order(simdataL$school, simdataL$child), ]
    
    ## Reshape to wide format
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "p_sdq", "c_dep"), idvar = "c_id", 
        timevar = "wave", direction = "wide")
    
    ## Remove unwanted variables
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "napscore_z.4", 
        "napscore_z.6", "p_sdq.3", "p_sdq.5", "p_sdq.7", "c_dep.3", "c_dep.5", "c_dep.7", 
        "child")]
    
    ## Set the number of imputations and iterations
    M <- 20
    
    ## Create the predictor matrix: this serves to specify the predictors of each
    ## imputation model
    pred <- make.predictorMatrix(simdataw)
    
    # A value of 1 means that the column variable is used as a predictor for the
    # target block (in the rows) A value of 2 indicates that the column variable is a
    # predictor with FIXED AND RANDOM effects for row variable A value of -2
    # indicates the cluster variable Other cells are set to zero
    pred[!rownames(pred) %in% c("school"), "school"] <- -2
    pred[!rownames(pred) %in% c("c_dep.2", "c_dep.4", "c_dep.6", "c_ses"), ] <- 0
    pred[, "c_id"] <- 0
    diag(pred) <- 0
    
    ## Imputation method
    meth <- rep("", ncol(pred))
    meth[substr(row.names(pred), 1, 5) %in% c("c_dep")] <- "2l.pan"
    meth[substr(row.names(pred), 1, 5) %in% c("c_ses")] <- "2l.pan"
    
    ## Perform the imputations
    imp4 <- mice(simdataw, m = M, maxit = 10, predictorMatrix = pred, method = meth)
    
    ## analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        # 1. Extract the ith imputed data set
        imputed <- complete(imp4, m)
        
        # 2. Remove unwanted variables
        imputed <- imputed[, !names(imputed) %in% c("p_sdq.2", "p_sdq.4", "p_sdq.6")]
        
        # 3. Rename variables
        colnames(imputed)[colnames(imputed) == "c_dep.2"] <- "prev_dep.3"
        colnames(imputed)[colnames(imputed) == "c_dep.4"] <- "prev_dep.5"
        colnames(imputed)[colnames(imputed) == "c_dep.6"] <- "prev_dep.7"
        
        # 4. Reshape to long format
        datL <- reshape(imputed, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7")), idvar = "c_id", v.names = c("prev_dep", 
            "napscore_z"), times = c(3, 5, 7), direction = "long")
        
        datL <- datL[order(datL$school, datL$c_id), ]
        
        # 5. Save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + prev_dep * time + c_age + c_gender + 
            c_nap1_z + c_ses + (1 | school/c_id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Store the estimates
    FCSmlwide_results.est[, i] <- MI_est$estimates[2:8, 1]
    FCSmlwide_results.sd[, i] <- MI_est$estimates[2:8, 2]
    FCSmlwide_results.RE[1, i] <- sqrt(MI_est$var.comp[2, 1])
    FCSmlwide_results.RE[2, i] <- sqrt(MI_est$var.comp[1, 1])
    FCSmlwide_results.RE[3, i] <- sqrt(MI_est$var.comp[3, 1])
    
    
    FCSmlwide_results.ICC[1, i] <- MI_est$var.comp[2, 1]/(MI_est$var.comp[2, 1] + 
        MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    FCSmlwide_results.ICC[2, i] <- (MI_est$var.comp[2, 1] + MI_est$var.comp[1, 1])/(MI_est$var.comp[2, 
        1] + MI_est$var.comp[1, 1] + MI_est$var.comp[3, 1])
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    FCSmlwide_results.CI <- cbind(FCSmlwide_results.CI, CI)
    
}

rownames(FCSmlwide_results.est) <- rows
colnames(FCSmlwide_results.est) <- c(seq(1:length(data)))

rownames(FCSmlwide_results.sd) <- rows
colnames(FCSmlwide_results.sd) <- c(seq(1:length(data)))

rownames(FCSmlwide_results.RE) <- c("level 3", "level 2", "level 1")
colnames(FCSmlwide_results.RE) <- c(seq(1:length(data)))

rownames(FCSmlwide_results.ICC) <- c("level 3", "level 2")
colnames(FCSmlwide_results.ICC) <- c(seq(1:length(data)))

rownames(FCSmlwide_results.CI) <- rows
colnames(FCSmlwide_results.CI) <- c(rep(1:length(data), each = 2))

# Save the results
write.xlsx(FCSmlwide_results.est, "FCSmlwide_results.est.xlsx")
write.xlsx(FCSmlwide_results.sd, "FCSmlwide_results.sd.xlsx")
write.xlsx(FCSmlwide_results.RE, "FCSmlwide_results.RE.xlsx")
write.xlsx(FCSmlwide_results.CI, "FCSmlwide_results.CI.xlsx")
write.xlsx(FCSmlwide_results.ICC, "FCSmlwide_results.ICC.xlsx")

###############################################
###              JM-1L-DI-wide              ###
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

# MAR-CATS, 40 school clsuters
set.seed(234562)
seed <- sample(1e+08, 5, replace = F)[datnum]

# MAR-CATS, 10 school clsuters
# set.seed(6718888)
# seed=sample.int(1e8,5,replace = F)[datnum] 

# MAR-inflated, 40 school clsuters
# set.seed(79010093)
# seed=sample.int(1e8,5,replace = F)[datnum] 

# MAR-inflated, 10 school clsuters
# set.seed(863331672)
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
    
    ## Create school dummy indicators (change to 10 for 10 school clusters)
    school_DI <- data.frame(model.matrix(simdataw$c_id ~ as.factor(simdataw$school) - 
        1, simdataw))
    names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x, y) paste(x, 
        seq(1, y), sep = "_"), "schoo_Ind", 40))
    school_DI <- school_DI[, 1:39]  #9 DIs if 10 school clusters
    
    ## Create a dataframe with complete variables
    datacomp <- cbind(Intercept = rep(1, nrow(simdataw)), simdataw[, names(simdataw) %in% 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7", "c_age", "c_gender", "c_nap1_z", 
            "p_sdq.2", "p_sdq.4", "p_sdq.6")], school_DI)
    
    ## Perform imputations without the random effects (SL imputation)
    imp1 <- jomo1con(Y = dataimp, X = datacomp, nimp = M, nburn = nburn, nbetween = NB)
    
    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        # 1. extract the mth imputed dataset
        datw <- imp1[imp1$Imputation == m, ]
        
        # 2. Attach the school variable and remove school indicator variables
        datw <- cbind(datw, simdataw$school)
        
        # 3. remove unwanted variables
        datw <- datw[, names(datw) %in% c("napscore_z.3", "napscore_z.5", "napscore_z.7", 
            "c_age", "c_gender", "c_dep.2", "c_dep.4", "c_dep.6", "c_ses", "c_nap1_z", 
            "id", "simdataw$school")]
        names(datw)[names(datw) == "simdataw$school"] <- "school"
        
        # 4. rename depression variables for reshape
        colnames(datw)[colnames(datw) == "c_dep.2"] <- "prev_dep.3"
        colnames(datw)[colnames(datw) == "c_dep.4"] <- "prev_dep.5"
        colnames(datw)[colnames(datw) == "c_dep.6"] <- "prev_dep.7"
        
        # 5. Reshape to long
        datL <- reshape(datw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
            c("napscore_z.3", "napscore_z.5", "napscore_z.7")), idvar = "id", v.names = c("prev_dep", 
            "napscore_z"), times = c(3, 5, 7), direction = "long")
        
        datL <- datL[order(datL$school, datL$id), ]
        
        # 6. save the dataset in a list
        mylist[[m]] <- datL
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + prev_dep * time + c_age + c_gender + 
            c_nap1_z + c_ses + (1 | school/id), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    
    # Compute CIs
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

# Save the results
write.xlsx(MVNIslwide_results.est, paste0("MVNIslwide_results", datnum, ".est.xlsx"))
write.xlsx(MVNIslwide_results.sd, paste0("MVNIslwide_results", datnum, ".sd.xlsx"))
write.xlsx(MVNIslwide_results.RE, paste0("MVNIslwide_results", datnum, ".RE.xlsx"))
write.xlsx(MVNIslwide_results.CI, paste0("MVNIslwide_results", datnum, ".CI.xlsx"))
write.xlsx(MVNIslwide_results.ICC, paste0("MVNIslwide_results", datnum, ".ICC.xlsx"))

###############################################
###               SMC-JM-2L-DI              ###
###############################################

## commands used for running the simulations in parallel
args <- commandArgs(trailingOnly = TRUE)  #These arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum <- as.numeric(args[1])

T1 <- list.files(pattern = "*.csv")

s2 <- seq(20, 1000, by = 20)
s1 <- s2 - 19

temp <- T1[s1[datnum]:s2[datnum]]

# setting random seeds for the parallely run replications

# MAR-CATS, 40 school clsuters
set.seed(6718920)
seed <- sample(1e+08, 50, replace = F)[datnum]

# MAR-CATS, 10 school clsuters
# set.seed(890015672)
# seed=sample.int(1e8,50,replace = F)[datnum] 

# MAR-inflated, 40 school clsuters
# set.seed(7902167)
# seed=sample.int(1e8,50,replace = F)[datnum] 

# MAR-inflated, 10 school clsuters
# set.seed(89995672)
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
    simdataL$c_id <- as.numeric(simdataL$c_id)
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Rearrange the data set for imputations
    
    ##Create the previous wave SDQ value
    simdataL <- slide(simdataL, Var = "p_sdq", GroupVar = "c_id", slideBy = -1)
    simdataL <- simdataL[, !names(simdataL) %in% c("p_sdq")]
    colnames(simdataL)[colnames(simdataL) == "p_sdq-1"] <- "prev_sdq"
    
    ##Rename the depression values at waves 2, 4 and 6 for reshape
    simdataw <- reshape(simdataL, v.names = c("napscore_z", "prev_sdq", "c_dep"), 
        idvar = "c_id", timevar = "wave", direction = "wide")
    
    colnames(simdataw)[colnames(simdataw) == "c_dep.2"] <- "prev_dep.3"
    colnames(simdataw)[colnames(simdataw) == "c_dep.4"] <- "prev_dep.5"
    colnames(simdataw)[colnames(simdataw) == "c_dep.6"] <- "prev_dep.7"
    
    simdataw <- simdataw[, !names(simdataw) %in% c("napscore_z.2", "prev_sdq.2", 
        "napscore_z.4", "prev_sdq.4", "napscore_z.6", "prev_sdq.6", "c_dep.3", "c_dep.5", 
        "c_dep.7")]
    
    ## Reshape back to long
    simdataL <- reshape(simdataw, varying = list(c("prev_dep.3", "prev_dep.5", "prev_dep.7"), 
        c("napscore_z.3", "napscore_z.5", "napscore_z.7"), c("prev_sdq.3", "prev_sdq.5", 
            "prev_sdq.7")), idvar = "c_id", v.names = c("prev_dep", "napscore_z", 
        "prev_sdq"), times = c(3, 5, 7), direction = "long")
    
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Recode sex as a factor
    simdataL$c_gender <- as.factor(simdataL$c_gender)
    
    ## Set number of imputations and number of burn-in iterations
    M <- 20
    nburn <- 500
    NB <- 10
    
    ## create school dummy indicators (change to 10 for 10 school clusters)
    school_DI <- data.frame(model.matrix(simdataL$c_id ~ as.factor(simdataL$school) - 
        1, simdataL))
    names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x, y) paste(x, 
        seq(1, y), sep = "_"), "schoo_Ind", 40))
    school_DI <- school_DI[, 1:39]  ##9 DIs if 10 school clusters
    
    ## Define a dataframe with all variables required
    data_jomo <- cbind(simdataL[!names(simdataL) %in% c("child", "school")], school_DI)
    
    ## Specify the levels of each variable in the imputation model
    mylevel <- c(2, 2, 2, 2, 2, 1, 1, 1, 1, rep(2, times = 39))  ##9 DIs if 10 school clusters
    
    ## Formula of the substantive lmer model (change accordingly for 10 school
    # clusters)
    formula <- as.formula(napscore_z ~ prev_dep + time + prev_dep * time + c_age + 
        c_gender + c_ses + c_nap1_z + schoo_Ind_1 + schoo_Ind_2 + schoo_Ind_3 + schoo_Ind_4 + 
        schoo_Ind_5 + schoo_Ind_6 + schoo_Ind_7 + schoo_Ind_8 + schoo_Ind_9 + schoo_Ind_10 + 
        schoo_Ind_11 + schoo_Ind_12 + schoo_Ind_13 + schoo_Ind_14 + schoo_Ind_15 + 
        schoo_Ind_16 + schoo_Ind_17 + schoo_Ind_18 + schoo_Ind_19 + schoo_Ind_20 + 
        schoo_Ind_21 + schoo_Ind_22 + schoo_Ind_23 + schoo_Ind_24 + schoo_Ind_25 + 
        schoo_Ind_26 + schoo_Ind_27 + schoo_Ind_28 + schoo_Ind_29 + schoo_Ind_30 + 
        schoo_Ind_31 + schoo_Ind_32 + schoo_Ind_33 + schoo_Ind_34 + schoo_Ind_35 + 
        schoo_Ind_36 + schoo_Ind_37 + schoo_Ind_38 + schoo_Ind_39 + prev_sdq + (1 | 
        c_id))
    
    ## Running the imputations
    imp5 <- jomo.lmer(formula, data_jomo, level = mylevel, nimp = M, nburn = nburn, 
        nbetween = NB)
    
    ## Analysis
    mylist <- list()
    
    for (m in 1:M) {
        
        # 1. Extract the ith imputed data set and attach the school variable
        imputed <- imp5[imp5$Imputation == m, ]
        imputed <- imputed[, names(imputed) %in% c("prev_dep", "c_age", "c_gender", 
            "c_ses", "c_nap1_z", "time", "napscore_z", "clus")]
        imputed$school <- simdataL$school
        
        # 2. save the dataset in a list for analysis
        mylist[[m]] <- imputed
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + time + prev_dep * time + c_age + c_gender + 
            c_nap1_z + c_ses + (1 | school/clus), data = d)
    })
    
    ## Combine the estimates
    MI_est <- testEstimates(mods, var.comp = TRUE)
    
    ## Compute CIs
    CI <- matrix(NA, 7, 2)
    Conf <- confint(MI_est)
    CI[, 1] <- as.vector(Conf[2:8, 1])
    CI[, 2] <- as.vector(Conf[2:8, 2])
    SMCJMDI_results.CI <- cbind(SMCJMDI_results.CI, CI)
    
    ## Store the estimates
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

# Save the results
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

# Setting random seeds for the parallely run replications

# MAR-CATS, 40 school clsuters
set.seed(890015672)
seed <- sample(1:1e+05, 10)[datnum]

# MAR-CATS, 10 school clsuters
# set.seed(10992)
# seed=sample.int(1e8,10,replace = F)[datnum] 

# MAR-inflated, 40 school clsuters
# set.seed(79019093)
# seed=sample.int(1e8,10,replace = F)[datnum] 

# MAR-inflated, 10 school clsuters
# set.seed(8001672)
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
    simdataL$c_id <- as.numeric(simdataL$c_id)
    simdataL <- simdataL[order(simdataL$school, simdataL$c_id), ]
    
    ## Create the previous wave deperssion symptoms
    simdataL <- slide(simdataL, Var = "c_dep", GroupVar = "c_id", slideBy = -1)
    
    colnames(simdataL)[colnames(simdataL) == "c_dep-1"] <- "prev_dep"
    
    simdataL$c_dep <- NULL
    
    ## Create previous wave SDQ variable
    simdataL <- slide(simdataL, Var = "p_sdq", GroupVar = "c_id", slideBy = -1)
    
    colnames(simdataL)[colnames(simdataL) == "p_sdq-1"] <- "prev_sdq"
    
    simdataL$p_sdq <- NULL
    
    ## Remove unwanted waves and variables
    simdataL <- subset(simdataL, wave != 2 & wave != 4 & wave != 6)
    simdataL <- simdataL[, !names(simdataL) %in% c("child")]
    
    ## Recode categorical and binary varaible with numeric codes 0,1,2,...
    simdataL$c_gender <- ifelse(simdataL$c_gender == "male", 1, 0)
    
    simdataL$school <- simdataL$school - 1
    
    ## Model specification
    iter <- 2000
    burnin <- 1000  ##burn in for 1000 iteration and one imputed dataset is
    Nimp <- 20  ##saved every 100th iteration=20imps
    
    ## Substantive model formula
    Y_formula <- napscore_z ~ prev_dep + wave + prev_dep * wave + c_age + c_ses + 
        c_gender + c_nap1_z + prev_sdq + school + (1 | c_id)
    
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
    mod1 <- mdmb::frm_fb(simdataL, dep, ind, iter = iter, burnin = burnin, Nimp = 20, 
        aggregation = TRUE)
    
    ## Convert output into list of imputed datasets
    datlist <- mdmb::frm2datlist(mod1)
    
    ## Analysis
    mylist <- list()
    
    for (j in 1:Nimp) {
        
        # 1. Extract the mth imputed dataset
        datL <- datlist[[j]]
        
        # 2. Reattach the school cluster indicator
        datL$school <- simdataL$school
        
        datL <- datL[, names(datL) %in% c("c_id", "wave", "napscore_z", "c_age", 
            "c_gender", "c_nap1_z", "prev_dep", "prev_sdq", "c_ses", "school")]
        
        # 3. Save the dataset in a list
        mylist[[j]] <- datL
        
    }
    
    ## Fit the analysis of interest on the imputed datasets
    mods <- lapply(mylist, function(d) {
        lmer(napscore_z ~ prev_dep + wave + prev_dep * wave + c_age + as.factor(c_gender) + 
            c_nap1_z + c_ses + (1 | school/c_id), data = d)
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
colnames(SMCSMDI_results.est) <- c(seq(1:length(data)))

rownames(SMCSMDI_results.sd) <- rows
colnames(SMCSMDI_results.sd) <- c(seq(1:length(data)))

rownames(SMCSMDI_results.RE) <- c("level 3", "level 2", "level 1")
colnames(SMCSMDI_results.RE) <- c(seq(1:length(data)))

rownames(SMCSMDI_results.ICC) <- c("level 3", "level 2")
colnames(SMCSMDI_results.ICC) <- c(seq(1:length(data)))

rownames(SMCSMDI_results.CI) <- rows
colnames(SMCSMDI_results.CI) <- c(rep(1:length(data), each = 2))

# Save the results
write.xlsx(SMCSMDI_results.est, paste0("SMCSMDI_results.est", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.sd, paste0("SMCSMDI_results.sd", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.RE, paste0("SMCSMDI_results.RE", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.CI, paste0("SMCSMDI_results.CI", datnum, ".xlsx"))
write.xlsx(SMCSMDI_results.ICC, paste0("SMCSMDI_results.ICC", datnum, ".xlsx"))
