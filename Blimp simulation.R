#################################################################################
#                                                                               #
# Running simulations for SMC-JM-3L in Blimp via the External simulation method #                                                                      #
#  This script accompanys the paper:                                             #
# "Evaluation of approaches for accommodating interactions and non-linear terms #
#   in multiple imputation of incomplete three-level data"                      #
#                                                                               #
# Rushani Wijesuriya, 11 October 2020                                           #
#                                                                               #
#################################################################################


#more instructions on running simulations in Blimp externally can be found in the Blimp user guide(v 1.0.6) 

## Running a simulation via External Blimp method.
# Note no spaces should be in the directory paths
# All output from Blimp Is saved in a list: outputs
# Does not work with 'sep' subcommand


# Set a Seed (MAR-CATS, 40 school clusters)
seed <- 37298372
#Seeds used for the different simulation scenarios

#Analysis model 1:Interaction between the time-varying exposure and time
# seed <- 37298372#MAR-CATS, 10 school clusters
# seed <- 38298372 #MAR-inflated, 40 school clusters
# seed <- 39298372  #MAR-inflated, 10 school clusters

#Analysis model 2: Interaction between the time-varying exposure and a time-fixed baseline variable 
# seed <- 77258472 #MAR-CATS, 40 school clusters
# seed <- 88210726 #MAR-CATS, 10 school clusters
# seed <- 9865241 #MAR-inflated, 40 school clusters
# seed <- 34528910  #MAR-inflated, 10 school clusters

#*Analysis model 3: Quadratic term in the exposure
# seed <- 32498352 #MAR-CATS, 40 school clusters
# seed <- 982671478 #MAR-CATS, 10 school clusters
# seed <- 6831790 #MAR-inflated, 40 school clusters
# seed <- 24378198  #MAR-inflated, 10 school clusters

# Path to Blimp
blimpPath <- "~/Desktop/Blimp/blimp.exe"
# Specify Path to Syntax File
inputPath <- "~/blimpimps/blimpimp2.imp"

# Specify Path to Data Folder
# Only the data files should be in folder.
dataPath <- "~/Blimpdata"
# Specify Path to Output Folder
outputPath <- "~/blimpimps"
# Specify Names of Imputation Data
# Use a * to represent where the name of data file.
# E.g., Data1.csv will give you impData1.csv
impsData <- "imp*.csv"

#######################################################
## PROGRAM BEGINS HERE
# Get file names
dataFiles <- list.files(path =dataPath)
dataFilesPath <-list.files(path = dataPath, full.names = T)
# Calculate total number of reps.
repNumber <- length(dataFiles)
# Set seed
set.seed(seed)
# Generate list of seeds
seeds <- sample.int(1e10, repNumber,replace = F)
# Execute Simulation
outputs <- lapply(seq_along(seeds),function(x){
  repla <- gsub("\\..+","",dataFiles[x],perl=T)
  fileP <- gsub(".+\\.","",dataFilesPath[x],perl=T)
  outfile <- paste0(fileP,gsub("\\*",repla,impsData,perl=T))
 
  out= system(paste(blimpPath,inputPath,'-o',outfile,'-s',
                    seeds[x],'-d',dataFilesPath[x]),intern = T)
  
  return(list(dataFiles[x],out))
})