TITLE OF MANUSCRIPT: EVALUATION OF APPROACHES FOR ACCOMMODATING INTERACTIONS AND NON-LINEAR TERMS IN MULTIPLE IMPUTATION OF INCOMPLETE THREE-LEVEL DATa

Authors: Rushani Wijesuriya*, Margarita Moreno-Betancur, John B. Carlin, Anurika P. De Silva, Katherine J. Lee
*responsible for writing the code
e-mail address: rushani.wijesuriya@mcri.edu.au
Configurations of software: R version 3.6.1 for Windows (64-bit x86-64) and Blimp studio version 1.0.6 
Special hardware or software requirements:  None

How to execute the code (Work-flow):

******************Step 1: Complete and missing data generation

data_generation.do contains the Stata program ("simul") for simulating the complete and missing data. 
This can be used to simulate the data under the  12 different simulation scenarios considered in the manuscript (3 substantive analysis models x 2 cluster sizes x 2 missingness mechanisms). 
Generate 12 folders to store the data (i.e. for each simulation scenario). This could be generated in a factorial manner as follows.  

For example:

1. For the simulation scenario: analysis model 1, missing data mechanism MAR-CATS, 40 clusters
/Interaction
./MAR-CATS
../Clus40
.../Data

2. For the simulation scenario: analysis model 1, missing data mechanism MAR-inflated, 40 clusters
/Interaction
./MAR-inflated
../Clus40
.../Data

Under each simulation scenario, the user should change the file path in the Stata script when generating the data to store them in a relevant folder. Currently it is cd "C:\temp". 
Data sets (1 or 2) were omitted under some simulation scenarios due to non-convergence of the imputation model in SMC-JM-3L. The iteration relevant to the omitted datasets are included as comments in the Stata script.  

*******************Step 2: Application of the MI methods on the simulated data

*************************
* METHODS EXECUTED IN R *
*************************
The MI_methods_analysis1.R, MI_methods_analysis2.R, and MI_methods_analysis3.R contains the R code for the MI approaches (executed in R) considered under the three substantive analysis models. 

While the MI methods executed under each analysis model are given sequentially in the R scripts for convenience, each MI approach were run in parallel on the 1000 simulated datasets in the simulation study using different seeds as given in the script. 
For time efficiency some MI approaches were also run on fractions of all repetitions in parallel. The relevant fractions and the random seeds used for the parallel runs of simulations under each simulation scenario are given as comments  under each MI approach. 

Using the R scripts, run each of the MI approach for each simulation scenario and save the parameter estimates of interest obtained for the 1000 replications (simulated datasets). 

****************************
* METHOD EXECUTED IN BLIMP *
****************************

*NOTE:There is no direct way to execute simulations on the Blimp application but this can be done externally using R.
The Blimp_simulation.R contains the R code which can be used for this. The R script make changes to the the input file path and the output file path in the Blimp syntax file blimpimputation.imp via command line arguments. Therefore these files should be used together. 

To run the simulations: 

**Step (i)- Recode data to avoid import fails

	Using the Blimp_prep.R file recode all missing values(as 999) and categorical variables as numerical variables to generate a new set of 1000 datasets for each simulation scenario. This is because special characters in the data such as NA values and string values can cause the data import to fail in Blimp.
	The input path should be changed to the folder where the original data for each simulation scenario is saved. The output path should be specified to a new folder to save the newly generated datasets after transformation. For example a new folder called 'Blimpdata' can be generated in each folder generated for the different simulation scenarios in Step 1.

	For example:

	For the simulation scenario: analysis model 1, missing data mechanism MAR-CATS, 40 clusters
	/Interaction
	./MAR-CATS
	../Clus40
	.../Blimpdata

 NOTE: Make sure that this new folder only contains the newly generated datasets when running steps(ii) and (iii) below
 
**Step (ii)- Change MODEL command in the Blimp syntax file to specify the correct imputation model according to the substantive analysis 

	Under each analysis model, change the MODEL command line  in the blimpimputation.imp file accordingly. Save this file in a new folder where the imputations generated will be saved. The new folder can be generated in each folder for the different simulation scenarios generated in Step 1.

	For example:

	For the simulation scenario: analysis model 1, missing data mechanism MAR-CATS, 40 clusters
	/Interaction
	./MAR-CATS
	../Clus40
	.../blimpimps

	The multiple imputations from a single dataset will be saved in a single file in a stacked format (as csvimpdata1.csv,..., csvimpdata1000.csv). 

**Step (iii) - Specify the file paths in the Blimp_simulation.R script to run the simulations

	In the R script Blimp_simulation.R: 
	
	Change the "blimpPath" to where the blimp application is saved locally in the computer
	Change the "inputPath" to where the Blimp syntax file from Step(ii) is saved (according to the above example, to the "blimpimps" folder)
	Change the "dataPath" to where the transformed data in step (i) is saved (according to the above example, to the "Blimpdata" folder). 
	Finally change the "outputPath" to the new folder created in Step (ii) (according to the above example, to the "blimpimps" folder). 

	Seeds used for each simulation scenario are given as comments. Use these to reproduce the results in the manuscript. 
    Save this R file in the  same folder as where imputations will be saved (blimpimps folder) and run the R script from beginning to end. 

**Step (iv)- Analyze and pool results from the multiply imputed datasets

After the simulations are run, use the Blimp_post_analysis.R script to conduct the analysis on the imputed datasets and pool the results to obtain the parameter estimates of interest for the 1000 replications.
In the Blimp_analysis.R script, set the  working directory to where the imputations are saved(for example: blimpimps folder). 

*Note: Execution of all MI methods on the different simulation scenarios as detailed in step 2 takes very long to run (at least several weeks). Therefore, intermediate results from Step 2 are available in the intermediate_results folder. These results should be used for the remaining steps.


*******************Step 3: Summarizing the performance of the different approaches under different simulation scenarios and produce tables for the manuscript

To summarize the results obtained for each MI approach and generate the tables S2- S13 in the supplementary file, use the R script Results_summarize.R. The user should change the file path to the folder where the Code_and_Data folder is saved. Currently it is setwd(file.path("..","/Code_and_Data/Intermediate_results/",Analysis[[k]], MD.mech[[i]], Clus[[j]]))
Run the script from beginning to end, the tables will be saved  in the clus10 folder under each missing data mechanism under each analysis model both in csv and pdf format. 

*******************Step 4: Produce figures for the manuscript

To generate the Figures 1-3 and S1-S3 use the R script Figures.R. The user should change the file path to the folder where the Code_and_Data folder is saved. Currently it is setwd(file.pat("..","/Code_and_Data/Intermediate_results/", Analysis[[k]], MD.mech[[i]], Clus[[j]])).The figures tables will be saved  in the clus10 folder under each missing data mechanism under each analysis model


The CATS dataset analyzed as the case study in section 5 is not publicly available due to ethics requirements. The simulation study in section 4 was designed to closely mimic the CATS data.