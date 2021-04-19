#######################################################################################
#   Filename    :   Figures.R 
#                   
#   Project     :   This script will generate the Figures in the Main text (Figures 1-3)and in
#                   the supplementary section (S1-S3) of
#                   BiomJ article "Evaluation of approaches for accommodating interactions and 
#                   non-linear terms in multiple imputation of incomplete three-level data"  
#   Author      :   Rushani Wijesuriya                                                                
#   Date        :   07.09.2020
########################################################################################

#************************************Notes for running this script**********************

# 1. Insert the path to the  Code_and_Data folder in line 73 
# 2. Run the script from begining to end, the tables will be saved  in the clus10 folder under 
#    each missing data mechanism under each analysis model
#    For example tables S2 and S3 in folder Interaction-> MAR-CATS -> Clus10 and tables S4 and S5 
#    in the folder Interaction -> MAR-inflated -> Clus10 etc.


rm(list = ls())

library(readxl)
library(WriteXLS)
library(rsimsum)
library(ggplot2)
library(tidyverse)
library(gtable)  
library(grid) 
library(ggpubr)
require(gridExtra)
library(stringr)


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
    
    datasets=list()
    dataSE.Clus=list()
    datasets.RB=list()
    
    for (j in 1:length(Clus)) {
      
      
      ## Set path to the working directory with results (See README file)
      setwd(file.path("..", "/Code_and_Data/Intermediate_results/",
                      Analysis[[k]], MD.mech[[i]], Clus[[j]]))
      
      ## Import all the results files in the folder
      temp <- list.files(pattern = glob2rx("*.xlsx"))
      
      ############################################### Regression coefficient estimates plot ###
      
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
      
  
      # Extract the main effect and interaction effect estimates
      dep=c()
      inter=c()
      
      for (m in 1:length(temp1)){
        data=read_excel(temp1[m])
        data=data[-1]
        bias1=as.numeric(as.vector(as.matrix(data[1,])))-beta1[k]
        dep=c(dep,bias1)
        
        bias2=as.numeric(as.vector(as.matrix(data[7,])))-beta3[k]
        inter=c(inter,bias2) }
      
      Biases=c(dep,inter)
      
      Method=rep(Method.list[[k]], each=1000,times=2)
      
      Dataset=rep(1:1000,length(temp1)*2)
      
      Estimate=rep(c(paste("Main effect (",beta1[k],")") ,paste("Interaction effect (",beta3[k],")")),each=length(temp1)*1000)
      
      Clus.size=rep(Clus[j])
      reg_data=data.frame(Method,Dataset,Estimate,Biases,Clus.size)
      
      datasets[[j]]=reg_data
      
      ###############################SE Plot#####################################
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
      
      
      # 1. Empirical SE
      emp <- cbind(get_data(S1, stats = c("empse")), order = od)
      
      # 2. Model based SE
      mod <- cbind(get_data(S1, stats = c("modelse")), order = od)
      
      prevdeplist <- list(emp, mod)
      sorted.est1 <- lapply(prevdeplist, function(d) {
        d[order(d$order), ]
      })
      
      ## For the interaction effect estimate
      S2 <- simsum(data = MSIM_inter, estvarname = "b", true = beta3[k], se = "se", 
                   methodvar = "Method", ref = "JM-1L-DI-wide")
      
      # 1. Empirical SE
      emp1 <- cbind(get_data(S2, stats = c("empse")), order = od)
      
      # 2. Model based SE
      mod1 <- cbind(get_data(S2, stats = c("modelse")), order = od)
      
      interlist <- list(emp1, mod1)
      sorted.est2 <- lapply(interlist, function(d) {
        d[order(d$order), ]
      })
      
      
      Method=rep(Method.list[[k]])
      
      Clus.size=rep(Clus[j])
      
      ##for prev_dep
      empSE1=round(sorted.est1[[1]][,2],6)
      modSE1=round(sorted.est1[[2]][,2],6)
      MCSE1=round(sorted.est1[[1]][,3],6)
      
      
      ##for interaction term 
      empSE2=round(sorted.est2[[1]][,2],6)
      modSE2=round(sorted.est2[[2]][,2],6)
      MCSE2=round(sorted.est2[[1]][,3],6)
      
      Empse=c(empSE1,empSE2)
      Modse=c(modSE1,modSE2)
      MCSE=c(MCSE1,MCSE2)
      Estimate=rep(c(paste("Main effect (",beta1[k],")") ,paste("Interaction effect (",beta3[k],")")),each=length(temp1))
      
      SE_data=data.frame(Method,Estimate,Empse,Modse,MCSE,Clus.size)
      dataSE.Clus[[j]]=SE_data
      
      ####Variance component estimates 
      
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
      L1 <- L2 <- L3 <- c()
      
      
      
      for (m in 1:length(temp3)) {
        data <- read_excel(temp3[m], col_types = "numeric")
        data <- data[-1]
        
        bias1=(as.vector(as.matrix(data[3,]))^2)-VC1
        bias2=(as.vector(as.matrix(data[2,]))^2)-VC2
        bias3=(as.vector(as.matrix(data[1,]))^2)-VC3
        
        L1=c(L1,bias1)
        L2=c(L2,bias2)
        L3=c(L3,bias3)
      }
      
      Biases= c(L3,L2,L1)
      
      Method=rep(Method.list[[k]],each=1000,times=3)
      
      Clus.size=rep(Clus[j])
      
      Estimate=rep(c("Level 3 VC (0.05)","Level 2 VC (0.45)",
                     "Level 1 VC (0.50)"),each=length(temp1)*1000)
      
      Dataset=rep(1:1000,length(temp1)*3)
      
      MSIM_var=data.frame(Method,Dataset,Estimate,Biases,Clus.size)
      
      datasets.RB[[j]]=MSIM_var
     
    }
    
    
    
    #######Plot for the reg coefficients
    df=rbind(datasets[[1]],datasets[[2]])
    
    df$Method=factor(df$Method,levels=Method.list[[k]])
    
    
    levels(df$Clus.size)=c("40 school clusters","10 school clusters")
    
    
    df$Estimate=factor(df$Estimate,levels=c(paste("Main effect (",beta1[k],")") ,paste("Interaction effect (",beta3[k],")")))

    
    df$Estimate = as.factor(str_wrap(df$Estimate, width = 11))
    
    P_est=ggplot(df, aes(x=Method,y=Biases,fill=Method))+
      geom_boxplot() +facet_grid(Estimate~Clus.size,scales="free_y")+
      geom_hline(yintercept = 0, linetype=1, color = "black", size= .5)+theme(axis.text.x=element_blank())+
      ylab("Estimate-True value")+ theme( axis.ticks.x=element_blank())+theme(axis.title.x = element_blank())+
      theme(strip.text.x = element_text(size = 10, face = "bold"),strip.text.y = element_text(size = 10,face = "bold"))
    
    
    ##Plot for the SEs
    df1=rbind(dataSE.Clus[[1]],dataSE.Clus[[2]])
    
    ###generate the uppper and lower CI based on MCE(1.95*mce)
    df1$LCI=df1$Empse-(1.95*df1$MCSE)
    df1$UCI=df1$Empse+(1.95*df1$MCSE)
    
    df1$Method=factor(df1$Method,levels=Method.list[[k]]) 
    
    df1$Estimate=factor(df1$Estimate,levels=c(paste("Main effect (",beta1[k],")") ,paste("Interaction effect (",beta3[k],")")))
    
    levels(df1$Clus.size)=c("40 school clusters","10 school clusters")
    
    df1$Estimate = str_wrap(df1$Estimate, width = 11)
    
    P_se <- ggplot(df1, aes(x=Method, y=Empse,ymin = LCI, ymax = UCI, color=Method )) + geom_pointrange()+
      facet_grid(Estimate~Clus.size,scales="free_y")+theme(axis.text.x=element_blank())+
      geom_point(aes(x=as.numeric(Method)+0.3, y=Modse),data = df1, fill="white",shape=21)+ylab("Standard errors") +
      theme(axis.title.x = element_blank())+
      theme(strip.text.x = element_text(size = 10, face = "bold"),strip.text.y = element_text(size = 10, face = "bold"))
    
    
    ##Plot for the VCs
    df2=rbind(datasets.RB[[1]],datasets.RB[[2]])
    
    df2$Method=factor(df2$Method,levels=Method.list[[k]])
    
    
    levels(df2$Clus.size)=c("40 school clusters","10 school clusters")
    
    df2$Estimate=factor(df2$Estimate,levels=c("Level 3 VC (0.05)","Level 2 VC (0.45)",
                                            "Level 1 VC (0.50)"))
    
    df2$Estimate = str_wrap(df2$Estimate, width= 8)
    
    P_vc= ggplot(df2, aes(x=Method,y=Biases,fill=Method))+
      geom_boxplot() +facet_grid(Estimate~Clus.size, scales="free_y")+
      geom_hline(yintercept = 0, linetype=1, color = "black", size= .5)+theme(axis.text.x=element_blank())+
      ylab("Estimate-True value")+ theme( axis.ticks.x=element_blank())+theme(axis.title.x = element_blank())+
      theme(strip.text.x = element_text(size = 10, face = "bold"),strip.text.y = element_text(size = 10, face = "bold"))
    
    
    #combine plots
    Plot <- ggarrange(P_est,P_se,P_vc, nrow=3, labels = c("a","b","c"))
  
    tiff(paste0("Figure",k,".",MD.mech[[i]],".tiff"), width= 7000, height= 7658, units="px", res=800)
    plot(Plot)
    dev.off()
    
    
  }
  
}

