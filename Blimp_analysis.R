#################################################################################
#                                                                               #
#  Post imputation analysis and pooling for SMC-JM-3L in Blimp                  #
#  for the analysis models in the paper:                                        #
# "Evaluation of approaches for accommodating interactions and non-linear terms #
#   in multiple imputation of incomplete three-level data"                      #
#                                                                               #
# Rushani Wijesuriya, 11 October 2020                                           #
#                                                                               #
#################################################################################

# Required packages
library(mitml)  #for fitting the lmer model
library(lme4)   # for pooling the results
library(writexl) #for writing/exporting the results

##clear the workspace
rm(list = ls())

##load the simulated data 
temp = list.files(pattern=glob2rx("*.csv"))

data=lapply(temp,read.csv,header=FALSE)

blimp_results.est=matrix(NA,nrow=7,ncol=length(temp))
blimp_results.sd=matrix(NA,nrow=7,ncol=length(temp))
blimp_results.RE=matrix(NA,nrow=3,ncol=length(temp))
blimp_results.ICC=matrix(NA,nrow=2,ncol=length(temp))
blimp_results.CI=c()


for (i in 1:length(temp)){
  
  impdata=data[[i]]
  names(impdata)= c("impno","c_id","time","napscore_z","school",
                    "c_age","c_gender","c_ses","c_nap1_z","prev_dep","prev_sdq")
  
  impdata$c_gender=as.factor(impdata$c_gender)
  impdata$c_id=as.factor(impdata$c_id)
  
  mylist=list()
  
  for (m in 1:20){
    
    mylist[[m]]<-impdata[impdata$impno==m,]
    
  }
  
  
  #fit the analysis of interest on the imputed datasets 
  
  #(i)for the analysis model involving an interaction between the time-varying exposure and time
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*time+c_age+
                                             as.factor(c_gender)+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  #(ii)for the analysis model involving an interaction between the time-varying exposure and a time-fixed confounder
  # mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
  #                                            as.factor(c_gender)+c_nap1_z+c_ses
  #                                          +(1|school/c_id), data = d)} )
  
  #(iii)for the analysis model involving a quadratic effect of the the time-varying exposure
  # mods <- lapply(mylist,function(d) {d$prev_dep2=d$prev_dep*d$prev_dep
  #                                   lmer( napscore_z~prev_dep+time+c_age+
  #                                   as.factor(c_gender)+c_nap1_z+c_ses+prev_dep2
  #                                    +(1|school/c_id), data = d)
  #                                     } )
 
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE,df.com=NULL)
  
  
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  blimp_results.CI=cbind(blimp_results.CI,CI)
  
  
  #store the estimates
  blimp_results.est[,i]=MI_est$estimates[2:8,1]
  blimp_results.sd[,i]=MI_est$estimates[2:8,2]
  blimp_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  blimp_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  blimp_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  blimp_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  blimp_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
  
}

rownames(blimp_results.est)=c("prev_dep","time","c_age","c_gender","c_nap1_z","c_ses","inter")
colnames(blimp_results.est)=c(seq(1:length(temp)))

rownames(blimp_results.sd)=c("prev_dep","time","c_age","c_gender","c_nap1_z","c_ses","inter")
colnames(blimp_results.sd)=c(seq(1:length(temp)))

rownames(blimp_results.RE)=c("level 3","level 2","level 1")
colnames(blimp_results.RE)=c(seq(1:length(temp)))

rownames(blimp_results.ICC)=c("level 3","level 2")
colnames(blimp_results.ICC)=c(seq(1:length(temp)))

rownames(blimp_results.CI)=c("prev_dep","time","c_age","c_gender","c_nap1_z","c_ses","inter")
colnames(blimp_results.CI)=c(rep(1:length(temp), each=2))

##save the results
write_xlsx(as.data.frame(cbind(rownames(blimp_results.est),blimp_results.est)),"BLIMP_results.est.xlsx")
write_xlsx(as.data.frame(cbind(rownames(blimp_results.sd),blimp_results.sd)),"BLIMP_results.sd.xlsx")
write_xlsx(as.data.frame(cbind(rownames(blimp_results.RE),blimp_results.RE)),"BLIMP_results.RE.xlsx")
write_xlsx(as.data.frame(cbind(rownames(blimp_results.ICC),blimp_results.ICC)),"BLIMP_results.ICC.xlsx")
write_xlsx(as.data.frame(cbind(rownames(blimp_results.CI),blimp_results.CI)),"BLIMP_results.CI.xlsx")