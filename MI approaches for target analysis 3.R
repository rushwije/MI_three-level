#################################################################################
#                                                                               #
#  MI approaches for the analysis model involving a quadratic effect of the     #
#  the time-varying exposure evaluated in the paper                             #
# "Evaluation of approaches for accommodating interactions and non-linear terms #
#   in multiple imputation of incomplete three-level data"                      #
#                                                                               #
# Rushani Wijesuriya, 07 september 2020                                         #
#                                                                               #
#################################################################################

##clear the workspace
rm(list = ls())

##load the required packages
require(lme4)   #for fitting the lmer model
require(jomo)   # for single-level and two-level JM and two-level SMC JM
require(mitml)  # for pooling the results
require(mice)   #for single and two-level FCS
require(xlsx)   #for writing/exporting the results
require(DataCombine) #for generating the lagged wave variables
require(mdmb)     #for two-level SMC-SM

##set working dir to where data is stored
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##load the simulated data 
data=lapply(list.files(pattern=glob2rx("*.csv")),read.csv)

##############################################################################
#                      JM-1L-DI-wide                                         #
##############################################################################

MVNIslwide_results.est=matrix(NA,nrow=7,ncol=length(data))
MVNIslwide_results.sd=matrix(NA,nrow=7,ncol=length(data))
MVNIslwide_results.RE=matrix(NA,nrow=3,ncol=length(data))
MVNIslwide_results.ICC=matrix(NA,nrow=2,ncol=length(data))
MVNIslwide_results.CI=c()


for (i in 1:length(data)){
  
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
  
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3",
                                           "c_dep.5","c_dep.7")]
  
  ##Set number of imputations and burn in and thining intervals 
  M=20
  nburn=1000
  NB=100
  
  ##create a dataframe with variables to be imputed
  myvars <- names(simdataw) %in% c("c_dep.2","c_dep.4","c_dep.6","c_ses") 
  dataimp=simdataw[myvars]
  
  ##create school dummy indicators
  school_DI=data.frame(model.matrix(simdataw$c_id~as.factor(simdataw$school)-1,
                                    simdataw))
  names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), 
                                                       "schoo_Ind",40))
  school_DI=school_DI[,1:39]
  
  ##create a dataframe with complete variables
  datacomp= cbind(Intercept=rep(1,nrow(simdataw)),
                  simdataw[,names(simdataw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                                 "c_age","c_gender",
                                                 "c_nap1_z","p_sdq.2","p_sdq.4","p_sdq.6"
                  )],
                  school_DI)
  
  ##perform imputations without the random effects (SL imputation)
  set.seed(2946)
  imp1<-jomo1con(Y=dataimp, X=datacomp, nimp=M,nburn=nburn,nbetween=NB)
  
  # Examine convergence
  # impCheck<-jomo1con.MCMCchain(Y=dataimp, X=datacomp, nburn=nburn)
  # plot(c(1:nburn),impCheck$collectbeta[1,1,1:nburn],type="l")
  
  ##Analysis
  mylist=list()
  
  for(m in 1:M)
  {
    datw<-imp1[imp1$Imputation==m,]
    
    #1. attach the school variable and remove school indicator variables
    datw=cbind(datw,simdataw$school)
    
    #2. remove unwanted variables
    datw=datw[,names(datw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                "c_ses","c_nap1_z","id","simdataw$school")]
  
    names(datw)[names(datw) == "simdataw$school"] <- "school"
    
    #3. rename depression variables for reshape
    colnames(datw)[colnames(datw)=="c_dep.2"] <- "prev_dep.3"
    colnames(datw)[colnames(datw)=="c_dep.4"] <- "prev_dep.5"
    colnames(datw)[colnames(datw)=="c_dep.6"] <- "prev_dep.7"
    
    #4. Reshape to long
    datL=reshape(datw,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                    c("napscore_z.3","napscore_z.5","napscore_z.7")),idvar="id", 
                 v.names=c("prev_dep","napscore_z"), times=c(3,5,7),direction= "long")
    
    datL <- datL[order(datL$school,datL$id),]    
    
    #5. save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  MVNIslwide_results.CI=cbind(MVNIslwide_results.CI,CI)
  
  
  #store the estimates
  MVNIslwide_results.est[,i]=MI_est$estimates[2:8,1]
  MVNIslwide_results.sd[,i]=MI_est$estimates[2:8,2]
  MVNIslwide_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  MVNIslwide_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  MVNIslwide_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  MVNIslwide_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  MVNIslwide_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
}

rows=c("prev_dep","time","c_age","c_gender","c_nap1_z","c_ses","interaction")

rownames(MVNIslwide_results.est)=rows
colnames(MVNIslwide_results.est)=c(seq(1:length(data)))

rownames(MVNIslwide_results.sd)=rows
colnames(MVNIslwide_results.sd)=c(seq(1:length(data)))

rownames(MVNIslwide_results.RE)=c("level 3","level 2","level 1")
colnames(MVNIslwide_results.RE)=c(seq(1:length(data)))

rownames(MVNIslwide_results.ICC)=c("level 3","level 2")
colnames(MVNIslwide_results.ICC)=c(seq(1:length(data)))

rownames(MVNIslwide_results.CI)=rows
colnames(MVNIslwide_results.CI)=c(rep(1:length(data),each=2))

##save the results
write.xlsx(MVNIslwide_results.est,"MVNIslwide_results.est.xlsx")
write.xlsx(MVNIslwide_results.sd,"MVNIslwide_results.sd.xlsx")
write.xlsx(MVNIslwide_results.RE,"MVNIslwide_results.RE.xlsx")
write.xlsx(MVNIslwide_results.CI,"MVNIslwide_results.CI.xlsx")
write.xlsx(MVNIslwide_results.ICC,"MVNIslwide_results.ICC.xlsx")

##############################################################################
#                      JM-1L-DI-wide-JAV                                     #
##############################################################################

MVNIslwideJAV_results.est=matrix(NA,nrow=7,ncol=length(data))
MVNIslwideJAV_results.sd=matrix(NA,nrow=7,ncol=length(data))
MVNIslwideJAV_results.RE=matrix(NA,nrow=3,ncol=length(data))
MVNIslwideJAV_results.ICC=matrix(NA,nrow=2,ncol=length(data))
MVNIslwideJAV_results.CI=c()

for (i in 1:length(data)){
  
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
 
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3","c_dep.5","c_dep.7")]
  
  ##Set number of imputations 
  M=20
  nburn=1000
  NB=100

  ##generate the interaction terms
  simdataw$c_dep.ses.2=simdataw$c_dep.2*simdataw$c_ses
  simdataw$c_dep.ses.4=simdataw$c_dep.4*simdataw$c_ses
  simdataw$c_dep.ses.6=simdataw$c_dep.6*simdataw$c_ses
  
  ##create a dataframe with variables to be imputed
  myvars <- names(simdataw) %in% c("c_dep.2","c_dep.4","c_dep.6","c_ses","c_dep.ses.2","c_dep.ses.4","c_dep.ses.6") 
  dataimp=simdataw[myvars]
  
  
  ##create school dummy indicators
  school_DI=data.frame(model.matrix(simdataw$c_id~as.factor(simdataw$school)-1,
                                    simdataw))
  names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), 
                                                       "schoo_Ind",40))
  school_DI=school_DI[,1:39]
  
  ##create a dataframe with complete variables
  datacomp= cbind(Intercept=rep(1,nrow(simdataw)),
                  simdataw[,names(simdataw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                                 "c_age","c_gender",
                                                 "c_nap1_z","p_sdq.2","p_sdq.4","p_sdq.6"
                  )],
                  school_DI)
  
  ##perform imputations without the random effects (SL imputation)
  set.seed(2946)
  imp1<-jomo1con(Y=dataimp, X=datacomp, nimp=M,nburn=nburn,nbetween=NB)
  
  #Examine convergence
  # impCheck<-jomo1con.MCMCchain(Y=dataimp, X=datacomp, nburn=nburn)
  # plot(c(1:nburn),impCheck$collectbeta[1,1,1:nburn],type="l")
  
  ##Analysis
  mylist=list()
  
  for(m in 1:M)
  {
    datw<-imp1[imp1$Imputation==m,]
    
    #1. attach the school variable and remove school indicator variables
    datw=cbind(datw,simdataw$school)
  
    #2. remove unwanted variables
    datw=datw[,names(datw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                "c_dep.ses.2","c_dep.ses.4","c_dep.ses.6",
                                "c_ses","c_nap1_z","id","simdataw$school")]
    colnames(datw)[colnames(datw)=="simdataw$school"] <- "school"
    
    #3. rename depression variables for reshape
    colnames(datw)[colnames(datw)=="c_dep.2"] <- "prev_dep.3"
    colnames(datw)[colnames(datw)=="c_dep.4"] <- "prev_dep.5"
    colnames(datw)[colnames(datw)=="c_dep.6"] <- "prev_dep.7"
    
    colnames(datw)[colnames(datw)=="c_dep.ses.2"] <- "inter.3"
    colnames(datw)[colnames(datw)=="c_dep.ses.4"] <- "inter.5"
    colnames(datw)[colnames(datw)=="c_dep.ses.6"] <- "inter.7"
    
    #4. Reshape to long
    
    datL=reshape(datw,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                "napscore_z.5","napscore_z.7"),
                                    c("inter.3","inter.5","inter.7"))
                 ,idvar="id", 
                 v.names=c("prev_dep","napscore_z","inter"), times=c(3,5,7),direction= "long")
    
    
    
    datL <- datL[order(datL$school,datL$id),]    
    
    #5. save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+c_age+
                                             c_gender+c_nap1_z+c_ses+inter
                                           +(1|school/id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  MVNIslwideJAV_results.CI=cbind(MVNIslwideJAV_results.CI,CI)
  
  #store the estimates
  MVNIslwideJAV_results.est[,i]=MI_est$estimates[2:8,1]
  MVNIslwideJAV_results.sd[,i]=MI_est$estimates[2:8,2]
  MVNIslwideJAV_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  MVNIslwideJAV_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  MVNIslwideJAV_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  MVNIslwideJAV_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  MVNIslwideJAV_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])

}

rownames(MVNIslwideJAV_results.est)=rows
colnames(MVNIslwideJAV_results.est)=c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.sd)=rows
colnames(MVNIslwideJAV_results.sd)=c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.RE)=c("level 3","level 2","level 1")
colnames(MVNIslwideJAV_results.RE)=c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.ICC)=c("level 3","level 2")
colnames(MVNIslwideJAV_results.ICC)=c(seq(1:length(data)))

rownames(MVNIslwideJAV_results.CI)=rows
colnames(MVNIslwideJAV_results.CI)=c(rep(1:length(data),each=2))

##save the results
write.xlsx(MVNIslwideJAV_results.est,"MVNIslwideJAV_results.est.xlsx")
write.xlsx(MVNIslwideJAV_results.sd,"MVNIslwideJAV_results.sd.xlsx")
write.xlsx(MVNIslwideJAV_results.RE,"MVNIslwideJAV_results.RE.xlsx")
write.xlsx(MVNIslwideJAV_results.CI,"MVNIslwideJAV_results.CI.xlsx")
write.xlsx(MVNIslwideJAV_results.ICC,"MVNIslwideJAV_results.ICC.xlsx")

##############################################################################
#                      FCS-1L-DI-wide                                        #
##############################################################################

FCSslwide_results.est=matrix(NA,nrow=7,ncol=length(data))
FCSslwide_results.sd=matrix(NA,nrow=7,ncol=length(data))
FCSslwide_results.RE=matrix(NA,nrow=3,ncol=length(data))
FCSslwide_results.ICC=matrix(NA,nrow=2,ncol=length(data))
FCSslwide_results.CI=c()

for (i in 1:length(data)){
  
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
  
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7",
                                           "c_dep.3","c_dep.5","c_dep.7")]
  
  ##Set number of imputations and number of burn-in iterations
  M<-20
  
  ##create school dummy indicators
  simdataw$school=as.factor(simdataw$school)
  
  ##Perform imputations
  set.seed(2946)
  imp2=mice(data=simdataw[,!names(simdataw)%in%c("c_id","child")],m=M,method="norm",maxit=10)
  
  # Examine convergence
  #plot(imp2)
  
  ##Analysis
  mylist=list()
  
  for(m in 1:M){
    
    #extract the mth imputed dataset
    datw<-complete(imp2,m)
    
    
    #2. remove unwanted variables
    datw=datw[,names(datw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                "c_ses","c_nap1_z","id","school")]
    
    #3. rename depression variables for reshape
    colnames(datw)[colnames(datw)=="c_dep.2"] <- "prev_dep.3"
    colnames(datw)[colnames(datw)=="c_dep.4"] <- "prev_dep.5"
    colnames(datw)[colnames(datw)=="c_dep.6"] <- "prev_dep.7"
    
    #4. Reshape to long
    
    datL=reshape(datw,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),
                                    c("napscore_z.3","napscore_z.5","napscore_z.7")),idvar="id", 
                 v.names=c("prev_dep","napscore_z"), times=c(3,5,7),direction= "long")
    
    datL <- datL[order(datL$id),]    
    
    #5. save the dataset in a list
    mylist[[m]]= datL
    
    
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*time+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  FCSslwide_results.CI=cbind(FCSslwide_results.CI,CI)
  
  #store the estimates
  FCSslwide_results.est[,i]=MI_est$estimates[2:8,1]
  FCSslwide_results.sd[,i]=MI_est$estimates[2:8,2]
  FCSslwide_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  FCSslwide_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  FCSslwide_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  FCSslwide_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  FCSslwide_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
}

rownames(FCSslwide_results.est)=rows
colnames(FCSslwide_results.est)=c(seq(1:length(data)))

rownames(FCSslwide_results.sd)=rows
colnames(FCSslwide_results.sd)=c(seq(1:length(data)))

rownames(FCSslwide_results.RE)=c("level 3","level 2","level 1")
colnames(FCSslwide_results.RE)=c(seq(1:length(data)))

rownames(FCSslwide_results.ICC)=c("level 3","level 2")
colnames(FCSslwide_results.ICC)=c(seq(1:length(data)))

rownames(FCSslwide_results.CI)=rows
colnames(FCSslwide_results.CI)=c(rep(1:length(data), each=2))

##save the results
write.xlsx(FCSslwide_results.est,"FCSslwide_results.est.xlsx")
write.xlsx(FCSslwide_results.sd,"FCSslwide_results.sd.xlsx")
write.xlsx(FCSslwide_results.RE,"FCSslwide_results.RE.xlsx")
write.xlsx(FCSslwide_results.ICC,"FCSslwide_results.ICC.xlsx")
write.xlsx(FCSslwide_results.CI,"FCSslwide_results.CI.xlsx")

##############################################################################
#                      FCS-1L-DI-wide-Passive_c                              #
##############################################################################

FCSslwidepassivec_results.est=matrix(NA,nrow=7,ncol=length(data))
FCSslwidepassivec_results.sd=matrix(NA,nrow=7,ncol=length(data))
FCSslwidepassivec_results.RE=matrix(NA,nrow=3,ncol=length(data))
FCSslwidepassivec_results.ICC=matrix(NA,nrow=2,ncol=length(data))
FCSslwidepassivec_results.CI=c()

for (i in 1:length(data)){
 
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
 
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3","c_dep.5","c_dep.7")]
  
  ##Set number of imputations and number of burn-in iterations
  M<-20
  
  simdataw$school=as.factor(simdataw$school)
  
  ##make new variables for the passive product terms
  simdataw=data.frame(simdataw,Nap3.SES=NA,Nap5.SES=NA,Nap7.SES=NA,Dep2.Nap3=NA,Dep4.Nap5=NA,Dep6.Nap7=NA)
  
  ##specify the predictor matrix
  pred=mice::make.predictorMatrix(simdataw)
  pred[,]<-0
  codes<-c(rep(1,times=14))
  pred["c_dep.2",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap3.SES")]<-codes
  
  pred["c_dep.4",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap5.SES")]<-codes
  
  pred["c_dep.6",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",     
                   "napscore_z.5","p_sdq.6","napscore_z.7",
                   "Nap7.SES")]<-codes
  
  pred["c_ses",c("school","c_age","c_gender","c_nap1_z","p_sdq.2","c_dep.2",
                 "napscore_z.3","p_sdq.4","c_dep.4",    
                 "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                 "Dep2.Nap3","Dep4.Nap5","Dep6.Nap7")]<-c(codes,1,1)
  
  ##specify the methods for the missing variables
  meth=make.method(simdataw)
  meth[substr(row.names(pred),1,5)%in%c("c_dep")]="norm"
  meth["c_ses"]="norm"
  
  ##Derive interactions
  meth["Nap3.SES"] <- "~ I(napscore_z.3 * c_ses)"
  meth["Nap5.SES"] <- "~ I(napscore_z.5 * c_ses)"
  meth["Nap7.SES"] <- "~ I(napscore_z.7 * c_ses)"
  meth["Dep2.Nap3"] <- "~ I(c_dep.2 * napscore_z.3)"
  meth["Dep4.Nap5"] <- "~ I(c_dep.4 * napscore_z.5)"
  meth["Dep6.Nap7"] <- "~ I(c_dep.6 * napscore_z.7)"
  
  ##visit sequence #updates the derived variables that depend on the target variable 
  #(defines the order in which the variables are to be imputed)
  visit<-c("c_ses","Nap3.SES","Nap5.SES","Nap7.SES",
           "c_dep.2","Dep2.Nap3","c_dep.4","Dep4.Nap5","c_dep.6","Dep6.Nap7")
  
  ##Perform imputations
  set.seed(2946)
  imp2=mice(data=simdataw,m=M,predictorMatrix = pred,method=meth,visit=visit,maxit=10,allow.na=TRUE)
  
  #Examine convergence
  #plot(imp2)
  
  ##Analysis
  mylist=list()
  
  for(m in 1:M)
  {
    
    #extract the mth imputed dataset
    datw<-complete(imp2,m)
    
    #remove unwanted variables
    datw=datw[,names(datw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                "c_ses","c_nap1_z","c_id","school")]
    
    #rename depression variables for reshape
    colnames(datw)[colnames(datw)=="c_dep.2"] <- "prev_dep.3"
    colnames(datw)[colnames(datw)=="c_dep.4"] <- "prev_dep.5"
    colnames(datw)[colnames(datw)=="c_dep.6"] <- "prev_dep.7"
    
    #Reshape to long
    datL=reshape(datw,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                "napscore_z.5","napscore_z.7")),idvar="id", 
                 v.names=c("prev_dep","napscore_z"), times=c(3,5,7),direction= "long")
    
    datL <- datL[order(datL$id),]    
    
    #save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  FCSslwidepassivec_results.CI=cbind(FCSslwidepassivec_results.CI,CI)
  
  #store the estimates
  FCSslwidepassivec_results.est[,i]=MI_est$estimates[2:8,1]
  FCSslwidepassivec_results.sd[,i]=MI_est$estimates[2:8,2]
  FCSslwidepassivec_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  FCSslwidepassivec_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  FCSslwidepassivec_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  FCSslwidepassivec_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  FCSslwidepassivec_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
  
}
rownames(FCSslwidepassivec_results.est)=rows
colnames(FCSslwidepassivec_results.est)=c(seq(1:length(data)))

rownames(FCSslwidepassivec_results.sd)=rows
colnames(FCSslwidepassivec_results.sd)=c(seq(1:length(data)))

rownames(FCSslwidepassivec_results.RE)=c("level 3","level 2","level 1")
colnames(FCSslwidepassivec_results.RE)=c(seq(1:length(data)))

rownames(FCSslwidepassivec_results.ICC)=c("level 3","level 2")
colnames(FCSslwidepassivec_results.ICC)=c(seq(1:length(data)))

rownames(FCSslwidepassivec_results.CI)=rows
colnames(FCSslwidepassivec_results.CI)=c(rep(1:length(data), each=2))

##save the results
write.xlsx(FCSslwidepassivec_results.est,"FCSslwidepassivec_results.est.xlsx")
write.xlsx(FCSslwidepassivec_results.sd,"FCSslwidepassivec_results.sd.xlsx")
write.xlsx(FCSslwidepassivec_results.RE,"FCSslwidepassivec_results.RE.xlsx")
write.xlsx(FCSslwidepassivec_results.ICC,"FCSslwidepassivec_results.ICC.xlsx")
write.xlsx(FCSslwidepassivec_results.CI,"FCSslwidepassivec_results.CI.xlsx")

##############################################################################
#                      FCS-1L-DI-wide-Passive_all                            #
##############################################################################

FCSslwidepassiveall_results.est=matrix(NA,nrow=7,ncol=length(data))
FCSslwidepassiveall_results.sd=matrix(NA,nrow=7,ncol=length(data))
FCSslwidepassiveall_results.RE=matrix(NA,nrow=3,ncol=length(data))
FCSslwidepassiveall_results.ICC=matrix(NA,nrow=2,ncol=length(data))
FCSslwidepassiveall_results.CI=c()

for (i in 1:length(temp)){
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$child),]
  
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3","c_dep.5","c_dep.7")]
  
  ##Set number of imputations and number of burn-in iterations
  M<-20
  
  ##create school dummy indicators
  simdataw$school=as.factor(simdataw$school)
  
  ##make new variables for the product terms
  simdataw=data.frame(simdataw,Nap3.SES=NA,Nap5.SES=NA,Nap7.SES=NA,Dep2.Nap3=NA,
                      Dep4.Nap5=NA,Dep6.Nap7=NA)
  
  ##specify the predictor matrix
  pred=mice::make.predictorMatrix(simdataw)
  pred[,]<-0
  codes<-c(rep(1,times=16))
  pred["c_dep.2",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap3.SES","Nap5.SES","Nap7.SES")]<-codes
  
  pred["c_dep.4",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap3.SES","Nap5.SES","Nap7.SES")]<-codes
  
  pred["c_dep.6",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",     
                   "napscore_z.5","p_sdq.6","napscore_z.7",
                   "Nap3.SES","Nap5.SES","Nap7.SES")]<-codes
  
  pred["c_ses",c("school","c_age","c_gender","c_nap1_z","p_sdq.2","c_dep.2",
                 "napscore_z.3","p_sdq.4","c_dep.4",    
                 "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                 "Dep2.Nap3","Dep4.Nap5","Dep6.Nap7")]<-codes
  
  ##specify the methods for the missing variables
  meth=make.method(simdataw)
  meth[substr(row.names(pred),1,5)%in%c("c_dep")]="norm"
  meth["c_ses"]="norm"
  
  ##Derive interactions
  meth["Nap3.SES"] <- "~ I(napscore_z.3 * c_ses)"
  meth["Nap5.SES"] <- "~ I(napscore_z.5 * c_ses)"
  meth["Nap7.SES"] <- "~ I(napscore_z.7 * c_ses)"
  meth["Dep2.Nap3"] <- "~ I(c_dep.2 * napscore_z.3)"
  meth["Dep4.Nap5"] <- "~ I(c_dep.4 * napscore_z.5)"
  meth["Dep6.Nap7"] <- "~ I(c_dep.6 * napscore_z.7)"
  
  ##visit sequence #updates the derived variables that depend on the target variable 
  #(defines the order in which the variables are to be imputed)
  visit<-c("c_ses","Nap3.SES","Nap5.SES","Nap7.SES",
           "c_dep.2","Dep2.Nap3","c_dep.4","Dep4.Nap5","c_dep.6","Dep6.Nap7")
  
  ##Perform imputations
  set.seed(2946)
  imp2=mice(data=simdataw,m=M,predictorMatrix = pred,method=meth,visit=visit,maxit=10,allow.na=TRUE)
  
  #Examine convergence
  #plot(imp2)
  
  ##Analysis
  mylist=list()
  
  for(m in 1:M)
  {
    
    #extract the mth imputed dataset
    datw<-complete(imp2,m)
    
    #remove unwanted variables
    datw=datw[,names(datw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                "c_ses","c_nap1_z","c_id","school")]
  
    #rename depression variables for reshape
    colnames(datw)[colnames(datw)=="c_dep.2"] <- "prev_dep.3"
    colnames(datw)[colnames(datw)=="c_dep.4"] <- "prev_dep.5"
    colnames(datw)[colnames(datw)=="c_dep.6"] <- "prev_dep.7"
    
    #Reshape to long
    datL=reshape(datw,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                "napscore_z.5","napscore_z.7")),idvar="id", 
                 v.names=c("prev_dep","napscore_z"), times=c(3,5,7),direction= "long")
    
    datL <- datL[order(datL$id),]    
    
    #save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  FCSslwidepassiveall.CI=cbind(FCSslwidepassiveall.CI,CI)
  
  #store the estimates
  FCSslwidepassiveall.est[,i]=MI_est$estimates[2:8,1]
  FCSslwidepassiveall.sd[,i]=MI_est$estimates[2:8,2]
  FCSslwidepassiveall.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  FCSslwidepassiveall.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  FCSslwidepassiveall.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  FCSslwidepassiveall.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  FCSslwidepassiveall.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
}

rownames(FCSslwidepassiveall_results.est)=rows
colnames(FCSslwidepassiveall_results.est)=c(seq(1:length(data)))

rownames(FCSslwidepassiveall_results.sd)=rows
colnames(FCSslwidepassiveall_results.sd)=c(seq(1:length(data)))

rownames(FCSslwidepassiveall_results.RE)=c("level 3","level 2","level 1")
colnames(FCSslwidepassiveall_results.RE)=c(seq(1:length(data)))

rownames(FCSslwidepassiveall_results.ICC)=c("level 3","level 2")
colnames(FCSslwidepassiveall_results.ICC)=c(seq(1:length(data)))

rownames(FCSslwidepassiveall_results.CI)=rows
colnames(FCSslwidepassiveall_results.CI)=c(rep(1:length(data), each=2))

##save the results
write.xlsx(FCSslwidepassiveall_results.est,"FCSslwidepassiveall_results.est.xlsx")
write.xlsx(FCSslwidepassiveall_results.sd,"FCSslwidepassiveall_results.sd.xlsx")
write.xlsx(FCSslwidepassiveall_results.RE,"FCSslwidepassiveall_results.RE.xlsx")
write.xlsx(FCSslwidepassiveall_results.ICC,"FCSslwidepassiveall_results.ICC.xlsx")
write.xlsx(FCSslwidepassiveall_results.CI,"FCSslwidepassiveall_results.CI.xlsx")

##############################################################################
#                      JM-2L-wide-JAV                                        #
##############################################################################

MVNImlwideJAV_results.est=matrix(NA,nrow=7,ncol=length(data))
MVNImlwideJAV_results.sd=matrix(NA,nrow=7,ncol=length(data))
MVNImlwideJAV_results.RE=matrix(NA,nrow=3,ncol=length(data))
MVNImlwideJAV_results.ICC=matrix(NA,nrow=2,ncol=length(data))
MVNImlwideJAV_results.CI=c()

for (i in 1:length(data)){
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
  ##Reshape to wide format
  
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3","c_dep.5","c_dep.7")]
  
  ##generate the interaction terms
  simdataw$c_dep.ses.2=simdataw$c_dep.2*simdataw$c_ses
  simdataw$c_dep.ses.4=simdataw$c_dep.4*simdataw$c_ses
  simdataw$c_dep.ses.6=simdataw$c_dep.6*simdataw$c_ses
  
  ##Set number of imputations and number of burn-in iterations
  M<-20
  nburn=1000
  NB=100
  
  ##create a dataframe with variables to be imputed
  myvars <- names(simdataw) %in% c("c_dep.2","c_dep.4","c_dep.6","c_ses","c_dep.ses.2",
                                   "c_dep.ses.4","c_dep.ses.6") 
  dataimp=simdataw[myvars]
  
  ##create a dataframe with complete variables
  datacomp= cbind(Intercept=rep(1,nrow(simdataw)),
                  simdataw[,names(simdataw)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                                 "c_age","c_gender",
                                                 "c_nap1_z","p_sdq.2","p_sdq.4","p_sdq.6" )])
  
  # Create a data frame  with column of 1's for random intercept
  datcompRE<-cbind(Intercept=rep(1,nrow(simdataw)))
  
  
  ##perform imputations without the random effects
  set.seed(2946)
  imp5<-jomo1rancon(Y=dataimp, X=datacomp, Z=datcompRE, clus=simdataw$school, nimp=M, nburn=nburn,nbetween = NB )
 
   #check convergence
  # impCheck<-jomo1rancon.MCMCchain(Y=dataimp, X=datacomp, Z=datcompRE, clus=simdataw$school,nburn=nburn)
  # plot(c(1:nburn),impCheck$collectbeta[1,1,1:nburn],type="l")
  
  
  mylist=list()
  
  for(m in 1:M) {
    
    #extract the ith imputed data set
    imputed<-imp5[imp5$Imputation==m,]
    imputed=imputed[,!names(imputed)%in%c("p_sdq.2","p_sdq.4","p_sdq.6", "Intercept","Imputation")]
    
    #rename depression for reshape
    colnames(imputed)[colnames(imputed)=="c_dep.2"] <- "prev_dep.3"
    colnames(imputed)[colnames(imputed)=="c_dep.4"] <- "prev_dep.5"
    colnames(imputed)[colnames(imputed)=="c_dep.6"] <- "prev_dep.7"
    
    colnames(imputed)[colnames(imputed)=="c_dep.ses.2"] <- "inter.3"
    colnames(imputed)[colnames(imputed)=="c_dep.ses.4"] <- "inter.5"
    colnames(imputed)[colnames(imputed)=="c_dep.ses.6"] <- "inter.7"
    
    #reshape to long format
    datL=reshape(imputed,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                   "napscore_z.5","napscore_z.7"),
                                       c("inter.3","inter.5","inter.7")),idvar="id", 
                 v.names=c("prev_dep","napscore_z","inter"), times=c(3,5,7),direction= "long")
    
    
    datL <- datL[order(datL$clus,datL$id),]
    
    #save the dataset in a list
    mylist[[m]]= datL
  }
  
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+c_age+
                                             c_gender+c_nap1_z+c_ses+inter
                                           +(1|clus/id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #store the estimates
  MVNImlwideJAV_results.est[,i]=MI_est$estimates[2:8,1]
  MVNImlwideJAV_results.sd[,i]=MI_est$estimates[2:8,2]
  MVNImlwideJAV_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  MVNImlwideJAV_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  MVNImlwideJAV_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  MVNImlwideJAV_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  MVNImlwideJAV_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  MVNImlwideJAV_results.CI=cbind(MVNImlwideJAV_results.CI,CI)  
}

rownames(MVNImlwideJAV_results.est)=rows
colnames(MVNImlwideJAV_results.est)=c(seq(1:length(temp)))

rownames(MVNImlwideJAV_results.sd)=rows
colnames(MVNImlwideJAV_results.sd)=c(seq(1:length(temp)))

rownames(MVNImlwideJAV_results.RE)=c("level 3","level 2","level 1")
colnames(MVNImlwideJAV_results.RE)=c(seq(1:length(temp)))

rownames(MVNImlwideJAV_results.ICC)=c("level 3","level 2")
colnames(MVNImlwideJAV_results.ICC)=c(seq(1:length(temp)))

rownames(MVNImlwideJAV_results.CI)=rows
colnames(MVNImlwideJAV_results.CI)=c(rep(1:length(temp),each=2))

##save the results
write.xlsx(MVNImlwideJAV_results.est,"MVNImlwideJAV_results.est.xlsx")
write.xlsx(MVNImlwideJAV_results.sd,"MVNImlwideJAV_results.sd.xlsx")
write.xlsx(MVNImlwideJAV_results.RE,"MVNImlwideJAV_results.RE.xlsx")
write.xlsx(MVNImlwideJAV_results.CI,"MVNImlwideJAV_results.CI.xlsx")
write.xlsx(MVNImlwideJAV_results.ICC,"MVNImlwideJAV_results.ICC.xlsx")

##############################################################################
#                      FCS-2L-wide-Passive_c                                 #
##############################################################################

FCSmlwidepassivec_results.est=matrix(NA,nrow=7,ncol=length(data))
FCSmlwidepassivec_results.sd=matrix(NA,nrow=7,ncol=length(data))
FCSmlwidepassivec_results.RE=matrix(NA,nrow=3,ncol=length(data))
FCSmlwidepassivec_results.ICC=matrix(NA,nrow=2,ncol=length(data))
FCSmlwidepassivec_results.CI=c()

for (i in 1:length(data)){
  
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$child),]
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3","c_dep.5","c_dep.7","child")]
  
  ##make new variables for the product terms
  simdataw=data.frame(simdataw,Nap3.SES=NA,Nap5.SES=NA,Nap7.SES=NA,Dep2.Nap3=NA,
                      Dep4.Nap5=NA,Dep6.Nap7=NA)
  
  ##specify the predictor matrix
  pred=mice::make.predictorMatrix(simdataw)
  pred[,]<-0
  
  codes<-c(-2,rep(1,times=13))
  pred["c_dep.2",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap3.SES")]<-codes
  
  pred["c_dep.4",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap5.SES")]<-codes
  
  pred["c_dep.6",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",     
                   "napscore_z.5","p_sdq.6","napscore_z.7",
                   "Nap7.SES")]<-codes
  
  pred["c_ses",c("school","c_age","c_gender","c_nap1_z","p_sdq.2","c_dep.2",
                 "napscore_z.3","p_sdq.4","c_dep.4",    
                 "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                 "Dep2.Nap3","Dep4.Nap5","Dep6.Nap7")]<-c(codes,1,1)
  
  ##specify the methods for the missing variables
  meth=make.method(simdataw)
  meth[substr(row.names(pred),1,5)%in%c("c_dep")]="2l.pan"
  meth["c_ses"]="2l.pan"
  
  ##Derive interactions
  meth["Nap3.SES"] <- "~ I(napscore_z.3 * c_ses)"
  meth["Nap5.SES"] <- "~ I(napscore_z.5 * c_ses)"
  meth["Nap7.SES"] <- "~ I(napscore_z.7 * c_ses)"
  meth["Dep2.Nap3"] <- "~ I(c_dep.2 * napscore_z.3)"
  meth["Dep4.Nap5"] <- "~ I(c_dep.4 * napscore_z.5)"
  meth["Dep6.Nap7"] <- "~ I(c_dep.6 * napscore_z.7)"
  
  ##visit sequence #updates the derived variables that depend on the target variable 
  #(defines the order in which the variables are to be imputed)
  visit<-c("c_ses","Nap3.SES","Nap5.SES","Nap7.SES",
           "c_dep.2","Dep2.Nap3","c_dep.4","Dep4.Nap5","c_dep.6","Dep6.Nap7")
  
  ##set the number of imputations and iterations
  M=20
  
  #perform the imputations
  imp4<-mice(simdataw,m=M, maxit=10,predictorMatrix=pred, method=meth,seed=2384,visit=visit,allow.na=TRUE)
  
  #Examine convergence
  #plot(imp4)
  
  ##analysis
  mylist=list()
  
  for(m in 1:M)
  {
    #extract the ith imputed data set
    imputed=complete(imp4,m)
    imputed=imputed[,names(imputed)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                         "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                         "c_ses","c_nap1_z","c_id","school")]
    
    colnames(imputed)[colnames(imputed)=="c_dep.2"] <- "prev_dep.3"
    colnames(imputed)[colnames(imputed)=="c_dep.4"] <- "prev_dep.5"
    colnames(imputed)[colnames(imputed)=="c_dep.6"] <- "prev_dep.7"
    
    #reshape to long format
    datL=reshape(imputed,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                   "napscore_z.5","napscore_z.7")),idvar="c_id", 
                 v.names=c("prev_dep","napscore_z"), times=c(3,5,7),direction= "long")
    
    datL <- datL[order(datL$school,datL$c_id),]
    
    #save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #store the estimates
  FCSmlwidepassivec_results.est[,i]=MI_est$estimates[2:8,1]
  FCSmlwidepassivec_results.sd[,i]=MI_est$estimates[2:8,2]
  FCSmlwidepassivec_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  FCSmlwidepassivec_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  FCSmlwidepassivec_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  FCSmlwidepassivec_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  FCSmlwidepassivec_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  FCSmlwidepassivec_results.CI=cbind(FCSmlwidepassivec_results.CI,CI)

}

rownames(FCSmlwidepassivec_results.est)=rows
colnames(FCSmlwidepassivec_results.est)=c(seq(1:length(data)))

rownames(FCSmlwidepassivec_results.sd)=rows
colnames(FCSmlwidepassivec_results.sd)=c(seq(1:length(data)))

rownames(FCSmlwidepassivec_results.RE)=c("level 3","level 2","level 1")
colnames(FCSmlwidepassivec_results.RE)=c(seq(1:length(data)))

rownames(FCSmlwidepassivec_results.ICC)=c("level 3","level 2")
colnames(FCSmlwidepassivec_results.ICC)=c(seq(1:length(data)))

rownames(FCSmlwidepassivec_results.CI)=rows
colnames(FCSmlwidepassivec_results.CI)=c(rep(1:length(data),each=2))

##save the results
write.xlsx(FCSmlwidepassivec_results.est,"FCSmlwidepassivec_results.est.xlsx")
write.xlsx(FCSmlwidepassivec_results.sd,"FCSmlwidepassivec_results.sd.xlsx")
write.xlsx(FCSmlwidepassivec_results.RE,"FCSmlwidepassivec_results.RE.xlsx")
write.xlsx(FCSmlwidepassivec_results.CI,"FCSmlwidepassivec_results.CI.xlsx")
write.xlsx(FCSmlwidepassivec_results.ICC,"FCSmlwidepassivec_results.ICC.xlsx")

##############################################################################
#                      FCS-2L-wide-Passive_all                               #
##############################################################################

FCSmlwidepassiveall_results.est=matrix(NA,nrow=7,ncol=length(data))
FCSmlwidepassiveall_results.sd=matrix(NA,nrow=7,ncol=length(data))
FCSmlwidepassiveall_results.RE=matrix(NA,nrow=3,ncol=length(data))
FCSmlwidepassiveall_results.ICC=matrix(NA,nrow=2,ncol=length(data))
FCSmlwidepassiveall_results.CI=c()

for (i in 1:length(data)){
  
  simdataL= data[[i]]
  simdataL <- simdataL[order(simdataL$school,simdataL$child),]
  
  ##Reshape to wide format
  simdataw=reshape(simdataL,v.names=c("napscore_z","p_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  ##remove unwanted variables
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","napscore_z.4","napscore_z.6",
                                           "p_sdq.3","p_sdq.5","p_sdq.7","c_dep.3","c_dep.5","c_dep.7","child")]
  
  ######Multilevel FCS with JAV for naplan scores###################
  
  ##make new variables for the product terms
  simdataw=data.frame(simdataw,Nap3.SES=NA,Nap5.SES=NA,Nap7.SES=NA,Dep2.Nap3=NA,
                      Dep4.Nap5=NA,Dep6.Nap7=NA)
  
  ##specify the predictor matrix
  pred=mice::make.predictorMatrix(simdataw)
  pred[,]<-0
  
  codes<-c(-2,rep(1,times=15))
  pred["c_dep.2",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap3.SES","Nap5.SES","Nap7.SES")]<-codes
  
  pred["c_dep.4",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4",    
                   "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                   "Nap3.SES","Nap5.SES","Nap7.SES")]<-codes
  
  pred["c_dep.6",c("school","c_age","c_gender","c_ses","c_nap1_z","p_sdq.2","c_dep.2",
                   "napscore_z.3","p_sdq.4","c_dep.4",     
                   "napscore_z.5","p_sdq.6","napscore_z.7",
                   "Nap3.SES","Nap5.SES","Nap7.SES")]<-codes
  
  pred["c_ses",c("school","c_age","c_gender","c_nap1_z","p_sdq.2","c_dep.2",
                 "napscore_z.3","p_sdq.4","c_dep.4",    
                 "napscore_z.5","p_sdq.6","c_dep.6","napscore_z.7",
                 "Dep2.Nap3","Dep4.Nap5","Dep6.Nap7")]<-codes
  
  ##specify the methods for the missing variables
  meth=make.method(simdataw)
  meth[substr(row.names(pred),1,5)%in%c("c_dep")]="2l.pan"
  meth["c_ses"]="2l.pan"
  
  ##Derive interactions
  meth["Nap3.SES"] <- "~ I(napscore_z.3 * c_ses)"
  meth["Nap5.SES"] <- "~ I(napscore_z.5 * c_ses)"
  meth["Nap7.SES"] <- "~ I(napscore_z.7 * c_ses)"
  meth["Dep2.Nap3"] <- "~ I(c_dep.2 * napscore_z.3)"
  meth["Dep4.Nap5"] <- "~ I(c_dep.4 * napscore_z.5)"
  meth["Dep6.Nap7"] <- "~ I(c_dep.6 * napscore_z.7)"
  
  ##visit sequence #updates the derived variables that depend on the target variable 
  #(defines the order in which the variables are to be imputed)
  visit<-c("c_ses","Nap3.SES","Nap5.SES","Nap7.SES",
           "c_dep.2","Dep2.Nap3","c_dep.4","Dep4.Nap5","c_dep.6","Dep6.Nap7")
  
  ##set the number of imputations and iterations
  M=20
  
  #perform the imputations
  imp4<-mice(simdataw,m=M, maxit=10,predictorMatrix=pred, method=meth,seed=2384,visit=visit,allow.na=TRUE)
  
  #Examine convergence
  #plot(imp4)
  
  ##analysis
  mylist=list()
  
  for(m in 1:M)
  {
    #extract the ith imputed data set
    imputed=complete(imp4,m)
    imputed=imputed[,names(imputed)%in%c("napscore_z.3","napscore_z.5","napscore_z.7",
                                         "c_age","c_gender","c_dep.2","c_dep.4","c_dep.6",
                                         "c_ses","c_nap1_z","c_id","school")]
    
    colnames(imputed)[colnames(imputed)=="c_dep.2"] <- "prev_dep.3"
    colnames(imputed)[colnames(imputed)=="c_dep.4"] <- "prev_dep.5"
    colnames(imputed)[colnames(imputed)=="c_dep.6"] <- "prev_dep.7"
   
     #reshape to long format
    datL=reshape(imputed,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                   "napscore_z.5","napscore_z.7")),idvar="c_id", 
                 v.names=c("prev_dep","napscore_z"), times=c(3,5,7),direction= "long")
    
    datL <- datL[order(datL$school,datL$c_id),]
  
    #save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #store the estimates
  FCSmlwidepassiveall_results.est[,i]=MI_est$estimates[2:8,1]
  FCSmlwidepassiveall_results.sd[,i]=MI_est$estimates[2:8,2]
  FCSmlwidepassiveall_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  FCSmlwidepassiveall_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  FCSmlwidepassiveall_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  FCSmlwidepassiveall_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  FCSmlwidepassiveall_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  FCSmlwidepassiveall_results.CI=cbind(FCSmlwidepassiveall_results.CI,CI)
}

rownames(FCSmlwidepassiveall_results.est)=rows
colnames(FCSmlwidepassiveall_results.est)=c(seq(1:length(data)))

rownames(FCSmlwidepassiveall_results.sd)=rows
colnames(FCSmlwidepassiveall_results.sd)=c(seq(1:length(data)))

rownames(FCSmlwidepassiveall_results.RE)=c("level 3","level 2","level 1")
colnames(FCSmlwidepassiveall_results.RE)=c(seq(1:length(data)))

rownames(FCSmlwidepassiveall_results.ICC)=c("level 3","level 2")
colnames(FCSmlwidepassiveall_results.ICC)=c(seq(1:length(data)))

rownames(FCSmlwidepassiveall_results.CI)=rows
colnames(FCSmlwidepassiveall_results.CI)=c(rep(1:length(data),each=2))

##save the results
write.xlsx(FCSmlwidepassiveall_results.est,"FCSmlwidepassiveall_results.est.xlsx")
write.xlsx(FCSmlwidepassiveall_results.sd,"FCSmlwidepassiveall_results.sd.xlsx")
write.xlsx(FCSmlwidepassiveall_results.RE,"FCSmlwidepassiveall_results.RE.xlsx")
write.xlsx(FCSmlwidepassiveall_results.CI,"FCSmlwidepassiveall_results.CI.xlsx")
write.xlsx(FCSmlwidepassiveall_results.ICC,"FCSmlwidepassiveall_results.ICC.xlsx")

##############################################################################
#                      SMC-JM-2L-DI                                          #
##############################################################################

SMCJMDI_results.est=matrix(NA,nrow=7,ncol=length(data))
SMCJMDI_results.sd=matrix(NA,nrow=7,ncol=length(data))
SMCJMDI_results.RE=matrix(NA,nrow=3,ncol=length(data))
SMCJMDI_results.ICC=matrix(NA,nrow=2,ncol=length(data))
SMCJMDI_results.CI=c()

for (i in 1:length(data)){
  
  simdataL= data[[i]]
  simdataL$c_id=as.numeric(simdataL$c_id)
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
  
  ##rearrange the data set for imputations
  
  #1. SDQ variable
  simdataL<- slide(simdataL, Var = "p_sdq", GroupVar = "c_id",
                   slideBy = -1)
  simdataL=simdataL[,!names(simdataL)%in%c("p_sdq")]
  colnames(simdataL)[colnames(simdataL)=="p_sdq-1"] <- "prev_sdq"
  
  #2. Depression values at waves 2, 4 and 6
  simdataw=reshape(simdataL,v.names=c("napscore_z","prev_sdq","c_dep"),idvar="c_id", 
                   timevar="wave", direction= "wide")
  
  colnames(simdataw)[colnames(simdataw)=="c_dep.2"] <- "prev_dep.3"
  colnames(simdataw)[colnames(simdataw)=="c_dep.4"] <- "prev_dep.5"
  colnames(simdataw)[colnames(simdataw)=="c_dep.6"] <- "prev_dep.7"
  
  simdataw=simdataw[,!names(simdataw)%in%c("napscore_z.2","prev_sdq.2","napscore_z.4","prev_sdq.4",
                                           "napscore_z.6","prev_sdq.6","c_dep.3","c_dep.5","c_dep.7")]
  
  ##reshape back to long
  simdataL=reshape(simdataw,varying =list(c("prev_dep.3","prev_dep.5","prev_dep.7"),c("napscore_z.3",
                                                                                      "napscore_z.5",  "napscore_z.7"),
                                          c("prev_sdq.3","prev_sdq.5","prev_sdq.7")),idvar="c_id", 
                   v.names=c("prev_dep","napscore_z","prev_sdq"), times=c(3,5,7),direction= "long")
  
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
  simdataL$c_gender=as.factor(simdataL$c_gender)
  
  ##Set number of imputations and number of burn-in iterations
  M<-20
  nburn=500
  NB=10
  
  ##create school dummy indicators
  school_DI=data.frame(model.matrix(simdataL$c_id~as.factor(simdataL$school)-1,
                                    simdataL))
  names(school_DI)[1:ncol(school_DI)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), 
                                                       "schoo_Ind",40))
  school_DI=school_DI[,1:39]
  
  
  ##define a dataframe with all variables required
  data_jomo=cbind(simdataL[!names(simdataL)%in%c("child","school")],school_DI)
  
  #specify the levels of each variable in the imputation model 
  mylevel<-c(2,2,2,2,2,1,1,1,1,rep(2,times=39))
  
  # formula of the substantive lmer model
  formula<-as.formula(napscore_z~prev_dep+time+prev_dep*c_ses+c_age+c_gender+c_ses+c_nap1_z+schoo_Ind_1+
                        schoo_Ind_2+schoo_Ind_3+schoo_Ind_4+schoo_Ind_5+schoo_Ind_6+schoo_Ind_7+schoo_Ind_8+
                        schoo_Ind_9+schoo_Ind_10+schoo_Ind_11+schoo_Ind_12+schoo_Ind_13+schoo_Ind_14+schoo_Ind_15+
                        schoo_Ind_16+schoo_Ind_17+schoo_Ind_18+schoo_Ind_19+schoo_Ind_20+schoo_Ind_21+schoo_Ind_22+
                        schoo_Ind_23+schoo_Ind_24+schoo_Ind_25+schoo_Ind_26+schoo_Ind_27+schoo_Ind_28+schoo_Ind_29+
                        schoo_Ind_30+schoo_Ind_31+schoo_Ind_32+schoo_Ind_33+schoo_Ind_34+schoo_Ind_35+schoo_Ind_36+schoo_Ind_37+
                        schoo_Ind_38+schoo_Ind_39+prev_sdq+(1|c_id))
  
  ##running the imputations
  set.seed(2946)
  imp5<-jomo.lmer(formula,data_jomo,level=mylevel,nimp=M,nburn=nburn,nbetween = NB)
  
  #Examine convergence
  # impcheck=jomo.lmer.MCMCchain(formula,data_jomo,level=mylevel,nburn=nburn)
  # plot(c(1:nburn),impcheck$collectbeta[1,1,1:nburn],type="l")
  
  mylist=list()
  
  for(m in 1:M) {
    
    #extract the ith imputed data set and attach the school variable
    imputed<-imp5[imp5$Imputation==m,]
    imputed=imputed[,names(imputed)%in%c("prev_dep","c_age", "c_gender",
                                         "c_ses","c_nap1_z","time",
                                         "napscore_z","clus")]
    imputed$school=simdataL$school
    #imputed$clus=as.integer(imputed$clus)
    
    #save the dataset in a list for analysis 
    mylist[[m]]= imputed
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+time+prev_dep*c_ses+c_age+
                                             c_gender+c_nap1_z+c_ses
                                           +(1|school/clus), data = d)} )
  
  #combine the estimates 
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  SMCJMDI_results.CI=cbind(SMCJMDI_results.CI,CI)
  
  #store the estimates
  SMCJMDI_results.est[,i]=MI_est$estimates[2:8,1]
  SMCJMDI_results.sd[,i]=MI_est$estimates[2:8,2]
  SMCJMDI_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  SMCJMDI_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  SMCJMDI_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  SMCJMDI_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  SMCJMDI_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
}

rownames(SMCJMDI_results.est)=rows
colnames(SMCJMDI_results.est)=c(seq(1:length(temp)))

rownames(SMCJMDI_results.sd)=rows
colnames(SMCJMDI_results.sd)=c(seq(1:length(temp)))

rownames(SMCJMDI_results.RE)=c("level 3","level 2","level 1")
colnames(SMCJMDI_results.RE)=c(seq(1:length(temp)))

rownames(SMCJMDI_results.ICC)=c("level 3","level 2")
colnames(SMCJMDI_results.ICC)=c(seq(1:length(temp)))

rownames(SMCJMDI_results.CI)=rows
colnames(SMCJMDI_results.CI)=c(rep(1:length(temp), each=2))

##save the results
write.xlsx(SMCJMDI_results.est,"SMCJMDI_results.est.xlsx")
write.xlsx(SMCJMDI_results.sd,"SMCJMDI_results.sd.xlsx")
write.xlsx(SMCJMDI_results.RE,"SMCJMDI_results.RE.xlsx")
write.xlsx(SMCJMDI_results.CI,"SMCJMDI_results.CI.xlsx")
write.xlsx(SMCJMDI_results.ICC,"SMCJMDI_results.ICC.xlsx")

##############################################################################
#                      SMC-SM-2L-DI                                          #
##############################################################################

SMCSMDI_results.est=matrix(NA,nrow=7,ncol=length(data))
SMCSMDI_results.sd=matrix(NA,nrow=7,ncol=length(data))
SMCSMDI_results.RE=matrix(NA,nrow=3,ncol=length(data))
SMCSMDI_results.ICC=matrix(NA,nrow=2,ncol=length(data))
SMCSMDI_results.CI=c()

for (i in 1:length(temp)){
  
  simdataL= data[[i]]
  simdataL$c_id=as.numeric(simdataL$c_id)
  simdataL <- simdataL[order(simdataL$school,simdataL$c_id),]
  
  #create the previous wave deperssion symptoms
  simdataL<- slide(simdataL, Var = "c_dep", GroupVar = "c_id",
                   slideBy = -1)
  colnames(simdataL)[colnames(simdataL)=="c_dep-1"] <- "prev_dep"
  simdataL$c_dep=NULL
  
  #create previous wave SDQ variable
  simdataL<- slide(simdataL, Var = "p_sdq", GroupVar = "c_id",
                   slideBy = -1)
  colnames(simdataL)[colnames(simdataL)=="p_sdq-1"] <- "prev_sdq"
  simdataL$p_sdq=NULL
  
  #remove unwanted waves and variables
  simdataL=subset(simdataL, wave!= 2 & wave!=4 & wave!=6)
  simdataL=simdataL[,!names(simdataL)%in%c("child")]
  
  #recode categorical and binary varaible with numeric codes 0,1,2,...
  simdataL$c_gender=ifelse(simdataL$c_gender=="male",1,0)
  
  simdataL$school=simdataL$school-1
  
  #model specification
  iter <-2000 ; burnin <- 1000 ##burn in for 1000 iteration and one imputed dataset is
  Nimp <- 20                               ##saved every 100th iteration=20imps
  
  #substantive model formula
  Y_formula <- napscore_z ~prev_dep+wave+prev_dep*c_ses+c_age+c_ses+c_gender+c_nap1_z+prev_sdq+school+(1|c_id)
  
  #model for the outcome/Y variable
  dep <- list("model"="mlreg", "formula"=Y_formula )
  
  ##covariate models
  #1.prev_dep
  X1_formula=prev_dep~wave+c_age+c_ses+c_gender+c_nap1_z+prev_sdq+school+(1|c_id)
  
  ind_x <- list( "model"="mlreg", "formula"=X1_formula,
                 sampling_level="c_id" )
  #2. SES
  X2_formula=c_ses~c_age+c_gender+c_nap1_z+school
  
  ind_x2<-list( "model"="linreg", "formula"=X2_formula,
                variable_level="c_id")
  
  #ordering of the conditional covariate models
  ind <- list(c_ses=ind_x2,prev_dep=ind_x)
  
  #estimate model
  mod1 <- mdmb::frm_fb(simdataL, dep, ind,iter=iter,burnin=burnin, Nimp =20,aggregation=TRUE)
  
  #convert output into list of imputed datasets
  datlist <- mdmb::frm2datlist(mod1)
  
  mylist=list()
  
  for(j in 1:Nimp){
    
    #extract the mth imputed dataset
    datL<-datlist[[j]]
    
    #reattach the school cluster indicator 
    datL$school=simdataL$school
    
    datL=datL[,names(datL)%in% c("c_id","wave","napscore_z","c_age","c_gender","c_nap1_z","prev_dep",
                                 "prev_sdq","c_ses","school")]
    
    #save the dataset in a list
    mylist[[j]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer( napscore_z~prev_dep+wave+prev_dep*c_ses+c_age+
                                             as.factor(c_gender)+c_nap1_z+c_ses
                                           +(1|school/c_id), data = d)} )
  
  MI_est=testEstimates(mods, var.comp=TRUE,df.com=NULL)
  
  #Compute CIs
  CI=matrix(NA,7,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:8,1])
  CI[,2]=as.vector(Conf[2:8,2])
  SMCSMDI_results.CI=cbind(SMCSMDI_results.CI,CI)
  
  
  #store the estimates
  SMCSMDI_results.est[,i]=MI_est$estimates[2:8,1]
  SMCSMDI_results.sd[,i]=MI_est$estimates[2:8,2]
  SMCSMDI_results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  SMCSMDI_results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  SMCSMDI_results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  SMCSMDI_results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  SMCSMDI_results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
}

rownames(SMCSMDI_results.est)=rows
colnames(SMCSMDI_results.est)=c(seq(1:length(temp)))

rownames(SMCSMDI_results.sd)=rows
colnames(SMCSMDI_results.sd)=c(seq(1:length(temp)))

rownames(SMCSMDI_results.RE)=c("level 3","level 2","level 1")
colnames(SMCSMDI_results.RE)=c(seq(1:length(temp)))

rownames(SMCSMDI_results.ICC)=c("level 3","level 2")
colnames(SMCSMDI_results.ICC)=c(seq(1:length(temp)))

rownames(SMCSMDI_results.CI)=rows
colnames(SMCSMDI_results.CI)=c(rep(1:length(temp), each=2))

##save the results
write.xlsx(SMCSMDI_results.est,"SMCSMDI_results.est.xlsx")
write.xlsx(SMCSMDI_results.sd,"SMCSMDI_results.sd.xlsx")
write.xlsx(SMCSMDI_results.RE,"SMCSMDI_results.RE.xlsx")
write.xlsx(SMCSMDI_results.CI,"SMCSMDI_results.CI.xlsx")
write.xlsx(SMCSMDI_results.ICC,"SMCSMDI_results.ICC.xlsx")