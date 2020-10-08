*************************************************************************************
*     Simulation of complete and missing data for the simulation scenarios          *
*     evaluated in the paper:                                                       *
*  "Evaluation of approaches for accommodating interactions and non-linear terms    *
*   in multiple imputation of incomplete three-level"                                *
*             Rushani Wijesuriya                                                    *
*            16th of September 2020                                                 *
*************************************************************************************


**NOTE: This script is only for data generation, all MI approaches evaluated in the 
*paper were implemented using R 

clear all
set more off
version 15.1
	
cap prog drop simul
prog define simul

	syntax, seed(integer) simno(integer) schoolclus(integer) clussize(integer) model(string) mdm(string)
	set seed `seed'			 				 								 
		  	  
 forvalues i=1(1)`simno'{

clear

**STEP (i):

*Generate the schools
set obs `schoolclus'
gen school= _n

*Simulate values for the random effects at school level (Level 3)
gen u_0i=rnormal(0,0.1)
gen b_0i=rnormal(0,sqrt(0.05))
gen v_0i=rnormal(0,0.6)

* Expand the schools to include students in a school with a total of 1200
expand `clussize'

bysort school: generate child = _n

**STEP (ii):
*Simulate values for  child's age
gen c_age=runiform(7,10)

**STEP (iii):
*Simulate values for child's sex (female=0,male=1)
gen uran= runiform(0,1)
gen c_gender=0 if uran<=0.5 
replace c_gender=1 if uran>0.5
drop uran 

**STEP (iv):
* Simulate values of child SES
gen c_ses=rnormal(0,1)

**STEP (v):
* Simulate values for wave 1 Naplan scores
gen e_nap1=rnormal(0,1)
gen c_nap1_z= -0.74+ 0.23*c_gender+0.07*c_age+0.22*c_ses+e_nap1
drop e_nap1

*********************************************************
*Simulating time dependent variables


*Simulate the values for the random effects at individual level (Level 2)
gen u_0ij=rnormal(0,0.9)
gen b_0ij=rnormal(0,sqrt(0.45))
gen v_0ij=rnormal(0,4.1)

*Expand the individuals to include 5 repeated observations (waves 2-6)
expand 6

bysort school child: generate wave = _n+1

*Generate a child ID
gen c_id=string(school)+"_"+string(child)

*STEP (vi):
* Simulate values for depression (mean centered) at waves 2-6
gen e_c_dep=rnormal(0,1.5)
gen c_dep= -0.7+0.10*c_age +(-0.46)*c_gender+ (-0.01)*c_nap1_z+(-0.22)*c_ses + 0.02*wave+ ///
u_0i+ u_0ij+e_c_dep if wave==2 |wave==4|wave==6
drop e_c_dep

*Create a variable for previous wave depression
bysort school child: gen prev_dep=c_dep[_n-1]
sort school child wave

*STEP (vii):
* Simulate the values for the NAPLAN scores at waves 3,5 and 7
gen e_ijk=rnormal(0,0.5)

if "`model'"=="T1" {
gen napscore_z=2.0+ (-0.07)*prev_dep+ (-0.01)*wave+ 0.013*prev_dep*wave+ ///
			(-0.20)*c_age+ 0.14*c_gender+0.71*c_nap1_z+ (-0.01)*c_ses + b_0i + b_0ij+e_ijk if ///
			wave==3 |wave==5|wave==7
			}
			
if "`model'"=="T2" {
gen napscore_z=2.0+ (-0.024)*prev_dep+ (-0.01)*wave+ 0.023*prev_dep*c.c_ses+ ///
			(-0.20)*c_age+ 0.14*c_gender+0.71*c_nap1_z+ (-0.01)*c_ses + b_0i + b_0ij+e_ijk if ///
			wave==3 |wave==5|wave==7
			}
					
if "`model'"=="T3" {
gen prev_dep2=prev_dep*prev_dep
gen napscore_z=2.0+ (-0.024)*prev_dep+ (-0.01)*wave+ (-0.009)*prev_dep2+ ///
			(-0.20)*c_age+ 0.14*c_gender+0.71*c_nap1_z+ (-0.01)*c_ses + b_0i + b_0ij+e_ijk if ///
			wave==3 |wave==5|wave==7
drop prev_dep2
			}

*STEP (viii):
* Simulate values for the SDQ variable
gen e_p_sdq=rnormal(0,2.8)
gen p_sdq=16.2+2.5*c_dep+(-0.1)*wave+v_0i+v_0ij+e_p_sdq if wave==2|wave==4|wave==6
drop e_p_sdq

bysort school child: gen prev_sdq=p_sdq[_n-1]

*remove unwanted variables 
drop  b_0i b_0ij u_0i u_0ij v_0i v_0ij  e_ijk

************************************************************
*generating missing values 

*generate post wave naplan score
bysort school child: gen post_nap=napscore_z[_n+1]

*generate missing data inidcator in waves 2, 4 and 6 (0 missing 1 observed)
if "`mdm'"=="MAR1" {
gen r_dep=runiform()<invlogit(-7.5 +0.4*post_nap+0.7*p_sdq) if wave==2
replace r_dep=runiform()<invlogit(-7.7+0.4*post_nap+0.7*p_sdq) if wave==4
replace r_dep=runiform()<invlogit(-9.0+0.4*post_nap+0.7*p_sdq) if wave==6	
}

if "`mdm'"=="MAR2" {
gen r_dep=runiform()<invlogit(-15.5 +1.1*post_nap+1.4*p_sdq) if wave==2
replace r_dep=runiform()<invlogit(-16.5+1.1*post_nap+1.4*p_sdq) if wave==4
replace r_dep=runiform()<invlogit(-19.0+1.1*post_nap+1.4*p_sdq) if wave==6	
}
							
drop post_nap

*replace all other wave value indicators with 1							
replace r_dep=1 if r_dep==.
*tab r_dep

reshape wide  napscore_z c_dep prev_dep p_sdq prev_sdq r_dep , i(c_id) j(wave)

**check the proportions of missing data in wave 2,4 and 6
*tab r_dep2
*tab r_dep4
*tab r_dep6

**generate missing values in SES at wave 1 
gen u1ran= runiform()
gen r_SES=0 if u1ran<= 0.1
replace r_SES=1 if u1ran>0.1
drop u1ran

*tab r_SES

reshape long  napscore_z c_dep prev_dep p_sdq prev_sdq r_dep , i(c_id) j(wave)
replace  c_dep=. if r_dep==0
replace  c_ses=. if r_SES==0
drop r_dep r_SES prev_dep prev_sdq

*reshape wide napscore_z c_dep  p_sdq , i(c_id) j(wave)

label define gender_lbl 0 "female" 1 "male"
label values c_gender gender_lbl
label var child "Child's ID" 
label var c_id "Child's unique ID"
label var school "Child's school indicator"
label var c_age "Child's age at wave 1"
label var c_gender  "Child's sex"
label var c_ses "Standardized SE value at wave 1"
label var c_nap1_z " Standardized NAPLAN numeracy score at wave 1"
label var c_dep "Child depressive symptom score at waves 2,4 and 6 "
label var napscore_z "Standadized NAPLAN scores at waves 3,5 and 7"

*save the data set
export delimited using data`i'.csv,replace

}


end 
******************************  data generation ***************************


*Analysis model 1:Interaction between the time-varying exposure and time

* change working directory for each simulation scenario to a folder to save the 
*simulated datasets
cd "C:\temp"

*MAR-CATS, 40 school clusters
simul, seed(10032020) simno(1000) schoolclus(40) clussize(30) model(T1) mdm(MAR1)

*MAR-CATS, 10 school clusters(omitted 2 datasets i=306, i=660 due to nonconvergence in SMC-JM-3L)
simul, seed(20305687) simno(1002) schoolclus(10) clussize(120) model(T1) mdm(MAR1)

*MAR-inflated, 40 school clusters(omitted 1 dataset i=646 due to nonconvergence in SMC-JM-3L)
simul, seed(32305687) simno(1001) schoolclus(40) clussize(30) model(T1) mdm(MAR2)

*MAR-inflated, 10 school clusters
simul, seed(40315687) simno(1000) schoolclus(10) clussize(120) model(T1) mdm(MAR2)

******************************************************************************************
*Analysis model 2: Interaction between the time-varying exposure and a time-fixed baseline variable 

*MAR-CATS, 40 school clusters(omitted 1 dataset i=488 due to nonconvergence in SMC-JM-3L)
simul, seed(86014685) simno(1001) schoolclus(40) clussize(30) model(T2) mdm(MAR1)

*MAR-CATS, 10 school clusters(omitted 2 datasets i=500, i=184 due to nonconvergence in SMC-JM-3L)
simul, seed(91014685) simno(1002) schoolclus(10) clussize(120) model(T2) mdm(MAR1)

*MAR-inflated, 40 school clusters
simul, seed(67014685) simno(1000) schoolclus(40) clussize(30) model(T2) mdm(MAR2)

*MAR-inflated, 10 school clusters(omitted 3 datasets i=63,i=102,i=92 due to nonconvergence in SMC-JM-3L)
simul, seed(17092020) simno(1003) schoolclus(10) clussize(120) model(T2) mdm(MAR2)

******************************************************************************************
*Analysis model 3: Quadratic term in the exposure

*MAR-CATS, 40 school clusters(omitted 1 dataset i=667 due to nonconvergence in SMC-JM-3L)
simul, seed(16119020) simno(1001) schoolclus(40) clussize(30) model(T3) mdm(MAR1)

*MAR-CATS, 10 school clusters
simul, seed(10985432) simno(1000) schoolclus(10) clussize(120) model(T3) mdm(MAR1)

*MAR-inflated, 40 school clusters
simul, seed(14117020) simno(1000) schoolclus(40) clussize(30) model(T3) mdm(MAR2)

*MAR-inflated, 10 school clusters(omitted 1 dataset i=861 due to nonconvergence in SMC-JM-3L)
simul, seed(45985432) simno(1001) schoolclus(10) clussize(120) model(T3) mdm(MAR2)
