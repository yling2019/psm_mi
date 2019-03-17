### @description
### this is the Step 2 of the simulation, which includes
### data generation, inducing missingness, application of MI methods with auxiliary variable fully observed
### @date: 10/24/2018
### @example command: Rscript step2.R 1000 2 0 0.5 MAR 0.5 0.2 passive within TRUE TRUE ../working_data/test.txt

# clear work space and set working directory
rm(list=ls())
options(warn=-1)
setwd("/home/users/yling/thesis/20181231/scripts")
#setwd("/Users/aling/Documents/PhD/the_last_push/weekly_progress/12312018/scripts/")

# takes input from command line
# nSimulation = 50; nBin = 2; nNor = 0; rho = 0.5; MDM = "MCAR"; missingLevel=0.5; caliper=0.2; intMethod="within"; impMethod="passive"; MIwithAux=TRUE; auxMissing=FALSE; b=2000;startingSeed=951; 
args <- commandArgs(trailingOnly = TRUE)
nSimulation<-as.numeric(args[1])
nBin<-as.numeric(args[2])
nNor<-as.numeric(args[3])
rho<-as.numeric(args[4])
corr<-"high"
MDM<-as.character(args[5])
missingLevel<-as.numeric(args[6])
caliper<-as.numeric(args[7])
impMethod<-as.character(args[8])
intMethod<-as.character(args[9])
MIwithAux<-TRUE
auxMissing<-as.logical(as.character(args[10]))
auxMDM<-as.character(args[11])
reverse<-as.logical(as.character(args[12]))
m<-5
startingSeed<-as.numeric(args[13])
outfilePath<-as.character(args[14])

b=2000
# remove existing files
system(paste("rm ", outfilePath, sep=""))

# load libraries
library("MASS")
library("Matching")
library("mice")
library("doParallel")
library("foreach")

# load helper functions
source("helperFunctions.R")
source("is.r")
source("internal.r")
source("mice.r")

# specify parameters for simulation
nSubj=2000
treatmentEffect=2
missingVariable<-"X2"
confounders=c("X1", "X2")
nSubsample=1000

registerDoParallel(cores=16)

# run simulation nSimulation times
for (i in 1:nSimulation){
	print(paste("Starting simulation ", i, sep=""))

	# set seed
	seed = startingSeed + i - 1
	set.seed(seed)

	# simulate data
	covariates <- simulateCovariate(seed=seed, nSubj=nSubj, nBin=nBin, nNor=nNor, rho=rho)
	treatment <- simulateTreatment(seed=seed, nSubj=nSubj, betas=c(2,2), covariates=covariates)
	outcome <- simulateContinousOutcome(seed=seed, nSubj=nSubj, betas=c(2,2), covariates=covariates, treatment=treatment, treatmentEffect=treatmentEffect, sigma=10)
	simulated.data<-as.data.frame(covariates)
	simulated.data$treatment = treatment
	simulated.data$outcome = outcome

	# format simulated data
	if (nBin==1) { simulated.data[,1]<-as.factor(simulated.data[,1]) }
	if (nBin==2) { 
		simulated.data[,1]<-as.factor(simulated.data[,1])
		simulated.data[,2]<-as.factor(simulated.data[,2])
	}

	# induce missingness
	missColname = "miss"
	aux = 1 + 10*as.numeric(as.character(simulated.data[,missingVariable])) + rnorm(nSubj, 0, 1) 

	if (MDM=="MCAR") {
		simulated.data$miss<-0
		simulated.data$miss[sample(1:nSubj, nSubj*missingLevel, replace = FALSE)]<-1
	} else if (MDM=="MAR1"){
		simulated.data$miss<-induceMissingness_mar1(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1))
	} else if (MDM=="MAR2"){
		simulated.data$miss<-induceMissingness_mar2(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1), aux=aux)
	} else if (MDM=="MNAR"){
		simulated.data$miss<-induceMissingness_mnar(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1))
	}
	  
	# multiple imputation
	complete.data <- simulated.data[,1:4]
	incomplete.data <- complete.data
	incomplete.data[which(simulated.data$miss==1),missingVariable]<-NA

    # add auxiliary variable
    missingVariables = missingVariable
	if (MIwithAux==TRUE){
		aux = 1 + 10*as.numeric(as.character(complete.data[,missingVariable])) + rnorm(nSubj, 0, 1) 
		if (auxMissing==TRUE) { 
			aux[sample(1:nSubj, nSubj*0.2, replace = FALSE)] <- NA 
			missingVariables <- c(missingVariable, "aux")
		}
		incomplete.data$aux <- aux
	}

	if (b==nSubj) {
		r = foreach( i = 1:nSubsample, .combine=rbind ) %dopar% {
			ids = sample(1:nSubj, nSubj, replace=TRUE)
			cd <- complete.data[ids,]
			id <- incomplete.data[ids,]
			MIwithPSmatching(complete.data=cd, simulated.data=id, impMethod=impMethod, intMethod=intMethod, MIwithAux=MIwithAux, missingVariables=missingVariables, confounders=confounders, nCovariates=nBin+nNor, outcomeBetas=0, nBin=nBin, nNor=nNor, caliper=caliper)
		}
	} else {
		r = foreach( i = 1:nSubsample, .combine=rbind ) %dopar% {
			ids = sample(1:nSubj, b, replace=FALSE)
			cd <- complete.data[ids,]
			id <- incomplete.data[ids,]
			MIwithPSmatching(complete.data=cd, simulated.data=id, impMethod=impMethod, intMethod=intMethod, MIwithAux=MIwithAux, missingVariables=missingVariables, confounders=confounders, nCovariates=nBin+nNor, outcomeBetas=0, nBin=nBin, nNor=nNor, caliper=caliper)
		}
	}


	for (x in 1:nSubsample){
		write(paste("s", seed, i, x, impMethod, intMethod, MIwithAux, paste(r[,1][x][[1]][1], collapse=","), sep=","), file=outfilePath, append=TRUE)
	}
}
