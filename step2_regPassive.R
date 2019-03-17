### @description
### this is the Step 2 of the simulation, which includes
### data generation, inducing missingness, application of MI methods with auxiliary variable fully observed
### @date: 10/24/2018
### @example command: Rscript step2_regPassive.R 100 3 1 0.5 high MAR 0.5 0.2 regPassive within TRUE TRUE 5 ../working_data/regPassive.txt

# clear work space and set working directory
rm(list=ls())
options(warn=-1)
setwd("/home/users/yling/thesis/20181231/scripts/")
#setwd("/Users/aling/Documents/PhD/the_last_push/weekly_progress/11212018/scripts/")

# takes input from command line
# nSimulation = 100; nBin = 2; nNor = 0; rho = 0.5; MDM = "MAR2"; missingLevel=0.5; caliper=0.2; reverse=FALSE;
# auxMDM="MCAR"; impMethod = "regPassive"; corr = "high"; intMethod="within"; MIwithAux=TRUE; auxMissing=TRUE; m=5; outfilePath="../working_data/derPassive_within_auxMissingMAR_n=100_exp1.txt"
args <- commandArgs(trailingOnly = TRUE)
nSimulation<-as.numeric(args[1])
nBin<-as.numeric(args[2])
nNor<-as.numeric(args[3])
rho<-as.numeric(args[4])
corr<-as.character(args[5])
MDM<-as.character(args[6])
missingLevel<-as.numeric(args[7])
caliper<-as.numeric(args[8])
impMethod<-as.character(args[9])
intMethod<-as.character(args[10])
MIwithAux<-as.logical(as.character(args[11]))
auxMissing<-as.logical(as.character(args[12]))
auxMDM<-as.character(args[13])
reverse<-as.logical(as.character(args[14]))
m<-as.numeric(args[15])
outfilePath<-as.character(args[16])

# load libraries
library("MASS")
library("Matching")
library("mice")

# load helper functions
source("helperFunctions.R")
source("is.r")
source("internal.r")
source("mice.r")

# specify parameters for simulation
startingSeed<-1
nSubj=2000
treatmentEffect=2
missingVariable<-"X2"
confounders=c("X1", "X2")

induceMissingness_regPassive<-function(seed=1, nSubj=2000, missingLevel=0.5, simulated.data, betas){
  simulated.data$psBinary<-ifelse(simulated.data$ps>=median(simulated.data$ps), 1, 0)

  missArray<-rep(0, nrow(simulated.data))
  missArray[which(simulated.data$psBinary==1 & simulated.data$treatment==1)]<-betas[1]
  missArray[which(simulated.data$psBinary==1 & simulated.data$treatment==0)]<-betas[2]
  missArray[which(simulated.data$psBinary==0 & simulated.data$treatment==1)]<-betas[3]
  missArray[which(simulated.data$psBinary==0 & simulated.data$treatment==0)]<-betas[4]

  # search for b0.missing that satisfies the missing level requirement
  b0.missings<-seq(-20, 20, 0.5)
  missingLevels<-{}
  for (b0.missing in b0.missings) {
    # the input for inv.log function is roughly (-6,6), if good outcome or on statin, almost ganrateed to be missing
    pi.miss.logit <- b0.missing + missArray
    pi.miss <- inv_logit(pi.miss.logit)
    miss <- rbinom(n=nSubj, size=1, p=pi.miss)
    missingLevels<-c(missingLevels, length(which(miss==1))/nSubj)
  }
  
  # find the optimal b0.missing
  b0.missing <- b0.missings[which.min(abs(missingLevels-missingLevel))]
  
  # induce missingness
  pi.miss.logit <- b0.missing + missArray
  pi.miss <- inv_logit(pi.miss.logit)
  miss <- rbinom(n=nSubj, size=1, p=pi.miss)
  
  # output
  return(miss)
}

induceMissingness_mar2_ps<-function(seed=1, nSubj=2000, missingLevel=0.5, simulated.data, betas, aux){
  simulated.data$goodOutcomeBinary<-ifelse(simulated.data$outcome>=median(simulated.data$outcome), 1, 0)
  simulated.data$highAux<-ifelse(aux>=median(aux), 1, 0)

  missArray<-rep(0, nrow(simulated.data))
  missArray[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==1)]<-betas[1]*as.numeric(as.character(simulated.data$highAux[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==1)]))
  missArray[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==0)]<-betas[2]*as.numeric(as.character(simulated.data$highAux[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==0)]))
  missArray[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==1)]<-betas[3]*as.numeric(as.character(simulated.data$highAux[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==1)]))
  missArray[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==0)]<-betas[4]*(1-as.numeric(as.character(simulated.data$highAux[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==0)])))

  # search for b0.missing that satisfies the missing level requirement
  b0.missings<-seq(-20, 20, 0.5)
  missingLevels<-{}
  for (b0.missing in b0.missings) {
    # the input for inv.log function is roughly (-6,6), if good outcome or on statin, almost ganrateed to be missing
    pi.miss.logit <- b0.missing + missArray
    pi.miss <- inv_logit(pi.miss.logit)
    miss <- rbinom(n=nSubj, size=1, p=pi.miss)
    missingLevels<-c(missingLevels, length(which(miss==1))/nSubj)
  }
  
  # find the optimal b0.missing
  b0.missing <- b0.missings[which.min(abs(missingLevels-missingLevel))]
  
  # induce missingness
  pi.miss.logit <- b0.missing + missArray
  pi.miss <- inv_logit(pi.miss.logit)
  miss <- rbinom(n=nSubj, size=1, p=pi.miss)
  
  # output
  return(miss)
}


# remove existing files
system(paste("rm ", outfilePath, sep=""))
# run simulation nSimulation times
for (i in 1:nSimulation){
	print(paste("Starting simulation ", i, sep=""))

	# set seed
	seed = startingSeed + i - 1
	set.seed(seed)

	# simulated.data
	covariates <- simulateCovariate(seed=seed, nSubj=nSubj, nBin=nBin, nNor=nNor, rho=rho)
	treatment <- simulateTreatment(seed=seed, nSubj=nSubj, betas=c(2,2), covariates=covariates)
	outcome <- simulateContinousOutcome(seed=seed, nSubj=nSubj, betas=c(2,2), covariates=covariates, treatment=treatment, treatmentEffect=treatmentEffect, sigma=10)
	simulated.data<-as.data.frame(covariates)
	simulated.data$treatment = treatment
	simulated.data$outcome = outcome

	simulated.data$ps<-glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data, family = binomial)$fitted
	X2aux = 1 + 10*as.numeric(as.character(simulated.data[,missingVariable])) + rnorm(nSubj, 0, 1) # correlation is 0.99 
	aux<-10*simulated.data$ps+rnorm(nSubj,0,1)

	# cor(aux, simulated.data$ps)
	
	#cor(simulated.data$X3, simulated.data$ps)
	#cor(simulated.data$X1, simulated.data$X2)
	#cor(simulated.data$X2, simulated.data$X3)
	#cor(simulated.data$X1, simulated.data$X3)

	#cor(simulated.data$X1, simulated.data$treatment)
	#cor(simulated.data$X2, simulated.data$treatment)
	#cor(simulated.data$X3, simulated.data$treatment)

	#cor(simulated.data$X1, simulated.data$outcome)
	#cor(simulated.data$X2, simulated.data$outcome)
	#cor(simulated.data$treatment, simulated.data$outcome)

	#cor(simulated.data$X1, simulated.data$ps)
	#cor(simulated.data$X2, simulated.data$ps)
	#cor(simulated.data$X3, simulated.data$ps)
	#cor(simulated.data$ps, simulated.data$treatment)

	#par(mfrow=c(1,2))
	#hist(simulated.data$ps[which(simulated.data$treatment==1)], main="Treated", xlab="PS")
	#hist(simulated.data$ps[which(simulated.data$treatment==0)], main="Controls", xlab="PS")

	# format simulated data
	if (nBin==1) { simulated.data[,1]<-as.factor(simulated.data[,1]) }
	if (nBin==2) { 
		simulated.data[,1]<-as.factor(simulated.data[,1])
		simulated.data[,2]<-as.factor(simulated.data[,2])
	}

	# induce missingness
	missColname = "miss"
	simulated.data$miss<-induceMissingness_mar2(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1), aux=aux)
	#simulated.data$miss<-0
	#simulated.data$miss[which(simulated.data$treatment==1 & simulated.data$outcome>=5)]<-1
	#simulated.data$miss[which(simulated.data$treatment==0 & simulated.data$outcome<=0)]<-1
	#simulated.data$miss[which(simulated.data$treatment==0 & simulated.data$outcome>=5)]<-1


	#simulated.data$miss<-induceMissingness_regPassive(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1))
	  
	# multiple imputation
	complete.data <- simulated.data[,c(1:4)]
	incomplete.data <- complete.data
	incomplete.data[which(simulated.data$miss==1),missingVariable]<-NA

    # add auxiliary variable
    missingVariables = missingVariable
	if (MIwithAux==TRUE){
		if (auxMissing==TRUE) { 
			if (auxMDM=="MCAR") { 
				aux[sample(1:nSubj, nSubj*0.2, replace = FALSE)] <- NA } else if (auxMDM=="MAR1") { 
				miss_aux<-induceMissingness_regPassive(seed=seed, nSubj=nSubj, missingLevel=0.2, simulated.data=simulated.data, betas=c(5,0,0,5)) 
				aux[which(miss_aux==1)]<-NA
			} else if (auxMDM=="MAR2") {
				miss_aux<-induceMissingness_regPassive(seed=seed, nSubj=nSubj, missingLevel=0.2, simulated.data=simulated.data, betas=c(0,5,5,0)) 
				aux[which(miss_aux==1)]<-NA
			} else if (auxMDM=="MAR3") {
				miss_aux<-induceMissingness_regPassive(seed=seed, nSubj=nSubj, missingLevel=0.2, simulated.data=simulated.data, betas=c(5,0,0,0)) 
				aux[which(miss_aux==1)]<-NA
			} else if (auxMDM=="MAR4") {
				miss_aux<-induceMissingness_regPassive(seed=seed, nSubj=nSubj, missingLevel=0.2, simulated.data=simulated.data, betas=c(0,0,0,5))
				aux[which(miss_aux==1)]<-NA 
			}
			missingVariables <- c(missingVariable, "aux")
		}
		incomplete.data$aux <- aux
	}

	# simulated.data=incomplete.data; nCovariates=nBin+nNor; outcomeBetas=0
	r<-MIwithPSmatching(complete.data=complete.data, simulated.data=incomplete.data, impMethod=impMethod, intMethod=intMethod, MIwithAux=MIwithAux, missingVariables=missingVariables, confounders=confounders, nCovariates=nBin+nNor, outcomeBetas=0, nBin=nBin, nNor=nNor, caliper=caliper, m=m, reverse=reverse)

	if (intMethod=="within"){
	    for (index in 1:m) {
	    write(paste("s", seed, i, impMethod, intMethod, MIwithAux, paste(r[[1]],collapse=","), paste(unlist(r[[2]][index]), collapse=","), paste(unlist(r[[3]][index]), collapse=","), paste(unlist(r[[4]][index]), collapse=","), paste(unlist(r[[5]][index]), collapse=","), paste(unlist(r[[6]][index]), collapse=","), sep=","), file=outfilePath, append=TRUE)
	    }  
	} else {
		write(paste("s", seed, i, impMethod, intMethod, MIwithAux, paste(r[[1]],collapse=","), paste(unlist(r[[2]]), collapse=","), paste(unlist(r[[3]]), collapse=","), paste(unlist(r[[4]]), collapse=","), paste(unlist(r[[5]]), collapse=","), paste(unlist(r[[6]]), collapse=","), sep=","), file=outfilePath, append=TRUE)
	} 
}








