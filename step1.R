### @description
### this is the Step 1 of the simulation, which includes
### data generation, inducing missingness, full data analysis, and application of commonly applied missing data methods
### @date: 10/24/2018
### @example command: Rscript step1.R 100 2 0 0.5 NMAR 0.5 0.2 NMAR.txt

# clear work space and set working directory
rm(list=ls())
options(warn=-1)
# setwd("/home/users/yling/thesis/20181024/scripts/")
setwd("/home/users/yling/thesis/20181231/scripts")

# takes input from command line
# nSimulation = 1000; nBin = 2; nNor = 0; rho = 0.5; MDM = "MCAR"; missingLevel=0.5; caliper=0.2
args <- commandArgs(trailingOnly = TRUE)
nSimulation<-as.numeric(args[1])
nBin<-as.numeric(args[2])
nNor<-as.numeric(args[3])
rho<-as.numeric(args[4])
MDM<-as.character(args[5])
missingLevel<-as.numeric(args[6])
caliper<-as.numeric(args[7])
outfilePath<-as.character(args[8])

# remove existing files
system(paste("rm ", outfilePath, sep=""))

# load libraries
library("MASS")
library("Matching")

# load helper functions
source("helperFunctions.R")

# specify parameters for simulation
startingSeed<-1
nSubj=2000
treatmentEffect=2
missingVariable<-"X2"
confounders=c("X1", "X2")

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
	if (MDM=="MCAR") {
		simulated.data$miss<-0
		simulated.data$miss[sample(1:nSubj, nSubj*missingLevel, replace = FALSE)]<-1
	} else if (MDM=="MAR1"){
		simulated.data$miss<-induceMissingness_mar1(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1))
	} else if (MDM=="MAR2"){
		aux = 1 + 10*as.numeric(as.character(simulated.data[,missingVariable])) + rnorm(nSubj, 0, 1) 
		simulated.data$miss<-induceMissingness_mar2(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1), aux=aux)
	} else if (MDM=="MNAR"){
		simulated.data$miss<-induceMissingness_mnar(seed=seed, nSubj=nSubj, missingLevel=missingLevel, simulated.data=simulated.data, betas=c(5,0,0,1))
	}

	# estimate treatment effect using true regression model on full data
	model1<-glm(outcome ~ . , data = simulated.data[,1:(ncol(covariates) + 2)])
  	ps.model1 <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data, family = binomial)
  	data.matched1 <- Match(Y=simulated.data$outcome, Tr=simulated.data$treatment, X=ps.model1$fitted, M=1, caliper = caliper, ties=FALSE, replace=FALSE)
  	balance1<-plotStandardizedDiff(simulated.data, data.matched1, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)
  	balance1.5<-plotStandardizedDiff(simulated.data, data.matched1, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)
  	#treatment.estimate1<-mean(simulated.data$outcome[data.matched1$index.treated]) - mean(simulated.data$outcome[data.matched1$index.control])
  	#var.estimate1<-varianceEstimate(simulated.data, data.matched1)
	post_matching<-simulated.data[c(data.matched1$index.treated, data.matched1$index.control),]
	post_matching$pair_id<-rep(1:length(data.matched1$index.treated), 2)
	dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")
 	r1<-c(summary(model1)$coefficients["treatment",1], summary(model1)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data[which(simulated.data$treatment==1),]), nrow(simulated.data[which(simulated.data$treatment==0),]), balance1)

	# complete case analysis (CC)
	simulated.data2<-simulated.data[which(simulated.data[,missColname]==0),]
	model2<-glm(outcome ~ . , data = simulated.data2[,1:(ncol(covariates) + 2)])
  
	ps.model2 <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data2, family = binomial)
	data.matched2 <- Match(Y=simulated.data2$outcome, Tr=simulated.data2$treatment, X=ps.model2$fitted, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
	#treatment.estimate2<-mean(simulated.data2$outcome[data.matched2$index.treated]) - mean(simulated.data2$outcome[data.matched2$index.control])
	#var.estimate2<-varianceEstimate(simulated.data2, data.matched2)
	balance2<-plotStandardizedDiff(simulated.data2, data.matched2, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)

	post_matching<-simulated.data2[c(data.matched2$index.treated, data.matched2$index.control),]
	post_matching$pair_id<-rep(1:length(data.matched2$index.treated), 2)
	dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")
	r2<-c(summary(model2)$coefficients["treatment",1], summary(model2)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data2[which(simulated.data2$treatment==1),]), nrow(simulated.data2[which(simulated.data2$treatment==0),]), balance2)

	# complete variable analysis (CVA) by excluding missing variables
	simulated.data3<-simulated.data # to avoid duplication of pr_scores
  	model3<-glm(outcome ~ . , data = simulated.data3[,1:(ncol(covariates) + 2)][,-which(colnames(simulated.data3)==missingVariable)])
  
	confounders_minus<-setdiff(confounders, missingVariable)
	ps.model3 <- glm(as.formula(paste("treatment ~", paste(confounders_minus, collapse="+"))), data=simulated.data3, family = binomial)
	data.matched3 <- Match(Y=simulated.data3$outcome, Tr=simulated.data3$treatment, X=ps.model3$fitted, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
	#treatment.estimate3<-mean(simulated.data3$outcome[data.matched3$index.treated]) - mean(simulated.data3$outcome[data.matched3$index.control])
	#var.estimate3<-varianceEstimate(simulated.data3, data.matched3)
	balance3<-plotStandardizedDiff(simulated.data3, data.matched3, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)
  	balance3.5<-plotStandardizedDiff(simulated.data, data.matched3, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)

  	post_matching<-simulated.data3[c(data.matched3$index.treated, data.matched3$index.control),]
	post_matching$pair_id<-rep(1:length(data.matched3$index.treated), 2)
	dr<-ols(outcome ~ treatment + X1, data=post_matching, cluster="pair_id")

	r3<-c(summary(model3)$coefficients["treatment",1], summary(model3)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data3[which(simulated.data3$treatment==1),]), nrow(simulated.data3[which(simulated.data3$treatment==0),]), balance3)
  	r3.5<-c(summary(model3)$coefficients["treatment",1], summary(model3)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data3[which(simulated.data3$treatment==1),]), nrow(simulated.data3[which(simulated.data3$treatment==0),]), balance3.5)

	# single imputation -- mean imputation
	simulated.data4<-simulated.data # to avoid duplication of pr_scores
	simulated.data4$X2<-as.numeric(as.character(simulated.data4$X2))
  
	simulated.data4[,missingVariable][which(simulated.data4[,missColname]==1)] <- mean(as.numeric(as.character(simulated.data4[which(simulated.data4[,missColname]==0),missingVariable])))
	model4<-glm(outcome ~ . , data = simulated.data4[,1:(ncol(covariates) + 2)])
  
	ps.model4 <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data4, family = binomial)
	data.matched4 <- Match(Y=simulated.data4$outcome, Tr=simulated.data4$treatment, X=ps.model4$fitted, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
	#treatment.estimate4<-mean(simulated.data4$outcome[data.matched4$index.treated]) - mean(simulated.data4$outcome[data.matched4$index.control])
	#var.estimate4<-varianceEstimate(simulated.data4, data.matched4)
	balance4<-plotStandardizedDiff(simulated.data4, data.matched4, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)
  	balance4.5<-plotStandardizedDiff(simulated.data, data.matched4, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=NULL, plot=FALSE, return=TRUE)

  	post_matching<-simulated.data4[c(data.matched4$index.treated, data.matched4$index.control),]
	post_matching$pair_id<-rep(1:length(data.matched4$index.treated), 2)
	dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")

	r4<-c(summary(model4)$coefficients["treatment",1], summary(model4)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data4[which(simulated.data4$treatment==1),]), nrow(simulated.data4[which(simulated.data4$treatment==0),]), balance4)
	r4.5<-c(summary(model4)$coefficients["treatment",1], summary(model4)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data4[which(simulated.data4$treatment==1),]), nrow(simulated.data4[which(simulated.data4$treatment==0),]), balance4.5)

  
	# single imputation -- create a missing indicator
	simulated.data5<-simulated.data # to avoid duplication of pr_scores
	  
	if (length(unique(simulated.data[,missingVariable]))==2) {
		simulated.data5[,1]<-as.numeric(as.character(simulated.data5[,1]))
		simulated.data5[,2]<-as.numeric(as.character(simulated.data5[,2]))
		simulated.data5[which(simulated.data5[,missColname]==1),missingVariable] <- 2
		simulated.data5[,1]<-as.factor(simulated.data5[,1])
		simulated.data5[,2]<-as.factor(simulated.data5[,2])

		model5<-glm(outcome ~ . , data = simulated.data5[,c(1:(ncol(covariates) + 2))])
		ps.model5 <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data5, family = binomial)
		data.matched5 <- Match(Y=simulated.data5$outcome, Tr=simulated.data5$treatment, X=ps.model5$fitted, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
		#treatment.estimate5<-mean(simulated.data5$outcome[data.matched5$index.treated]) - mean(simulated.data5$outcome[data.matched5$index.control])
		#var.estimate5<-varianceEstimate(simulated.data5, data.matched5)
		balance5<-plotStandardizedDiff(simulated.data, data.matched5, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
	
		data.matched5.missing<-list(index.treated=intersect(which(simulated.data$miss==1), data.matched5$index.treated), index.control=intersect(which(simulated.data$miss==1), data.matched5$index.control))
		data.matched5.nonmiss<-list(index.treated=intersect(which(simulated.data$miss==0), data.matched5$index.treated), index.control=intersect(which(simulated.data$miss==0), data.matched5$index.control))

		# note that the pre-matching balance in the below two calcualations are off. Should put in simulated.data[which(simulated.data$miss==1),0] for example
		balance5.missing<-plotStandardizedDiff(simulated.data, data.matched5.missing, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
		balance5.nonmiss<-plotStandardizedDiff(simulated.data, data.matched5.nonmiss, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)

		post_matching<-simulated.data5[c(data.matched5$index.treated, data.matched5$index.control),]
		post_matching$pair_id<-rep(1:length(data.matched5$index.treated), 2)
		dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")

	} else {
		simulated.data5[,missColname] <- as.factor(simulated.data5[,missColname])
		simulated.data5[,missingVariable][which(simulated.data5[,missColname]==1)] <- mean(simulated.data5[,missingVariable][which(simulated.data5[,missColname]==0)])

		model5<-glm(outcome ~ . , data = simulated.data5[,c(1:(ncol(covariates) + 2), grep(missColname, colnames(simulated.data5)))])
		ps.model5 <- glm(as.formula(paste("treatment ~", paste(c(confounders, missColname), collapse="+"))), data=simulated.data5, family = binomial)
		data.matched5 <- Match(Tr=simulated.data5$treatment, X=ps.model5$fitted, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
		#treatment.estimate5<-mean(simulated.data5$outcome[data.matched5$index.treated]) - mean(simulated.data5$outcome[data.matched5$index.control])
		#var.estimate5<-varianceEstimate(simulated.data5, data.matched5)
		balance5<-plotStandardizedDiff(simulated.data, data.matched5, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
		
		data.matched5.missing<-list(index.treated=intersect(which(simulated.data$miss==1), data.matched5$index.treated), index.control=intersect(which(simulated.data$miss==1), data.matched5$index.control))
		data.matched5.nonmiss<-list(index.treated=intersect(which(simulated.data$miss==0), data.matched5$index.treated), index.control=intersect(which(simulated.data$miss==0), data.matched5$index.control))

		# note that the pre-matching balance in the below two calcualations are off. Should put in simulated.data[which(simulated.data$miss==1),0] for example
		balance5.missing<-plotStandardizedDiff(simulated.data, data.matched5.missing, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
		balance5.nonmiss<-plotStandardizedDiff(simulated.data, data.matched5.nonmiss, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)

		post_matching<-simulated.data5[c(data.matched5$index.treated, data.matched5$index.control),]
		post_matching$pair_id<-rep(1:length(data.matched5$index.treated), 2)
		dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")
  	}  

	# gather results of direct modeling and PS matching results for indicator variable
	r5<-c(summary(model5)$coefficients["treatment",1], summary(model5)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data5[which(simulated.data5$treatment==1),]), nrow(simulated.data5[which(simulated.data5$treatment==0),]), balance5)
  	r6<-c(summary(model5)$coefficients["treatment",1], summary(model5)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data5[which(simulated.data5$treatment==1),]), nrow(simulated.data5[which(simulated.data5$treatment==0),]), balance5.missing)
	r7<-c(summary(model5)$coefficients["treatment",1], summary(model5)$coefficients["treatment",2], dr[1], dr[2], nrow(simulated.data5[which(simulated.data5$treatment==1),]), nrow(simulated.data5[which(simulated.data5$treatment==0),]), balance5.nonmiss)

	# output diagnosis statistics for propensity score matching2
	write(paste(c("s", seed, i, "fullData", r1), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "CC", r2), collapse=","), file=outfilePath, append=TRUE) # note that the number of patients pre-matching is 50% of nSubj
	write(paste(c("s", seed, i, "CVA", r3), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "CVA", r3.5), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "meanImputation", r4), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "meanImputation", r4.5), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "indicatorImputation", r5), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "indicatorImputation", r6), collapse=","), file=outfilePath, append=TRUE)
	write(paste(c("s", seed, i, "indicatorImputation", r7), collapse=","), file=outfilePath, append=TRUE)
}
