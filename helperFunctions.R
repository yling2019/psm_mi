# this file includes all internal helper functions, including
# data generation, propensity score matching, and various MI methods
simulateCovariate<-function(seed=1, nSubj=2000, nBin=1, nNor=1, rho=0.5){
	set.seed(seed)

  nCovariates = nBin + nNor
	data <- as.data.frame(mvrnorm(n=nSubj, mu = rep(0, nCovariates), Sigma=matrix(c(1, rho, rho, 1),nCovariates,nCovariates)))

	if (nBin==0) { out <- data[,1:nCovariates] } else {
    for (x in 1:nBin) {
      data[,(ncol(data)+1)]<-ifelse(data[,x]>=0, 1, 0)
    }
    if (nNor>0){
      out <- data[,c((nCovariates+1) : (nCovariates+nBin), (nCovariates-nNor+1):nCovariates)]
      } else {
        out <- data[,c((nCovariates+1) : (nCovariates+nBin))]
      }
  }

	colnames(out)<-paste("X", 1:nCovariates, sep="")
	return(out)
}

inv_logit <- function(x) { 
  out <- exp(x)/(1+exp(x))
  return(out) 
}

simulateTreatment<-function(seed=1, nSubj=2000, betas, covariates, propTreated=0.3){
  if (length(betas) != ncol(covariates)) { print("Beta and covariate dimensions do not match!") }
  
  set.seed(seed)
  pi.treatment.logit <- rep(0, nSubj)
  for (i in 1:ncol(covariates)){
    pi.treatment.logit <- pi.treatment.logit + betas[i]*covariates[,i]
  }

  # calculate beta0 based on proportion treated input argument
  range<-seq(-20,20,0.05)
  numTreated<-rep(0,length(range))
  for (j in 1:length(range)){
    tmp1 <- pi.treatment.logit + range[j]
    tmp2 <- inv_logit(tmp1)
    tmp3 <- rbinom(n=nSubj, size=1, p=tmp2)
    numTreated[j]<-length(which(tmp3==1))
  }
  beta0=range[which.min(abs(numTreated/nSubj-propTreated))] 
  
  pi.treatment <- inv_logit(pi.treatment.logit + beta0)
  treatment <- rbinom(n=nSubj, size=1, p=pi.treatment)
  
  return(treatment)
}


simulateContinousOutcome<-function(seed=1, nSubj=2000, betas, covariates, treatment, treatmentEffect, sigma, intercept = 0){
  set.seed(seed)
  
  outcome <- intercept + treatment*treatmentEffect + rnorm(nSubj, 0, sigma)
  for (i in 1:ncol(covariates)){
    outcome <- outcome + betas[i]*covariates[,i]
  }
  
  return(outcome)
}

# reference: https://diffuseprior.wordpress.com/2012/06/15/standard-robust-and-clustered-standard-errors-computed-in-r/
# robust cluster SE estimator
ols <- function(form, data, robust=FALSE, cluster=NULL,digits=4){
  r1 <- lm(form, data)
  if(length(cluster)!=0){
    data <- na.omit(data[,c(colnames(r1$model),cluster)])
    r1 <- lm(form, data)
  }
  X <- model.matrix(r1)
  n <- dim(X)[1]
  k <- dim(X)[2]
  if(robust==FALSE & length(cluster)==0){
    se <- sqrt(diag(solve(crossprod(X)) * as.numeric(crossprod(resid(r1))/(n-k))))
    res <- cbind(coef(r1),se)
  }
  if(robust==TRUE){
    u <- matrix(resid(r1))
    meat1 <- t(X) %*% diag(diag(crossprod(t(u)))) %*% X
    dfc <- n/(n-k)    
    se <- sqrt(dfc*diag(solve(crossprod(X)) %*% meat1 %*% solve(crossprod(X))))
    res <- cbind(coef(r1),se)
    }
  if(length(cluster)!=0){
    clus <- cbind(X,data[,cluster],resid(r1))
    colnames(clus)[(dim(clus)[2]-1):dim(clus)[2]] <- c(cluster,"resid")
    m <- dim(table(clus[,cluster]))
    dfc <- (m/(m-1))*((n-1)/(n-k))
    uclust  <- apply(resid(r1)*X,2, function(x) tapply(x, clus[,cluster], sum))
    se <- sqrt(diag(solve(crossprod(X)) %*% (t(uclust) %*% uclust) %*% solve(crossprod(X)))*dfc)   
    res <- cbind(coef(r1),se)
  }
  res <- cbind(res,res[,1]/res[,2],(1-pnorm(abs(res[,1]/res[,2])))*2)
  res1 <- matrix(as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),res)),nrow=dim(res)[1])
  rownames(res1) <- rownames(res)
  colnames(res1) <- c("Estimate","Std. Error","t value","Pr(>|t|)")
  return(res1[2,1:2])
}

# reference: Austin, Peter C. "Type I error rates, coverage of confidence intervals, and variance estimation in propensity-score matched analyses." The international journal of biostatistics 5.1 (2009).
# reference: https://www.quora.com/In-R-what-is-the-difference-between-dt-pt-and-qt-in-reference-to-the-student-t-distribution
varianceEstimate<-function(simulated.data, data.matched){
  ds<-simulated.data$outcome[data.matched$index.treated]-simulated.data$outcome[data.matched$index.control]
  te<-mean(simulated.data$outcome[data.matched$index.treated]) - mean(simulated.data$outcome[data.matched$index.control])
  stand.error<-sd(ds)/sqrt(length(ds))
  lower<-te -  qt(0.975, df=(length(ds)-1))*stand.error
  upper<-te +  qt(0.975, df=(length(ds)-1))*stand.error
  
  return(list(se=stand.error, lower=lower, upper=upper))
}


standardizedDifferenceB<-function(p1, p2){ #1 treatment 2 control
  d=(p1-p2)/sqrt((p1*(1-p1)+p2*(1-p2))/2)
  return(d)
}

standardizedDifferenceC<-function(mu1, mu2, var1, var2){ #1 treatment 2 control
  d=(mu1-mu2)/sqrt( (var1+var2)/2 )
  return(d)
}

plotStandardizedDiff<-function(data, data.matched, nBin, nNor, nCovariates, return=FALSE, confounders, outcomeBetas, main=NULL, plot=FALSE) {
  n=nBin + nNor

  before<-rep(NA, n)
  if (nBin>0 & nNor>0){
    for (i in 1:nBin){
      var = colnames(data)[i]
      before[i]<-standardizedDifferenceB(mean(as.numeric(as.character(data[which(data$treatment==1), var])), na.rm=TRUE), mean(as.numeric(as.character(data[which(data$treatment==0), var])), na.rm=TRUE))
    }
    for (i in (nBin+1):nCovariates){
      var = colnames(data)[i]
      before[i]<-standardizedDifferenceC(mean(data[which(data$treatment==1), var], na.rm=TRUE), mean(data[which(data$treatment==0), var], na.rm=TRUE), var(data[which(data$treatment==1), var], na.rm=TRUE), var(data[which(data$treatment==0), var], na.rm=TRUE))
    }
  } else if (nBin>0 & nNor==0){
    for (i in 1:nCovariates){
      var = colnames(data)[i]
      before[i]<-standardizedDifferenceB(mean(as.numeric(as.character(data[which(data$treatment==1), var])), na.rm=TRUE), mean(as.numeric(as.character(data[which(data$treatment==0), var])), na.rm=TRUE))
    }
  } else if (nBin==0 & nNor>0){
    for (i in 1:nCovariates){
      var = colnames(data)[i]
      before[i]<-standardizedDifferenceC(mean(data[which(data$treatment==1), var], na.rm=TRUE), mean(data[which(data$treatment==0), var], na.rm=TRUE), var(data[which(data$treatment==1), var], na.rm=TRUE), var(data[which(data$treatment==0), var], na.rm=TRUE))
    }
  }
  
  after<-rep(NA, n)
  if (nBin>0 & nNor>0){
    for (i in 1:nBin){
      var = colnames(data)[i]
      after[i]<-standardizedDifferenceB(mean(as.numeric(as.character(data[data.matched$index.treated, var])), na.rm=TRUE), mean(as.numeric(as.character(data[data.matched$index.control, var])), na.rm=TRUE))
    }
    for (i in (nBin+1):nCovariates){
      var = colnames(data)[i]
      after[i]<-standardizedDifferenceC(mean(data[data.matched$index.treated, var], na.rm=TRUE), mean(data[data.matched$index.control, var], na.rm=TRUE), var(data[data.matched$index.treated, var], na.rm=TRUE), var(data[data.matched$index.control, var], na.rm=TRUE))
    }
  } else if (nBin>0 & nNor==0){
    for (i in 1:nCovariates){
      var = colnames(data)[i]
      after[i]<-standardizedDifferenceB(mean(as.numeric(as.character(data[data.matched$index.treated, var])), na.rm=TRUE), mean(as.numeric(as.character(data[data.matched$index.control, var])), na.rm=TRUE))
    }
  } else if (nBin==0 & nNor>0){
    for (i in 1:nCovariates){
      var = colnames(data)[i]
      after[i]<-standardizedDifferenceC(mean(data[data.matched$index.treated, var], na.rm=TRUE), mean(data[data.matched$index.control, var], na.rm=TRUE), var(data[data.matched$index.treated, var], na.rm=TRUE), var(data[data.matched$index.control, var], na.rm=TRUE))
    }
  }
  
  if (plot==TRUE) {
    if (is.null(main)) { main="Matching Complete Dataset" }
    plot(0, bty='n', pch='', xlab='Standardized Difference', ylab='Variables', xlim=c(-1,1), ylim=c(0, (nBin+nNor)), main=main)
    for (i in which(outcomeBetas!=0)){
      arrows(x0=before[i], x1=after[i], y0=i, y1=i, xlim=range(c(before, after)), ylim=c(0, (nBin+nNor)), length=0.1)
    }
    abline(v=0.1, lty=3)
    abline(v=-0.1, lty=3)
  }
  
  if (return==TRUE) { return(list=c(length(data.matched$index.treated), length(data.matched$index.control), before=before, after=after)) }
  
}

induceMissingness_mar1<-function(seed=1, nSubj=2000, missingLevel=0.5, simulated.data, betas){
  simulated.data$goodOutcomeBinary<-ifelse(simulated.data$outcome>=median(simulated.data$outcome), 1, 0)

  missArray<-rep(0, nrow(simulated.data))
  missArray[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==1)]<-betas[1]
  missArray[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==0)]<-betas[2]
  missArray[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==1)]<-betas[3]
  missArray[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==0)]<-betas[4]

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

induceMissingness_mar2<-function(seed=1, nSubj=2000, missingLevel=0.5, simulated.data, betas, aux){
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

induceMissingness_mnar<-function(seed=1, nSubj=2000, missingLevel=0.5, simulated.data, betas){
  simulated.data$goodOutcomeBinary<-ifelse(simulated.data$outcome>=median(simulated.data$outcome), 1, 0)

  missArray<-rep(0, nrow(simulated.data))
  missArray[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==1)]<-betas[1]*as.numeric(as.character(simulated.data$X2[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==1)]))
  missArray[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==0)]<-betas[2]*as.numeric(as.character(simulated.data$X2[which(simulated.data$goodOutcomeBinary==1 & simulated.data$treatment==0)]))
  missArray[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==1)]<-betas[3]*as.numeric(as.character(simulated.data$X2[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==1)]))
  missArray[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==0)]<-betas[4]*(1-as.numeric(as.character(simulated.data$X2[which(simulated.data$goodOutcomeBinary==0 & simulated.data$treatment==0)])))

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


# multiple imputation functions
MIwithPSmatching<-function(complete.data, simulated.data, m=5, impMethod="active", intMethod="within", MIwithAux=FALSE, missingVariables, confounders, nCovariates, nBin, nNor, outcomeBetas, plot=FALSE, caliper=NULL, reverse=FALSE){

  # format imputation input
  auxColname<-"aux"
  
  ### step 1 imputation: active vs passive
  if (impMethod=="active" | impMethod=="redActive"){
    # estimate ps on observed data
    simulated.data$index<-seq(1, nrow(simulated.data), 1)
    simulated.data.partial<-simulated.data[complete.cases(simulated.data[,1:nCovariates]),]
    m_ps <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data.partial, family = binomial)
    
    simulated.data$PS<-NA
    simulated.data$PS[simulated.data.partial$index]<-m_ps$fitted
    
    # impute all missing covariates and PS at the same time
    if (MIwithAux==TRUE) { 
      impute.input<-simulated.data[,c(1:(nCovariates+2), which(colnames(simulated.data)=="PS"), which(colnames(simulated.data)==auxColname))]
    } else {
      impute.input<-simulated.data[,c(1:(nCovariates+2), which(colnames(simulated.data)=="PS"))]
    }
    #print(head(impute.input))
    imputed_Data <- mice(impute.input, m=m, maxit = 5, seed = 500, print=FALSE, remove_collinear = FALSE)
    #print(imputed_Data$pred)

    # format imputated values and output NA if MICE fails
    imp.data <- complete(imputed_Data, "long")
    if ( length(which(is.na(imp.data$PS)))>0 ) {
      print("Bypassing this iteration of simulation, empty pr values after imputation.")
      print(imputed_Data$loggedEvents)
      return(NA)
    }
    
    if (impMethod=="redActive"){
      for (i in 1:m){
        m_ps <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=imp.data[which(imp.data$.imp==i),], family = binomial)
        imp.data$PS[which(imp.data$.imp==i)]<-m_ps$fitted
      }
    }
  } else if (impMethod=="passive") {
    if (MIwithAux==TRUE) {
      impute.input<-simulated.data[,c(1:(nCovariates+2),which(colnames(simulated.data)==auxColname))]
    } else {
      impute.input<-simulated.data[,c(1:(nCovariates+2))]
    }
    #print(head(impute.input))

    #ini <- mice(impute.input, m=m, max = 0, print=FALSE)
    #pred <- ini$pred
    #for (i in 1:nrow(pred)) { pred[i,]<-1; pred[i,i]<-0 }

    if (reverse==TRUE) {
      vis<-c(which(colnames(impute.input)==missingVariables[2]), which(colnames(impute.input)==missingVariables[1]))
      names(vis)<-c(missingVariables[2], missingVariables[1])
    }
    
    imputed_Data <- mice(impute.input, m=m, maxit = 5, seed = 500, print=FALSE, remove_collinear = FALSE)
    #print(imputed_Data$pred)
    
    # format imputated values and output NA if MICE fails
    imp.data <- complete(imputed_Data, "long")
    imp.data$PS<-NA
    for (i in 1:m){
      m_ps <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=imp.data[which(imp.data$.imp==i),], family = binomial)
      imp.data$PS[which(imp.data$.imp==i)]<-m_ps$fitted
    }
  } else if (impMethod=="regPassive") {
    # initiate PS
    #simulated.data$index<-seq(1, nrow(simulated.data), 1)
    #simulated.data.partial<-simulated.data[complete.cases(simulated.data[,1:nCovariates]),]
    #m_ps <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=simulated.data.partial, family = binomial)
    
    simulated.data$PS<-NA
    #simulated.data$PS[simulated.data.partial$index]<-predict(m_ps, type = "response")
    
    # format MICE input
    impute.input<-simulated.data[,c(1:(nCovariates+2), which(colnames(simulated.data)==auxColname), which(colnames(simulated.data)=="PS"))]

 
    ini <- mice(impute.input, max = 0, print = FALSE, remove_collinear = FALSE) # initial run
    
    # change method of the variable to be passively imputed
    meth <- ini$meth
    meth["PS"] <- "~"
    
    # change pred
    pred <- ini$pred

    for (i in 1:length(missingVariables)){
      pred[missingVariables[i],]<-1
      pred[missingVariables[i], "PS"]<-0
      pred[missingVariables[i], missingVariables[i]]<-0
    }
    # estimate ps from confounders only
    pred["PS",]<-0
    pred["PS", which(colnames(pred) %in% confounders)] <- 1
    # impute non-PS variable using everything else
    pred[auxColname, ]<-1
    pred[auxColname, auxColname]<-0
    
    # change visit sequence
    #vis<-mice_regPassive(impute.input, print=FALSE, max = 0, remove_collinear = FALSE)$vis
    #vis<-vis[c(which(names(vis)==missingVariables[1]), which(names(vis)=="PS"), which(names(vis)==missingVariables[2]))]
    vis<-c(which(colnames(impute.input)==missingVariables[1]), which(colnames(impute.input)=="PS"), which(colnames(impute.input)==missingVariables[2]))
    names(vis)<-c(missingVariables[1], "PS", missingVariables[2])

    if (reverse==TRUE) {
      vis<-c(which(colnames(impute.input)==missingVariables[2]), which(colnames(impute.input)=="PS"), which(colnames(impute.input)==missingVariables[1]))
      names(vis)<-c(missingVariables[2], "PS", missingVariables[1])
    }
    # run MI again
    imputed_Data <- mice_regPassive(impute.input, method = meth, predictorMatrix = pred, vis = vis, m=m, seed=500, print = FALSE, remove_collinear = FALSE)
    
    # format output
    imp.data <- complete(imputed_Data, "long")
    
    if ( length(which(is.na(imp.data$PS)))>0 ) {
      print("Bypassing this iteration of simulation, empty pr values after imputation.")
      print(imputed_Data$loggedEvents)
      return(NA)
    }
  }

  ### step 2: PS integration
  if (intMethod=="within"){
    # for each imputated dataset, match, estimate treatment effect and variance
    mean.list <- rep(NA, m)
    variance.list <- rep(NA, m)
    stDiff1 <- list()
    stDiff2 <- list()
    stDiff3 <- list()
    stDiff4 <- list()
    oi <- list()
    
    for (i in 1:m) {
      # get each imputated dataset
      data=imp.data[imp.data$.imp == i,]
      data$index<-data$.id

      # matching
      data.matched <- Match(Y=data$outcome, Tr=data$treatment, X=data$PS, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
      
      # plot standardized difference
      stDiff1[[i]]<-plotStandardizedDiff(complete.data, data.matched, nBin=nBin, nNor=nNor, nCovariates=nCovariates, confounders=confounders, outcomeBetas=outcomeBetas, main=paste(impMethod, intMethod, MIwithAux, i, sep="-"), return=TRUE, plot=FALSE)
      stDiff2[[i]]<-plotStandardizedDiff(data[,-c(1:2)], data.matched, nBin=nBin, nNor=nNor, nCovariates=nCovariates, confounders=confounders, outcomeBetas=outcomeBetas, main=paste(impMethod, intMethod, MIwithAux, i, sep="-"), return=TRUE, plot=FALSE)
      stDiff3[[i]]<-rep(NA, length(stDiff1[[i]]))
      stDiff4[[i]]<-rep(NA, length(stDiff1[[i]]))
      # originial vs imputed values in matched samples
      oi[[i]]<-as.data.frame(table(is.na(simulated.data[data.matched$index.treated,missingVariables[1]]), is.na(simulated.data[data.matched$index.control,missingVariables[1]])))$Freq

      post_matching<-data[c(data.matched$index.treated, data.matched$index.control),]
      post_matching$pair_id<-rep(1:length(data.matched$index.treated), 2)
      dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")

      # estimate treatment effect
      mean.list[i] <- dr[1]
      # estimate variance
      variance.list[i] <- dr[2]^2
    }
    
    # calculate final mean and variance
    overall.mean<-mean(mean.list)
    overall.se<-sqrt((mean(variance.list) + var(mean.list)*(1+1/m)))
    
  } else if (intMethod=="across"){
    # average imputed PS to generated a final PS score
    impute.input2<-impute.input
    impute.input2$index<-seq(1, nrow(impute.input2), 1)
    for (i in 1:nrow(impute.input2)){
      impute.input2$finalPS[i]<-mean(imp.data$PS[which(imp.data$.id==impute.input2$index[i])])
      impute.input2$finalX2[i]<-mean(as.numeric(as.character(imp.data$X2[which(imp.data$.id==impute.input2$index[i])]))) # need to average covariates after imputation
    }
    
    # matching
    data.matched <- Match(Y=impute.input2$outcome, Tr=impute.input2$treatment, X=impute.input2$finalPS, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
    
    # plot standardized difference
    impute.input2<-impute.input2[,c("X1", "finalX2", "treatment", "outcome", "index")]
    colnames(impute.input2)[2]<-"X2"
    
    data.matched.missing<-list(index.treated=intersect(which(is.na(impute.input$X2)), data.matched$index.treated), index.control=intersect(which(is.na(impute.input$X2)), data.matched$index.control))
    data.matched.nonmiss<-list(index.treated=intersect(which(!is.na(impute.input$X2)), data.matched$index.treated), index.control=intersect(which(!is.na(impute.input$X2)), data.matched$index.control))

    stDiff1<-plotStandardizedDiff(complete.data, data.matched, nBin=nBin, nNor=nNor, nCovariates=nCovariates, confounders=confounders, outcomeBetas=outcomeBetas, main=paste(impMethod, intMethod, MIwithAux, i, sep="-"), return=TRUE, plot=FALSE)
    stDiff2<-plotStandardizedDiff(impute.input2, data.matched, nBin=nBin, nNor=nNor, nCovariates=nCovariates, confounders=confounders, outcomeBetas=outcomeBetas, main=paste(impMethod, intMethod, MIwithAux, i, sep="-"), return=TRUE, plot=FALSE)
    stDiff3<-plotStandardizedDiff(impute.input2, data.matched.nonmiss, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
    stDiff4<-plotStandardizedDiff(impute.input2, data.matched.missing, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
    
    oi<-as.data.frame(table(is.na(simulated.data[data.matched$index.treated,missingVariables[1]]), is.na(simulated.data[data.matched$index.control,missingVariables[1]])))$Freq # identical in m datasets

    post_matching<-impute.input[c(data.matched$index.treated, data.matched$index.control),]
    post_matching$pair_id<-rep(1:length(data.matched$index.treated), 2)
    dr<-ols(outcome ~ treatment, data=post_matching, cluster="pair_id")

    # estimate treatment effect
    overall.mean<-dr[1]
    
    # estimate variance
    overall.se<-dr[2]
  } else if (intMethod=="MIpar"){
    # calculate average PS coefficients
    coef<-matrix(data=NA, nrow=m, ncol=3)
    for (i in 1:m){
      m_ps <- glm(as.formula(paste("treatment ~", paste(confounders, collapse="+"))), data=imp.data[which(imp.data$.imp==i),], family = binomial)
      coef[i,]<-coef(summary(m_ps))[,1]
    }
    coef_combined<-colMeans(coef)

    # # need to average covariates after imputation
    impute.input2<-impute.input
    impute.input2$index<-seq(1, nrow(impute.input2), 1)
    for (i in 1:nrow(impute.input2)){
      impute.input2$finalX2[i]<-mean(as.numeric(as.character(imp.data$X2[which(imp.data$.id==impute.input2$index[i])]))) 
    }
    impute.input2$X1<-as.numeric(as.character(impute.input2$X1))
    
    # calculate PS
    PS<-coef_combined[1] + coef_combined[2]*impute.input2$X1 + coef_combined[3]*impute.input2$finalX2

    # matching
    data.matched <- Match(Y=impute.input2$outcome, Tr=impute.input2$treatment, X=PS, M=1, caliper=caliper, ties=FALSE, replace=FALSE) # Matching package
    
    # plot standardized difference
    impute.input2<-impute.input2[,c("X1", "finalX2", "treatment", "outcome", "index")]
    colnames(impute.input2)[2]<-"X2"
    
    data.matched.missing<-list(index.treated=intersect(which(is.na(impute.input$X2)), data.matched$index.treated), index.control=intersect(which(is.na(impute.input$X2)), data.matched$index.control))
    data.matched.nonmiss<-list(index.treated=intersect(which(!is.na(impute.input$X2)), data.matched$index.treated), index.control=intersect(which(!is.na(impute.input$X2)), data.matched$index.control))

    stDiff1<-plotStandardizedDiff(complete.data, data.matched, nBin=nBin, nNor=nNor, nCovariates=nCovariates, confounders=confounders, outcomeBetas=outcomeBetas, main=paste(impMethod, intMethod, MIwithAux, i, sep="-"), return=TRUE, plot=FALSE)
    stDiff2<-plotStandardizedDiff(impute.input2, data.matched, nBin=nBin, nNor=nNor, nCovariates=nCovariates, confounders=confounders, outcomeBetas=outcomeBetas, main=paste(impMethod, intMethod, MIwithAux, i, sep="-"), return=TRUE, plot=FALSE)
    stDiff3<-plotStandardizedDiff(impute.input2, data.matched.nonmiss, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
    stDiff4<-plotStandardizedDiff(impute.input2, data.matched.missing, nBin=nBin, nNor=nNor, nCovariates=nBin+nNor, confounders=confounders, outcomeBetas=outcomeBetas, plot=FALSE, return=TRUE)
    
    oi<-as.data.frame(table(is.na(simulated.data[data.matched$index.treated,missingVariables[1]]), is.na(simulated.data[data.matched$index.control,missingVariables[1]])))$Freq # identical in m datasets

    post_matching<-impute.input[c(data.matched$index.treated, data.matched$index.control),]
    post_matching$pair_id<-rep(1:length(data.matched$index.treated), 2)
    dr<-ols(outcome ~ treatment + X1 + X2, data=post_matching, cluster="pair_id")

    # estimate treatment effect
    overall.mean<-dr[1]
    
    # estimate variance
    overall.se<-dr[2]
  }
  
  # return values
  r<-c(overall.mean, overall.se, length(which(impute.input$treatment==1)), length(which(impute.input$treatment==0)))
  return(list(r, stDiff1, stDiff2, stDiff3, stDiff4, oi)) # treated vs untreated: oo, io, oi, ii
}
