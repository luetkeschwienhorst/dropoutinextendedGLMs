####### R formulas for the paper
####### Nonparametric estimation of mean and dispersion functions in extended generalized linear models
####### by Gibels, Prosdocimi and Claeskens
##### last modified February 2010

#### the formulas provide a tool to obtain flexible estimates for both the mean and the dispersion function
#### the method has been implemented for data coming from a normal, a poisson or a binomial 



library(splines)
library(MASS)
#library(gamlss)



############################################################################
################                Normal                      ################
############################################################################
############################################################################


#### you need to build the Bmu and Bgamma B-splines matrices before (via splines).....



## one step procedures (by one step I mean that the smoothing parameter selection is done in one step...)

normalBoth<-function(y,Bmu,Bgamma,alphaM=0,alphaG=0,lambdaM,lambdaG,maxitG=35,maxitT=1250,penordM=2,penordG=2,tol=10^-8){

# for given lambdaM and lambdaG the function normalBoth estimates the mean and dispersion function
# by performing maxitT iterations and during every iteration there are maxitG iterations allowed,
# to optimize the dispersion parameter, probably because the dispersion parameter is more difficult
# to estimate
# y			dependent variable
# Bmu 		spline matrix for mu
# Bgamma 	spline matrix for gamma
# alphaM 	coefficient for mu of mean model
# alphaG 	coefficient for gamma of dispersion model
# lambdaM 	tuning parameter of mean model
# lambdaG 	tuning parameter of dispersion model
# maxitG 	maximal number of iterations of the dispersion loop
# maxitT 	maximal number of iterations of the total loop
# penordM	order of penalties for mean model
# penordG	order of penalties for dispersion model
# tol 		tolerance which determines the convergence of the procedure
	
	# initialize alphaM and alphaG
	if (length(alphaM)==1){
		if(alphaM==0) alphaM<-rep(mean(y),l=ncol(Bmu))
		else alphaM<-rep(alphaM,l=ncol(Bmu))
	}
	if (length(alphaG)==1){
		if(alphaG==0) alphaG<-rep(0.01,l=ncol(Bgamma))
		else alphaG<-rep(alphaG,l=ncol(Bgamma))
	}

	# sample size
	n<-length(y)

	# empty list for the results	
	res<-list()

	# number of columns of Bmu and Bgamma
	k<-ncol(Bmu)
	kt<-ncol(Bgamma)
	
	# generate the penalty matrices for the mean and dispersion model
	DistMu<-diag(k)
	DistGam<-diag(kt)
	if(penordM>0){
		DistMu<-diff(DistMu,diff=penordM)
		DistMu
		}
	if(penordG>0){
		DistGam<-diff(DistGam,diff=penordG);DistGam
	}
	Pmu<-t(DistMu)%*%DistMu
	Pgam<-t(DistGam)%*%DistGam

	# counter for iterations and dummy for convergence
	nit<-0
	conv<-"false"

	# outer loop
	for(i in 1:maxitT){

		# set initial values for alphaM and alphaG
		alphaM0<-alphaM
		alphaG0<-alphaG

		# increase counter
		nit<-nit+1

		# diagonal matrix W for normal DEF where b''=1
		W<-diag(exp(c(-Bgamma%*%alphaG0)))

		# update for mean parameter based on lambdaM and the current alphaG0
		# IMPORTANT OBSERVATION: the update on alphaM does not depend on the
		# previous value for alphaM
		alphaM<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*(Pmu))%*%(t(Bmu)%*%W%*%y)

		# inner loop for dispersion parameter
		for(j in 1:maxitG){

			# q according to the paper
			q<-(Bgamma%*%alphaG0+(diag(exp(c(-Bgamma%*%alphaG0))))%*%((y-Bmu%*%alphaM)^2)-rep(1,l=length(y)))

			# update for alphaG
			alphaG<-0.5*(ginv(0.5*t(Bgamma)%*%Bgamma+(lambdaG*(Pgam))))%*%(t(Bgamma)%*%q)
			deltaeps<-sum((alphaG-alphaG0)^2)/sum(alphaG0^2)
			alphaG0<-alphaG
			if(deltaeps < tol) break
		}
		if (sum(((alphaM-alphaM0)/alphaM0)^2)<tol & sum(((alphaG-alphaG0)/alphaG0)^2)<tol){
			conv<-"true"
			break
		}
	}

	# fill the empty result list
	res$alphaM<-alphaM
	res$alphaG<-alphaG
	res$numberOfIteration<-nit
	res$conv<-conv
	W<-diag(c(exp(-Bgamma%*%alphaG)))
	res$hmean<-Bmu%*%ginv(t(Bmu)%*%W%*%Bmu+lambdaM*(Pmu))%*%t(Bmu)%*%W
	res$dfM<-dfM<-sum(diag(ginv(t(Bmu)%*%W%*%Bmu+lambdaM*(Pmu))%*%t(Bmu)%*%W%*%Bmu))
	res$hdisp<-0.5*Bgamma%*%ginv(0.5*t(Bgamma)%*%Bgamma+lambdaG*(Pgam))%*%t(Bgamma)
	res$dfG<-dfG<-sum(diag(ginv(0.5*t(Bgamma)%*%Bgamma+lambdaG*(Pgam))%*%t(Bgamma)%*%(0.5*Bgamma)))
	res$deviance<-Dev<-((y-Bmu%*%alphaM)^2)
	gammaest<-exp(Bgamma%*%alphaG)
	#### all possible sorts of 'simultaneous' information criterion
	res$GCV<-2*n*(sum(log(gammaest/Dev)+Dev/gammaest-1)+sum(Dev/gammaest))/((2*n-(dfM+dfG))^2)
	res$GCV1<-n*sum(log(gammaest/Dev)+Dev/gammaest)/((n-(dfM+dfG))^2)
	res$AIC<-sum(log(gammaest)+Dev/gammaest)+2*(res$dfM+res$dfG)
	res$AIC1<-sum(log(gammaest/Dev)+Dev/gammaest-1)+sum(Dev/gammaest)+2*(res$dfM+res$dfG)
	res
}

fisher_def_normal<-function(y,Bmu,Bgamma,alphaM=0,alphaG=0,phi=1,lambdaM,lambdaG,maxitG=35,maxitT=1250,penordM=2,penordG=2,tol=10^-8){
    
    # for given lambdaM and lambdaG the function normalBoth estimates the mean and dispersion function
    # by performing maxitT iterations and during every iteration there are maxitG iterations allowed,
    # to optimize the dispersion parameter, probably because the dispersion parameter is more difficult
    # to estimate
    # y			dependent variable
    # Bmu 		spline matrix for mu
    # Bgamma 	spline matrix for gamma
    # alphaM 	coefficient for mu of mean model
    # alphaG 	coefficient for gamma of dispersion model
    # lambdaM 	tuning parameter of mean model
    # lambdaG 	tuning parameter of dispersion model
    # maxitG 	maximal number of iterations of the dispersion loop
    # maxitT 	maximal number of iterations of the total loop
    # penordM	order of penalties for mean model
    # penordG	order of penalties for dispersion model
    # tol 		tolerance which determines the convergence of the procedure
    
    # initialize alphaM and alphaG
    if (length(alphaM)==1){
        if(alphaM==0) alphaM<-rep(mean(y),l=ncol(Bmu))
        else alphaM<-rep(alphaM,l=ncol(Bmu))
    }
    if (length(alphaG)==1){
        if(alphaG==0) alphaG<-rep(0.01,l=ncol(Bgamma))
        else alphaG<-rep(alphaG,l=ncol(Bgamma))
    }
    
    # sample size
    n<-length(y)
    
    # empty list for the results	
    res<-list()
    
    # number of columns of Bmu and Bgamma
    k<-ncol(Bmu)
    kt<-ncol(Bgamma)
    
    # generate the penalty matrices for the mean and dispersion model
    DistMu<-diag(k)
    DistGam<-diag(kt)
    if(penordM>0){
        DistMu<-diff(DistMu,diff=penordM)
        DistMu
    }
    if(penordG>0){
        DistGam<-diff(DistGam,diff=penordG);DistGam
    }
    Pmu<-t(DistMu)%*%DistMu
    Pgam<-t(DistGam)%*%DistGam
    
    # counter for iterations and dummy for convergence
    nit<-0
    conv<-"false"
    
    # outer loop
    for(i in 1:maxitT){
        
        # set initial values for alphaM and alphaG
        alphaM0<-alphaM
        alphaG0<-alphaG
        
        # increase counter
        nit<-nit+1
        
        # diagonal matrix W for normal DEF where b''=1
        W<-diag(exp(c(-Bgamma%*%alphaG0)))
        
        # update for mean parameter based on lambdaM and the current alphaG0
        # IMPORTANT OBSERVATION: the update on alphaM does not depend on the
        # previous value for alphaM
        alphaM<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*(Pmu))%*%(t(Bmu)%*%W%*%y)
        
        # inner loop for dispersion parameter
        for(j in 1:maxitG){
            
            # q according to the paper
            q<-(Bgamma%*%alphaG0+(diag(exp(c(-Bgamma%*%alphaG0))))%*%((y-Bmu%*%alphaM)^2)-rep(1,l=length(y)))
            
            # update for alphaG
            alphaG<-0.5*(ginv(0.5*t(Bgamma)%*%Bgamma+(lambdaG*(Pgam))))%*%(t(Bgamma)%*%q)
            deltaeps<-sum((alphaG-alphaG0)^2)/sum(alphaG0^2)
            alphaG0<-alphaG
            if(deltaeps < tol) break
        }
        if (sum(((alphaM-alphaM0)/alphaM0)^2)<tol & sum(((alphaG-alphaG0)/alphaG0)^2)<tol){
            conv<-"true"
            break
        }
    }
    
    # fill the empty result list
    res$alphaM<-alphaM
    res$alphaG<-alphaG
    res$meanest <- as.vector(Bmu%*%alphaM)
    res$dispest <- as.vector(phi*exp(Bgamma%*%alphaG)^{-1})
    res$numberOfIteration<-nit
    res$conv<-conv
    res
}

	

normalBGCV<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitG=35,maxitT=750,penordM=2,penordG=2,tol=10^-8){
	normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol)$GCV
}

normalBGCV1<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitG=35,maxitT=750,penordM=2,penordG=2,tol=10^-8){
	normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol)$GCV1
}

normalBAIC<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitG=35,maxitT=750,penordM=2,penordG=2,tol=10^-8){
	normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol)$AIC
}

normalBAIC1<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitG=35,maxitT=750,penordM=2,penordG=2,tol=10^-8){
	normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol)$AIC1
}



# normalBothCH<-function(y,Bmu,Bgamma,alphaM=0,alphaG=0,lambdaM,lambdaG,maxitG=35,maxitT=1250,maxitOUT=7,penordM=2,penordG=2,tol=10^-8,crit="GCV",minlM=0.0005,maxlM=8000,minlG=0.0005,maxlG=8000,tolOUT=10^-3){
# ##### this chooses 'optimal' value for the smoothing parameters and fits the model 
# ##### the optimal values for lambdaM and lambdaG are chosen simultaneously via numerical minimization
# ##### it doesn't perform very well (see paper)
# 	if (length(alphaM)==1){
# 		if(alphaM==0) alphaM<-rep(mean(y),l=ncol(Bmu))
# 		else alphaM<-rep(alphaM,l=ncol(Bmu))
# 	}
# 	if (length(alphaG)==1){
# 		if(alphaG==0) alphaG<-rep(0.01,l=ncol(Bgamma))
# 		else alphaG<-rep(alphaG,l=ncol(Bgamma))
# 	}
# 	qq<-normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdaM,lambdaG=lambdaG,maxitG=maxitG,maxitT=3,penordM=penordM,penordG=penordG,tol=tol)
# 	alphaM<-qq$alphaM
# 	alphaG<-qq$alphaG
# 	lambdaMinit<-lambdaM
# 	lambdaGinit<-lambdaG
# 	conv<-"false"
# 	res<-list()
# 	if(crit != "GCV" & crit != "AIC" & crit != "GCV1" & crit != "AIC1"){
# 		warning("unrecocognized criterion, use GCV instead")
# 		crit<-"GCV"
# 	}
# 	maxs<-c(maxlM,maxlG)
# 	mins<-c(minlM,minlG)
# 	nit<-0
# 	n<-length(y)
# 	for(i in 1:maxitOUT){
# 		nit<-nit+1
# 		lambdaM<-lambdaM+lambdaM*runif(1,-0.005,0.005)
# 		lambdaG<-lambdaG+lambdaG*runif(1,-0.005,0.005)
# 		if(crit=="AIC"){ sl<-optim(c(lambdaM,lambdaG),normalBAIC,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method= "L-BFGS-B",lower=mins,upper=maxs)
# 		if(sl$conv > 0.5){
# 			warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
# sl<-optim(c(lambdaMinit,lambdaGinit),normalBAIC,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method="L-BFGS-B",lower=mins,upper=maxs)
# }}
# if(crit=="AIC1"){ sl<-optim(c(lambdaM,lambdaG),normalBAIC1,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method= "L-BFGS-B",lower=mins,upper=maxs)
# if(sl$conv > 0.5){
# warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
# sl<-optim(c(lambdaMinit,lambdaGinit),normalBAIC1,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method="L-BFGS-B",lower=mins,upper=maxs)
# }}
# if(crit=="GCV"){ sl<-optim(c(lambdaM,lambdaG),normalBGCV,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method="L-BFGS-B",lower=mins,upper=maxs)
# if(sl$conv > 0.5){
# warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
# sl<-optim(c(lambdaMinit,lambdaGinit),normalBGCV,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method="L-BFGS-B",lower=mins,upper=maxs)
# }}
# if(crit=="GCV1"){ sl<-optim(c(lambdaM,lambdaG),normalBGCV1,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method="L-BFGS-B",lower=mins,upper=maxs)
# if(sl$conv > 0.5){
# warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
# sl<-optim(c(lambdaMinit,lambdaGinit),normalBGCV1,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,penordM=penordM,penordG=penordG,tol=tol,maxitG=maxitG,maxitT=maxitT,method="L-BFGS-B",lower=mins,upper=maxs)
# }}
# lambdaM<-sl$par[1];lambdaG<-sl$par[2]
# ee<-normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdaM,lambdaG=lambdaG,penordM=penordM,penordG=penordG,tol=tol)
# diffs<-max(abs(Bmu%*%alphaM/Bmu%*%ee$alphaM),abs(Bgamma%*%alphaG/Bgamma%*%ee$alphaG))
# # print(diffs)
# alphaM<-ee$alphaM;alphaG<-ee$alphaG
# if(diffs<=(1+tolOUT) & diffs>=(1 - tolOUT)){
# conv<-"true"
# break
# }}
# res$conv<-conv;res$numberOfIteration<-nit
# res$alphaM<-ee$alphaM;res$alphaG<-ee$alphaG
# res$lambdaM<-lambdaM;res$lambdaG<-lambdaG;res$vals<-c(lambdaM,lambdaG)
# W<-diag(c(exp(-Bgamma%*%ee$alphaG)))
# res$hmean<-ee$hmean;res$hdisp<-ee$hdisp
# res$dfM<-ee$dfM;res$dfG<-ee$dfG
# res$deviance<-Dev<-((y-Bmu%*%alphaM)^2);gammaest<-exp(Bgamma%*%alphaG)
# res$GCV<-2*n*(sum(log(gammaest/Dev)+Dev/gammaest-1)+sum(Dev/gammaest))/((2*n-(res$dfM+res$dfG))^2)
# res$AIC<-sum(log(gammaest)+Dev/gammaest)+2*(res$dfM+res$dfG)
# res
# }






############################
#### two steps procedures

NormalMean<-function(y,Bmu,lambdaM,alphaM=0,gammatrue,penordM=2,tol=10^-10){
# for a given lambdaM and a given dispersion function the function NormalMean
# performs one Fisher iteration for alphaM
    
    # sample size and empty list for results
	n<-length(y)
	rr<-list()
	
	# penalty matrix Pmu according to penordM
	DistMu<-diag(ncol(Bmu))
	if (penordM>0) DistMu<-diff(DistMu,diff=penordM)
	Pmu<-t(DistMu)%*%DistMu
	
	# W for given gammatrue
	W<-diag(as.vector(1/gammatrue))

	# Fisher update for alphaM
	matr<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
	alfa1<-matr%*%y

	# fill the empty result list
	rr$alphaM<-alfa1
	df<-sum(diag(matr%*%Bmu))
	rr$df<-df
	## information criteria
	rr$GCV<-sum(((y-Bmu%*%rr$alphaM)^2)/gammatrue)/((1-df/n)^2)
	rr$AIC<-sum(((y-Bmu%*%rr$alphaM)^2)/gammatrue)+2*df
	rr
}



NormalDisp<-function(y,Bgamma,alphaG=0,lambdaG,truemean,penordG=2,tol=10^-8,maxit=25){
# for a given lambdaG and a given mean function, the function NormalDisp
# performs up to 25 Fisher iterations for alphaG

	# initialize alphaG incase is has not been done yet
	if (length(alphaG)==1){
		if(alphaG==0) alphaG<-rep(0.01,l=ncol(Bgamma))
		else alphaG<-rep(alphaG,l=ncol(Bgamma))
	}

	# sample size and empty list for results
	n<-length(y)
	rr<-list()

	# penalty matrix Pgam according to penordG
	DistGam<-diag(ncol(Bgamma))
	if(penordG>0) DistGam<-diff(DistGam,diff=penordG)
	Pgam<-t(DistGam)%*%DistGam

	# hatmatrix, initial value for alphaG and eps, counter
	matr<-0.5*(ginv(0.5*t(Bgamma)%*%Bgamma+lambdaG*(Pgam)))%*%t(Bgamma)
	alphaG0<-alphaG
	eps<-1
	ind<-0

	# Fisher updates for alphaG
	while(eps > tol){
		ind<-ind+1
		q<-(Bgamma%*%alphaG0+(diag(exp(c(-Bgamma%*%alphaG0))))%*%((y-truemean)^2)-rep(1,l=length(y)))
		alphaG<-matr%*%q
		eps<-sum((alphaG-alphaG0)^2)/sum(alphaG0^2)
		alphaG0<-alphaG
		if (is.na(eps)) break
		if (ind>maxit) break
	}

	# fill the empty result list
	hdisp<-Bgamma%*%matr
	rr$df<-df<-sum(diag(matr%*%Bgamma))
	rr$alphaG<-as.vector(alphaG)
	rr$hdisp<-hdisp
	gammaest<-exp(Bgamma%*%alphaG)
	Dev<-(y-truemean)^2
	## information criteria
	rr$GCV<-n*(sum(log(gammaest/Dev)+Dev/gammaest-1)/((1-df/n)^2))
	rr$AIC<-sum(log(gammaest)+Dev/gammaest)+(2*df)
	rr
}




NormalMeanGCV<-function(y,Bmu,lambdasM,gammatrue,penordM=2,tol=10^-10){
# for all values in lambdasM the function NormalMeanGCV performs a Fisher
# iteration on alphaM and then computes the GCV, given the true dispersion function
	
	# initialization and empty list for results
	GCV<-NULL
	AIC<-NULL
	n<-length(y)
	co<-NULL
	lambdaseq<-NULL
	rr<-list()

	# penalty matrix Pmu according to penordM
	DistMu<-diag(ncol(Bmu))
	if (penordM>0) DistMu<-diff(DistMu,diff=penordM)
	Pmu<-t(DistMu)%*%DistMu

	# weight matrix W according to gammatrue
	W<-diag(as.vector(1/gammatrue))

	# Fisher updates for all values in lambdasM
	for( i in 1:length(lambdasM)){
		lambdaM<-lambdasM[i]
		matr<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
		alfa1<-matr%*%y
		#hmean<-Bmu%*%matr
		df<-sum(diag(matr%*%Bmu))
		GCV[i]<-sum(((gammatrue^-1)*(y-Bmu%*%alfa1)^2)/((1-df/n)^2))
		AIC[i]<-sum(((y-Bmu%*%alfa1)^2)/gammatrue)+2*df
		lambdaseq<-c(lambdaseq,lambdaM)
		#rm(alfa1)
	}
	rr$GCV<-GCV
	rr$AIC<-AIC
	rr$lambdaseq<-lambdaseq
	rr
}


nmAIC<-function(lambdaM,y,Bmu,alphaM=0,gammatrue,penordM=2,tol=10^-10){
	rr<-NormalMean(y=y,Bmu=Bmu,lambdaM=lambdaM,alphaM=alphaM,gammatrue=gammatrue,penordM=penordM,tol=tol)$AIC;rr
}
nmGCV<-function(lambdaM,y,Bmu,alphaM=0,gammatrue,penordM=2,tol=10^-10){
	rr<-NormalMean(y=y,Bmu=Bmu,lambdaM=lambdaM,alphaM=alphaM,gammatrue=gammatrue,penordM=penordM,tol=tol)$GCV;rr
}

NormalDispGCV<-function(y,Bgamma,alphaGinit=0,lambdasG,truemean,penordG=2,tol=10^-10,maxit=25){
# for all values in lambdasG the function NormalDispGCV calls the function NormalDisp, i.e. it performs up to
# 25 Fisher iterations on alphaGinit for some given truemean, and then computes the GCV
	
	# initialization
	GCV<-AIC<-NULL
	n<-length(y)
	lambdaseq<-NULL
	rr<-list()

	# call NormalDisp for all values in lambdasG
	for( i in 1:length(lambdasG)){
		lambdaG<-lambdasG[i]
		dd<-NormalDisp(y=y,Bgamma=Bgamma,alphaG=alphaGinit,lambdaG=lambdaG,truemean=truemean,penordG=penordG,tol=tol,maxit=maxit)
		GCV[i]<-dd$GCV
		AIC[i]<-dd$AIC
		lambdaseq<-c(lambdaseq,lambdaG)
	}
	rr$GCV<-GCV
	rr$AIC<-AIC
	rr$lambdaseq<-lambdaseq
	rr
}

nvGCV<-function(lambdaG,y,Bgamma,alphaG=0,truemean,penordG=2,tol=10^-8,maxit=25){
rr<-NormalDisp(y=y,Bgamma=Bgamma,alphaG=alphaG,lambdaG=lambdaG,truemean=truemean,penordG=penordG,tol=tol,maxit=maxit)$GCV;rr
}
nvAIC<-function(lambdaG,y,Bgamma,alphaG=0,truemean,penordG=2,tol=10^-8,maxit=25){
rr<-NormalDisp(y=y,Bgamma=Bgamma,alphaG=alphaG,lambdaG=lambdaG,truemean=truemean,penordG=penordG,tol=tol,maxit=maxit)$AIC;rr
}




NormalBothGCV<-function(y,Bmu,Bgamma,seqM,seqG,lenM,lenG,alphaM=0,alphaG=0,penordM=2,
	penordG=2,maxitG=45,tol=10^-10,tolSTOP=10^-2,plotM=FALSE,plotG=FALSE,minlM=0.005,
	maxlM=8000,minlG=0.005,maxlG=8000,uppars=TRUE,crit="GCV"){
##### this chooses 'optimal' value for the smoothing parameters in two steps and fits
##### the model the smoothing parameter choice is done via a grid search

# y        	dependent variable
# Bmu      	spline matrix for mu
# Bgamma	spline matrix for gamma
# seqM 		min and max of lambdas for the mean function
# seqG 		min and max of lambdas for the dispersion function
# lenM		number of grid points between min and max of lambdas for mean function
# lenG		number of grid points between min and max of lambdas for dispersion function

	# initialize list of lambdas for mu and gamma
	llM<-seq(min(seqM),max(seqM),length=lenM)
	llG<-seq(min(seqG),max(seqG),length=lenG)
	diffsMean<-NULL
	diffsDisp<-NULL

	# penalty matrix for mu
	n<-length(y)
	DistMu<-diag(ncol(Bmu))
	if (penordM>0) DistMu<-diff(DistMu,diff=penordM)
	Pmu<-t(DistMu)%*%DistMu

	# initialize alphaM and alphaG
	if (length(alphaM)==1){
		if(alphaM==0) alphaM<-rep(mean(y),l=ncol(Bmu))
		else alphaM<-rep(alphaM,l=ncol(Bmu))
	}
	if (length(alphaG)==1){
		if(alphaG==0) alphaG<-rep(0.01,l=ncol(Bgamma))
		else alphaG<-rep(alphaG,l=ncol(Bgamma))
	}

	# get rough initial estimate of mean and dispersion function, by performing
	# 3 iterations
	qq<-normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,
		lambdaM=mean(llM),lambdaG=mean(llG),maxitG=maxitG,maxitT=3,penordM=penordM,
		penordG=penordG,tol=tol)
	
	# use initial estimate to set initial values for alphaM, alphaG, dispest and meanest
	alphaM<-qq$alphaM
	alphaG<-qq$alphaG
	dispest<-as.vector(exp(Bgamma%*%alphaG))
	meanest<-as.vector(Bmu%*%alphaM)

	# initial values to start from for every iteration of the loop
	alphaMinit<-alphaM
	alphaGinit<-alphaG
	
	# empty lists to fill
	gcvMmins<-NULL
	gcvGmins<-NULL
	res<-NULL
	ntot<-0
	lambdaM<-NULL
	lambdaG<-NULL
	extremesM<-NULL
	extremesG<-NULL
	numsM<-NULL
	numsG<-NULL

	# outer loop
	for(j in 2:8){
		#print(paste("Outer loop: ",sep="",j))
		alphaM<-alphaMinit
		alphaG<-alphaGinit
		ntot<-ntot+1
		checkM<-NULL
		checkM[1]<-10^-5
		checkG<-NULL
		checkG[1]<-10^-5
		for(i in 1:5){
		    #print(paste("Mean loop: ",sep="",i))
		    
			# get GCVs for all values in llM after performing one Fisher iteration
			gcv<-NormalMeanGCV(y=y,Bmu=Bmu,lambdasM=llM,gammatrue=dispest,
				penordM=2,tol=10^-10)
			if(crit=="AIC"){gcv$GCV<-gcv$AIC}
			if (plotM==TRUE) plot(llM,gcv$GCV,type="l")

			# get index of minimal gcv value and select the corresponding lambda
			ss<-min(seq(1,lenM)[gcv$GCV==min(gcv$GCV)])
			lM<-llM[ss]
			checkM[i+1]<-lM

			# check if the new minimum and the previous minimum are sufficiently close
			# to each other, if yes, then start new grid on the interval [lM*0.2,lM*5]
			# with length lenM
			if(checkM[i]/checkM[i+1]<1.05 & checkM[i]/checkM[i+1]>0.95){
				print("GCV-breaking criterion fulfilled.")
				llM<-sort(seq(lM*0.2,lM*5,length=lenM))
				# save the index of lambda with minimal GCV
				extremesM<-c(extremesM,as.character(ss))
				break
			# if no, then
			}else{
				paceM<-llM[2]-llM[1]
				#lM<-llM[ss]
				if(ss==lenM){
					print("Check flatness of GCV curve of mean.")
					### this checks whether the GCV curve is flat
					### if the curve is flat we don't take the maximum lambda
					### value as chosen lambda....
					diffGCV1<-NULL
					diffGCV2<-NULL
					diffGCV3<-NULL
					diffGCV4<-NULL
					diffGCV5<-NULL
					diffsTOT<-NULL
					diffGCVFin<-NULL
					eps<-min(gcv$GCV)*10^-3
					sum5<-NULL
					for( k in 1:(lenM-5)){
						diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
						diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
						diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
						diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
						diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
						diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenM])<eps)
						diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],
							diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))
					}
					ns<-min(ss,seq(1,length(llM))[diffsTOT>5.5])
					if(ns != ss){
						lM<-llM[ns]
						extremesM<-c(extremesM,"max - cutted")
						llM<-sort(seq(from=max(minlM,min(llM[1]*0.8,lM-6*paceM)),by=1.8*paceM,l=lenM))
						break
					}
					if(ns==ss){
						lM<-llM[ss]
						llM<- sort(seq(from=max(minlM,lM-10*paceM),by=1.2*paceM,l=lenM))
						extremesM<-c(extremesM,"max")
					}
				}
				if(ss==1){
					if(lM==minlM){
						llM<-sort(seq(minlM,by=1.2*paceM,l=lenM))
					}else{
						llM<-sort(seq(max(minlM,lM-(0.8*lenM*paceM)),by=1.2*paceM,l=lenM))
					}
					extremesM<-c(extremesM,"min")
				}
				if(ss!=1 & ss!=lenM){
					llM<-sort(seq(from=max(minlM,lM-6*paceM),to=max(minlM+10*paceM,lM+6*paceM),l=lenM))
					extremesM<-c(extremesM,as.character(ss))
				}
			}
		}
		extremesM<-c(extremesM," ")
		if(lM>maxlM){
			lM<-maxlM
			warning("lM higher then maxlM forced to be maxlM")
		}
		W<-diag(1/dispest)
		matr<-ginv(t(Bmu)%*%W%*%Bmu+lM*Pmu)%*%t(Bmu)%*%W
		alphaM1<-matr%*%y
		alphaMdiff<-sum(((alphaM1-alphaM)/(alphaM))^2)
		alphaM<-alphaM1
		meanest<-as.vector(Bmu%*%alphaM)
		rm(gcv,ss)
		numsM<-c(numsM,i)
		for(i in 1:5){
		    #print(paste("Dispersion loop: ",sep="",i))
			gcv<-NormalDispGCV(y=y,Bgamma=Bgamma,alphaGinit=alphaGinit,lambdasG=llG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG)
			if(crit=="AIC") {gcv$GCV<-gcv$AIC}
			gcvG<-min(gcv$GCV)	
			if(plotG==TRUE) plot(llG,gcv$GCV,type="l")
			ss<-min(seq(1,lenG)[gcv$GCV==min(gcv$GCV)])
			lG<-llG[ss]
			checkG[i+1]<-lG#min(gcv$GCV)
			if(checkG[i]/checkG[i+1]<1.05 & checkG[i]/checkG[i+1]>.95){
				llG<-sort(seq(0.5*lG,lG*2,length=lenG))
				extremesG<-c(extremesG,as.character(ss))
				break
			}else{
				paceG<-llG[2]-llG[1]
				if(ss==lenG){
					print("Check flatness of GCV curve of dispersion!")
					diffGCV1<-NULL
					diffGCV2<-NULL
					diffGCV3<-NULL
					diffGCV4<-NULL
					diffGCV5<-NULL
					diffsTOT<-NULL
					diffGCVFin<-NULL
					eps<-min(gcv$GCV)*10^-2
					for( k in 1:(lenG-5)){
						diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
						diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
						diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
						diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
						diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
						diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenG])<eps)
						diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))
					}
					sum6<-diffsTOT
					ns<-min(ss,seq(1,length(llG))[sum6>5.5])
					rm(sum6)
					if(ns != ss){
						lG<-llG[ns]
						extremesG<-c(extremesG,"max - cutted")
						llG<-sort(seq(from=max(minlG,min(llG[1]*0.8,lG-6*paceG)),by=1.8*paceG,l=lenG))
						break
					}
					if(ns==ss){
						lG<-llG[ss]
						llG<- sort(seq(from=max(minlG,lG-10*paceG),by=1.2*paceG,l=lenG))
						extremesG<-c(extremesG,"max")
					}
				}
				if(ss==1){
					if(lG==0){
						llG<-sort(seq(max(minlG,lG-(0.8*lenG*paceG)),by=1.2*paceG,l=lenG))
					}else{
						llG<-sort(seq(from=max(minlG,lG-(0.8*lenG*paceG)),by=1.2*paceG,l=lenG))
					}
					extremesG<-c(extremesG,"min")
				}
				if(ss!=1 & ss!=lenG){
					llG<-sort(seq(from=max(minlG,lG-6*paceG),to=max(minlG+10*paceG,lG+6*paceG),l=lenG))
					extremesG<-c(extremesG,as.character(ss))
				}
			}
		}
		extremesG<-c(extremesG," ")
		if(lG>maxlG){
			lG<-maxlG
			warning("lG higher then maxlG forced to be maxlG")
		}
		numsG<-c(numsG,i)
		zz<-NormalDisp(y=y,Bgamma=Bgamma,alphaG=alphaG,lambdaG=lG,truemean=meanest,penordG=penordG,tol=tol)
		dispest<-as.vector(exp(Bgamma%*%zz$alphaG))
		alphaGdiff<-sum(((zz$alphaG-alphaG)/(alphaG))^2)
		alphaG<-zz$alphaG
		rm(gcv,ss)
		lambdaM<-c(lambdaM,lM)
		lambdaG<-c(lambdaG,lG)
		if(uppars == TRUE) alphaGinit<-alphaG
		if(max(alphaGdiff,alphaMdiff)<tolSTOP){
			break
		}
	}
	res$vals<-c(lM,lG)
	res$ntri<-j-1
	res$sequencesM<-extremesM
	res$sequencesG<-extremesG
	res$meanest<-as.vector(meanest)
	res$dispest<-as.vector(dispest)
	res$deviance<-(y-meanest)^2
	res$lambdasM<-lambdaM
	res$lambdasG<-lambdaG
	res$alphaM<-as.vector(alphaM)
	res$alphaG<-as.vector(zz$alphaG)
	#res$hmean<-Bmu%*%matr;res$hdisp<-zz$hdisp
	res$dfM<-sum(diag(matr%*%Bmu))
	res$dfG<-zz$df
	res$vecNumM<-numsM
	res$vecNumG<-numsG
	res
	}



Normal2sNUM<-function(y,Bmu,Bgamma,lambdaM,lambdaG,alphaM=0,alphaG=0,penordM=2,penordG=2,tol=10^-10,tolSTOP=10^-2,minlM=0.005,maxlM=8000,minlG=0.005,maxlG=8000,uppars=TRUE,crit="GCV",maxitT=120,maxitG=45){
##### this chooses 'optimal' value for the smoothing parametersin two steps and fits the model 
##### the smoothing parameter choice is done via numerical minimization of GCV or AIC
if(crit != "GCV" & crit != "AIC"){
warning("unrecognized criterion, use GCV instead")
crit<-"GCV"
}
n<-length(y)
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(mean(y),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.01,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bgamma))
}
# qq<-normalBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdaM,lambdaG=lambdaG,maxitG=maxitG,maxitT=3,penordM=penordM,penordG=penordG,tol=tol)
# alphaM<-qq$alphaM;alphaG<-qq$alphaG
alphaMinit<-alphaM;alphaGinit<-alphaG
lambdaMinit<-lambdaM;lambdaGinit<-lambdaG
dispest<-as.vector(exp(Bgamma%*%alphaG));meanest<-as.vector(Bmu%*%alphaM)
res<-list();ntot<-0
for(j in 2:8){
alphaM<-alphaMinit;alphaG<-alphaGinit
ntot<-ntot+1
if(crit=="AIC"){
sl<-optim(lambdaM,nmAIC,y=y,Bmu=Bmu,alphaM=alphaM,gammatrue=dispest,penordM=penordM,tol=tol,method="L-BFGS-B",lower=minlM,upper=maxlM)
if(sl$conv > 0.5){
warning(paste("at iteration",j,"MEAN optim did not converge, I restart from initial lambda"))
sl<-optim(lambdaMinit,nmAIC,y=y,Bmu=Bmu,alphaM=alphaM,gammatrue=dispest,penordM=penordM,tol=tol,method="L-BFGS-B",lower=minlM,upper=maxlM)
}
if(sl$par <=(lambdaM+0.00001*lambdaM) & sl$par >=(lambdaM-0.00001*lambdaM)  & sl$par != minlM  & sl$par != maxlM){
warning(paste("at iteration",j,"MEAN sp did not change, I change it a bit"))
if(sl$par >25) ll<-lambdaM/runif(1,14,22)
if(sl$par <25) ll<-lambdaM*runif(1,4,8)
sl<-optim(ll,nmAIC,y=y,Bmu=Bmu,alphaM=alphaM,gammatrue=dispest,penordM=penordM,tol=tol,method="L-BFGS-B",lower=minlM,upper=maxlM)
}}
if(crit=="GCV"){
sl<-optim(lambdaM,nmGCV,y=y,Bmu=Bmu,alphaM=alphaM,gammatrue=dispest,penordM=penordM,tol=tol,method="L-BFGS-B",lower=minlM,upper=maxlM)
if(sl$conv > 0.5){
warning(paste("at iteration",j,"MEAN optim did not converge, I restart from initial lambda"))
sl<-optim(lambdaMinit,nmGCV,y=y,Bmu=Bmu,alphaM=alphaM,gammatrue=dispest,penordM=penordM,tol=tol,method="L-BFGS-B",lower=minlM,upper=maxlM)
}
if(sl$par <=lambdaM-0.00001*lambdaM & lambdaM-0.00001*lambdaM  & sl$par != minlM  & sl$par != maxlM){
warning(paste("at iteration",j,"MEAN sp did not change, I change it a bit"))
if(sl$par >25) ll<-lambdaM/runif(1,14,22)
if(sl$par <25) ll<-lambdaM*runif(1,4,8)
sl<-optim(ll,nmGCV,y=y,Bmu=Bmu,alphaM=alphaM,gammatrue=dispest,penordM=penordM,tol=tol,method="L-BFGS-B",lower=minlM,upper=maxlM)
}}
lambdaM<-sl$par
mm<-NormalMean(y=y,Bmu=Bmu,lambdaM=lambdaM,gammatrue=dispest,penordM=penordM,tol=tol)
alphaM<-mm$alphaM;meanest<-Bmu%*%alphaM
if(crit=="AIC"){
sl<-optim(lambdaG,nvAIC,y=y,Bgamma=Bgamma,alphaG=alphaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG,method="L-BFGS-B",lower=minlG,upper=maxlG)
if(sl$conv > 0.5){
warning(paste("at iteration",j,"DISP optim did not converge, I restart from initial lambda"))
sl<-optim(lambdaGinit,nvAIC,y=y,Bgamma=Bgamma,alphaG=alphaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG,method="L-BFGS-B",lower=minlG,upper=maxlG)
}
if(sl$par <=lambdaG-0.00001*lambdaG & lambdaG-0.00001*lambdaG & sl$par != minlG  & sl$par != maxlG){
warning(paste("at iteration",j,"DISP sp did not change, I change it a bit"))
if(sl$par >25) ll<-lambdaG/runif(1,14,22)
if(sl$par <25) ll<-lambdaG*runif(1,4,8)
sl<-optim(ll,nvAIC,y=y,Bgamma=Bgamma,alphaG=alphaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG,method="L-BFGS-B",lower=minlG,upper=maxlG)
}}
if(crit=="GCV"){
sl<-optim(lambdaG,nvGCV,y=y,Bgamma=Bgamma,alphaG=alphaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG,method="L-BFGS-B",lower=minlG,upper=maxlG)
if(sl$conv > 0.5){
warning(paste("at iteration",j,"DISP optim did not converge, I restart from initial lambda"))
sl<-optim(lambdaGinit,nvGCV,y=y,Bgamma=Bgamma,alphaG=alphaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG,method="L-BFGS-B",lower=minlG,upper=maxlG)
}
if(sl$par <=lambdaG-0.00001*lambdaG & lambdaG-0.00001*lambdaG & sl$par != minlG  & sl$par != maxlG){
warning(paste("at iteration",j,"DISP sp did not change, I change it a bit"))
if(sl$par >25) ll<-lambdaG/runif(1,14,22)
if(sl$par <25) ll<-lambdaG*runif(1,4,8)
sl<-optim(ll,nvGCV,y=y,Bgamma=Bgamma,alphaG=alphaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG,method="L-BFGS-B",lower=minlG,upper=maxlG)
}}
lambdaG<-sl$par
vv<-NormalDisp(y=y,Bgamma=Bgamma,alphaG=alphaG,lambdaG=lambdaG,truemean=meanest,penordG=penordG,tol=tol,maxit=maxitG)
alphaG<-vv$alphaG;dispest<-exp(Bgamma%*%alphaG)
alphaGdiff<-sum(((alphaG-alphaGinit)/(alphaGinit))^2)
alphaMdiff<-sum(((alphaM-alphaMinit)/(alphaMinit))^2)
if(uppars == TRUE) {alphaMinit<-alphaM;alphaGinit<-alphaG}
if(max(alphaGdiff,alphaMdiff)<tolSTOP){
break
}}
res$vals<-c(lambdaM,lambdaG)
res$ntri<-j-1
res$meanest<-as.vector(meanest);res$dispest<-as.vector(dispest)
res$deviance<-(y-meanest)^2
res$alphaM<-as.vector(alphaM);res$alphaG<-as.vector(alphaG)
res$dfM<-mm$df;res$dfG<-vv$df
res
}


# lidar example
# library(SemiPar)
# data(lidar)
# x<-lidar$range
# y<-lidar$logratio
# dx<-(max(x)-min(x))/43
# kn<-seq(min(x)-3*dx,max(x)+3*dx,by=dx)
# basem<-splineDesign(knots=kn,x=x,ord=(3+1),0*x,outer=TRUE)
# dx<-(max(x)-min(x))/33
# kn<-seq(min(x)-3*dx,max(x)+3*dx,by=dx)
# based<-splineDesign(knots=kn,x=x,ord=(3+1),0*x,outer=TRUE)
# ll<-NormalBothGCV(y=y,Bmu=basem,Bgamma=based,seqM=c(1,2545),seqG=c(0.1,220),lenM=57,lenG=63,minlM=0.000005,maxlM=12500,minlG=0.0001,maxlG=8000,plotM=TRUE,plotG=TRUE)
# par(mfrow=c(1,2))
# plot(x,y,ylab="Logratio",xlab="Range",pch=20,cex=0.7,bty="l",main="(a)")
# points(x,ll$meanest,col=4,type="l",lwd=2)
# plot(x,(y-ll$meanest)^2,xlab="Range",ylab="Squared residuals",pch=20,cex=0.7,bty="l",main="(b)")
# points(x,ll$dispest,col=4,type="l",lwd=2)










############################################################################
################                Poisson                     ################
############################################################################
############################################################################

### ranPois generates samples from the double poisson
ranPois<-function(size=1,lambda,theta=1){
rho<-1/theta
pp<-NULL;ss<-NULL;
lambda<-rep(lambda,l=size)
rho<-rep(rho,l=size)
for(i in 1:size){
ll<-lambda[i];tt<-rho[i]
for(k in 1:141){
y<-(k-1)
pp[k]<-(tt^(1/2))*(exp(-y)*(y^y)/factorial(y))*(exp(1)*ll/y)^(tt*y)
probpos<-length(pp[pp!=Inf])
}
ss[i]<-sample(seq(0,(probpos-1)),size=1,p=pp[1:(probpos)])
}
ss
}

### the deviance function for poisson data
DevPois<-function(y,mm){
Dev<-2*as.vector(y*log(y)-y-y*log(mm)+mm)
ss<-seq(1,length(y))[y==0];Dev[ss]<-2*mm[ss]
Dev
}

#### notes for most functions 

#### you need to build the Bmu and Bgamma B-splines matrices before (via splines).....

##### type can be either "quasi" or "pseudo"
  ## with quasi the procedure described in the paper is performed: use the deviance residuals as a responce for the dispersion model
  ## with pseudo a pseudo likelihood procedure is performed (Davidian and Carrol 1986): use the pearson residuals as a responce for the dispersion model  
  
##### link can be either "nb" or "exp" and it corresponds to the link function for the dispersion function
  ## with "exp" we have gamma=exp(xi), the default in the paper
  ## with "NB"/"nb" we have gamma=1+exp(xi), as in a negative binomial (note how gamma can not be larger then 1 - only overdispersed data)



############################
#### one steps procedures


poissonBoth1step<-function(y,Bmu,Bgamma,alphaM=0,alphaG=0,lambdaM,lambdaG,maxitM=95,maxitG=95,maxitT=130,penordM=2,penordG=2,tol=10^-10,type="quasi",link="exp"){
##### for given lambdaM and lambdaG estimate both functions
if(link == "nb") link<-"NB"
if(link != "exp" & link != "NB"){ warning("which link? I use exp");link<-"exp"}
n<-length(y);res<-list();k<-ncol(Bmu);kt<-ncol(Bgamma)
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
if (length(alphaG)==1){	
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bgamma))
}
DistMu<-diag(k);DistGam<-diag(kt);nit<-0
if(penordM>0) DistMu<-diff(DistMu,diff=penordM)
if(penordG>0) DistGam<-diff(DistGam,diff=penordG)
Pmu<-t(DistMu)%*%DistMu;Pgam<-t(DistGam)%*%DistGam
Xmat<-matrix(0,ncol=(k+kt),nrow=(2*n))
Xmat[1:n,1:k]<-Bmu;Xmat[(1+n):(2*n),(1+k):(k+kt)]<-Bgamma
PenMat<-matrix(0,ncol=(k+kt),nrow=(k+kt))
PenMat[1:k,1:k]<-lambdaM*Pmu;PenMat[(1+k):(k+kt),(1+k):(k+kt)]<-lambdaG*Pgam
for(l in 1:10){
W<-diag(as.vector(exp(-Bgamma%*%alphaG)*exp(Bmu%*%alphaM)))
alphaM<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%(Bmu%*%alphaM+(y-exp(Bmu%*%alphaM))/c(exp(Bmu%*%alphaM)))
}
### this is to be sure we have a reasonable starting point (of course initial lambdaM needs to be reasonable...)
alphas<-c(alphaM,alphaG)
for(i in 1:maxitT){    ### here the iteration begins
nit<-nit+1;truemean<-exp(Bmu%*%alphaM)
if (type=="quasi"){
Dev<-DevPois(y,truemean)}
if(type=="pseudo"){
Dev<-((y-truemean)^2)/(truemean)}
if(type != "pseudo" & type != "quasi" ){
stop("the type of likelihood to use is not recognized")}
zvec<-c(Bmu%*%alphaM+(y-exp(Bmu%*%alphaM))*exp(-Bmu%*%alphaM));W<-diag(as.vector(exp(-Bgamma%*%alphaG)*exp(Bmu%*%alphaM)))
alphaM1<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%zvec
if(link=="exp"){
W<-diag(x=.5,length(y));zvec<-Bgamma%*%alphaG+((Dev-c(exp(Bgamma%*%alphaG)))/exp(Bgamma%*%alphaG))}
if(link=="NB"){
W<-diag(0.5*as.vector((exp(Bgamma%*%alphaG)/(1+exp(Bgamma%*%alphaG)))^2))
zvec<-Bgamma%*%alphaG+((Dev-1-exp(Bgamma%*%alphaG))/exp(Bgamma%*%alphaG))}
alphaG1<-ginv(t(Bgamma)%*%W%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%W%*%zvec
alphaGeps<-sum(((alphaG1-alphaG)/alphaG1)^2)
alphas1<-c(alphaM1,alphaG1)
alphaeps<-sum(((alphas1-alphas)/alphas)^2)
alphas<-alphas1;alphaM<-alphas[1:k];alphaG<-alphas[(k+1):(k+kt)] #alphaGeps<-10^{-j*3}
if (alphaeps <= tol) break
if (i == maxitT){
if(i != 1) warning("Iteration limit reached without convergence")}
}
res$meanest<-meanest<-exp(Bmu%*%alphaM)
W<-diag(as.vector(exp(-Bgamma%*%alphaG)*exp(Bmu%*%alphaM)))
dfM<-sum(diag(ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%Bmu));res$dfM<-dfM
res$deviance<-Dev<-DevPois(y,truemean)
if(link=="exp") {W<-diag(x=.5,length(y));res$gammaest<-gammaest<-exp(Bgamma%*%alphaG)}
if(link=="NB") {W<-diag(0.5*as.vector((exp(Bgamma%*%alphaG)/(1+exp(Bgamma%*%alphaG)))^2))
res$gammaest<-gammaest<-(1+exp(Bgamma%*%alphaG))}
dfG<-sum(diag(ginv(t(Bgamma)%*%W%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%W%*%Bgamma));res$dfG<-dfG
res$alphaM<-as.vector(alphaM);res$alphaG<-as.vector(alphaG)
res$numberOfIteration<-nit
## differnt information criteria
res$AIC<-sum(log(gammaest)+Dev/gammaest)+2*(dfM+dfG)
res$AIC3<-sum(log(gammaest/Dev)+Dev/gammaest-1)+sum(Dev/gammaest)+2*(dfM+dfG)
res$GCV1<-n*n*sum(log(gammaest)+Dev/gammaest-log(Dev))/((n-(dfM+dfG))^2)
res$GCV3<-2*n*n*(sum(log(gammaest/Dev)+Dev/gammaest-1)+sum(Dev/gammaest))/((2*n-(dfM+dfG))^2)
res
}

fisher_def_poisson<-function(y,Bmu,Bgamma,alphaM=0,alphaG=0,lambdaM,lambdaG,maxitM=95,maxitG=95,maxitT=130,penordM=2,penordG=2,tol=10^-10,type="quasi",link="exp"){
    ##### for given lambdaM and lambdaG estimate both functions
    conv <- FALSE
    if(link == "nb") link<-"NB"
    if(link != "exp" & link != "NB"){ warning("which link? I use exp");link<-"exp"}
    n<-length(y);res<-list();k<-ncol(Bmu);kt<-ncol(Bgamma)
    if (length(alphaM)==1){
        if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
        else alphaM<-rep(alphaM,l=ncol(Bmu))
    }
    if (length(alphaG)==1){	
        if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
        else alphaG<-rep(alphaG,l=ncol(Bgamma))
    }
    DistMu<-diag(k);DistGam<-diag(kt);nit<-0
    if(penordM>0) DistMu<-diff(DistMu,diff=penordM)
    if(penordG>0) DistGam<-diff(DistGam,diff=penordG)
    Pmu<-t(DistMu)%*%DistMu;Pgam<-t(DistGam)%*%DistGam
    Xmat<-matrix(0,ncol=(k+kt),nrow=(2*n))
    Xmat[1:n,1:k]<-Bmu;Xmat[(1+n):(2*n),(1+k):(k+kt)]<-Bgamma
    PenMat<-matrix(0,ncol=(k+kt),nrow=(k+kt))
    PenMat[1:k,1:k]<-lambdaM*Pmu;PenMat[(1+k):(k+kt),(1+k):(k+kt)]<-lambdaG*Pgam
    for(l in 1:10){
        W<-diag(as.vector(exp(-Bgamma%*%alphaG)*exp(Bmu%*%alphaM)))
        alphaM<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%(Bmu%*%alphaM+(y-exp(Bmu%*%alphaM))/c(exp(Bmu%*%alphaM)))
    }
    ### this is to be sure we have a reasonable starting point (of course initial lambdaM needs to be reasonable...)
    alphas<-c(alphaM,alphaG)
    for(i in 1:maxitT){    ### here the iteration begins
        nit<-nit+1;truemean<-exp(Bmu%*%alphaM)
        if (type=="quasi"){
            Dev<-DevPois(y,truemean)}
        if(type=="pseudo"){
            Dev<-((y-truemean)^2)/(truemean)}
        if(type != "pseudo" & type != "quasi" ){
            stop("the type of likelihood to use is not recognized")}
        zvec<-c(Bmu%*%alphaM+(y-exp(Bmu%*%alphaM))*exp(-Bmu%*%alphaM));W<-diag(as.vector(exp(-Bgamma%*%alphaG)*exp(Bmu%*%alphaM)))
        alphaM1<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%zvec
        if(link=="exp"){
            W<-diag(x=.5,length(y));zvec<-Bgamma%*%alphaG+((Dev-c(exp(Bgamma%*%alphaG)))/exp(Bgamma%*%alphaG))}
        if(link=="NB"){
            W<-diag(0.5*as.vector((exp(Bgamma%*%alphaG)/(1+exp(Bgamma%*%alphaG)))^2))
            zvec<-Bgamma%*%alphaG+((Dev-1-exp(Bgamma%*%alphaG))/exp(Bgamma%*%alphaG))}
        alphaG1<-ginv(t(Bgamma)%*%W%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%W%*%zvec
        alphaGeps<-sum(((alphaG1-alphaG)/alphaG1)^2)
        alphas1<-c(alphaM1,alphaG1)
        alphaeps<-sum(((alphas1-alphas)/alphas)^2)
        alphas<-alphas1;alphaM<-alphas[1:k];alphaG<-alphas[(k+1):(k+kt)] #alphaGeps<-10^{-j*3}
        if (alphaeps <= tol){
            conv <- TRUE
            break
        }
        if (i == maxitT){
            if(i != 1) warning("Iteration limit reached without convergence")}
    }
    # fill the empty result list
    res$alphaM<-alphaM
    res$alphaG<-alphaG
    res$meanest <- as.vector(exp(Bmu%*%alphaM))
    res$dispest <- as.vector(exp(Bgamma%*%alphaG)^{-1})
    res$numberOfIteration<-nit
    res$conv<-conv
    res
}


pb1sGCV3<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",link="exp"){
poissonBoth1step(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link)$GCV3
}
pb1sGCV1<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",link="exp"){
poissonBoth1step(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link)$GCV1
}
pb1sAIC<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",link="exp"){
poissonBoth1step(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link)$AIC
}
pb1sAIC3<-function(lambdas,y,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",link="exp"){
poissonBoth1step(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdas[1],lambdaG=lambdas[2],maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link)$AIC3
}



Poisson1sGCV<-function(y,Bmu,Bgamma,alphaM=0,alphaG=0,lambdaM,lambdaG,maxitM=85,maxitG=95,maxitT=125,maxitOUT=7,penordM=2,penordG=2,tol=10^-10,type="quasi",crit="GCV1",minlM=0.005,minlG=0.005,maxlM=10000,maxlG=10000,link="exp"){
##### this chooses 'optimal' value for the smoothing parameters and fits the model 
##### the optimal values for lambdaM and lambdaG are chosen simultaneously via numerical minimization
if(crit!="GCV1" & crit!="GCV3" & crit!="AIC" & crit!="AIC3"){
warning("unidentified criterion, use GCV3 instead")
crit<-"GCV3"
}
if(link == "nb") link<-"NB"
if(link != "exp" & link != "NB"){ warning("which link? I use exp");link<-"exp"}
n<-length(y);res<-list();k<-ncol(Bmu);kt<-ncol(Bgamma)
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bgamma))
}
DistMu<-diag(k);DistGam<-diag(kt);nit<-0
if(penordM>0) DistMu<-diff(DistMu,diff=penordM)
if(penordG>0) DistGam<-diff(DistGam,diff=penordG)
Pmu<-t(DistMu)%*%DistMu;Pgam<-t(DistGam)%*%DistGam
lambdaMinit<-lambdaM;lambdaGinit<-lambdaG
meanest<-exp(Bmu%*%alphaM);gammaest<-exp(Bgamma%*%alphaG)
for(i in 1:maxitOUT){
# lambdaM<-lambdaM+10
nit<-nit+1;alphaMinit<-alphaM;alphaGinit<-alphaG
if(crit=="AIC"){
ee<-optim(c(lambdaM,lambdaG),pb1sAIC,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),pb1sAIC,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
if(crit=="AIC3"){
ee<-optim(c(lambdaM,lambdaG),pb1sAIC3,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),pb1sAIC3,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
if(crit=="GCV1"){ee<-optim(c(lambdaM,lambdaG),pb1sGCV1,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),pb1sGCV1,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
if(crit=="GCV3"){ee<-optim(c(lambdaM,lambdaG),pb1sGCV3,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),pb1sGCV3,y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
lambdaM<-ee$par[1];lambdaG<-ee$par[2]
qq<-poissonBoth(y=y,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,lambdaM=lambdaM,lambdaG=lambdaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,link=link)
alphaM<-qq$alphaM;alphaG<-qq$alphaG
if (sum((alphaGinit-alphaG)^2)/sum(alphaG^2)<=tol & sum((alphaMinit-alphaM)^2)/sum(alphaM^2)<=tol  ) break
if (i == maxitOUT) warning("Iteration limit reached without convergence")
}
res$lM<-res$lambdaM<-lambdaM;res$lG<-res$lambdaG<-lambdaG;res$vals<-c(lambdaM,lambdaG)
truemean<-exp(Bmu%*%alphaM);Dev<-DevPois(y,truemean)
res$dfM<-qq$dfM;res$dfG<-qq$dfG
# W<-diag(as.vector(exp(-Bgamma%*%alphaG)*truemean))
# res$hmean<-Bmu%*%ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
# res$dfM<-dfM<-sum(diag(ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%(Bmu)))
# res$hdisp<-0.5*Bgamma%*%ginv(0.5*t(Bgamma)%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)
# res$dfG<-dfG<-sum(diag(ginv(0.5*t(Bgamma)%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%Bgamma))
res$alphaM<-as.vector(alphaM);res$alphaG<-as.vector(alphaG);res$numberOfIteration<-nit
res$numberOfIteration<-nit
res$gammaest<-qq$gammaest;res$meanest<-qq$meanest
res$AIC<-qq$AIC;res$GCV<-qq$GCV3
res
}








############################
#### two steps procedures

PoissonMean<-function(y,Bmu,thetatrue,alphaM=0,lambdaM,maxit=80,penordM=2,tol=10^-10){
#### fit, for a given lambdaM and a given dispersion function, the model for the mean of a normal
n<-length(y);res<-list();k<-ncol(Bmu)
DistMu<-diag(k);nit<-0
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
if(penordM>0) DistMu<-diff(DistMu,diff=penordM)
Pmu<-t(DistMu)%*%DistMu
for(l in 1:maxit) {
nit<-nit+1
W<-diag(as.vector((1/thetatrue)*exp(Bmu%*%alphaM)))
alphaM1<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%(Bmu%*%alphaM+(y-exp(Bmu%*%alphaM))/c(exp(Bmu%*%alphaM)))
alphaMeps<-sum((alphaM-alphaM1)^2)/sum(alphaM^2);alphaM<-alphaM1
if (alphaMeps<= tol) break
if (l == maxit){
if(l != 1) warning("Iteration limit reached without convergence")
}
}
W<-diag(as.vector((exp(Bmu%*%alphaM))/thetatrue))
truemean<-exp(Bmu%*%alphaM)
# res$hmean<-Bmu%*%ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
res$alphaM<-alphaM
res$numberOfIteration<-nit
res$df<-sum(diag(ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%Bmu))
res$deviance<-Dev<-DevPois(y,truemean)
res$AIC<-2*sum((truemean-y*log(truemean))/thetatrue)+2*res$df
res$GCV<-n*sum(Dev/thetatrue)/((n-res$df)^2)
res
}




PoissonMeanAIC<-function(lambda,y,Bmu,thetatrue,alphaM=0,maxit=80,penordM=2,tol=10^-10){
### given a lambda computes the AIC, function needed for numerical optimization
pp<-PoissonMean(y=y,Bmu=Bmu,thetatrue=thetatrue,alphaM=alphaM,lambdaM=lambda,maxit=maxit,penordM=penordM,tol=tol)$AIC
pp
}
PoissonMeanGCV<-function(lambda,y,Bmu,thetatrue,alphaM=0,maxit=80,penordM=2,tol=10^-10){
### given a lambda computes the GCV, function needed for numerical optimization
pp<-PoissonMean(y=y,Bmu=Bmu,thetatrue=thetatrue,alphaM=alphaM,lambdaM=lambda,maxit=maxit,penordM=penordM,tol=tol)$GCV
pp
}

PoissonMeanGCVgrid<-function(lambdagrid,y,Bmu,thetatrue,alphaM=0,maxit=80,penordM=2,tol=10^-10){
### computes the GCV for a certain grid of lambda values, function needed for grid search
pp<-list();GCV<-NULL
for(i in 1:length(lambdagrid)){
GCV[i]<-PoissonMean(y=y,Bmu=Bmu,thetatrue=thetatrue,alphaM=alphaM,lambdaM=lambdagrid[i],maxit=maxit,penordM=penordM,tol=tol)$GCV}
pp$GCV<-GCV;pp$lambdaseq<-lambdagrid;pp
}

PoissonMeanAICgrid<-function(lambdagrid,y,Bmu,thetatrue,alphaM=0,maxit=80,penordM=2,tol=10^-10){
### computes the AIC for a certain grid of lambda values, function needed for grid search
pp<-list();AIC<-NULL
for(i in 1:length(lambdagrid)){
AIC[i]<-PoissonMean(y=y,Bmu=Bmu,thetatrue=thetatrue,alphaM=alphaM,lambdaM=lambdagrid[i],maxit=maxit,penordM=penordM,tol=tol)$AIC}
pp$AIC<-AIC;pp$lambdaseq<-lambdagrid;pp
}




PoissonDisp<-function(y,Bgamma,truemean,alphaG=0,lambdaG,maxit=80,penordG=2,tol=10^-10,type="quasi",link="exp"){
### fit the dispersion model for a given lambdaG value
if(link == "nb") link<-"NB"
if(link != "exp" & link != "NB"){ warning("which link? I use exp");link<-"exp"}
n<-length(y);res<-list();kt<-ncol(Bgamma)
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bgamma))}
DistGam<-diag(kt);nit<-0
if(penordG>0) DistGam<-diff(DistGam,diff=penordG)
Pgam<-t(DistGam)%*%DistGam
if (type=="quasi"){
Dev<-DevPois(y,truemean)}
if(type=="pseudo"){
Dev<-((y-truemean)^2)/(truemean)}
if(type != "pseudo" & type != "quasi" ){
stop("the type of likelihood to use is not recognized")
}
for(j in 1:maxit){
nit<-nit+1
if(link=="exp"){
W<-diag(x=.5,length(y));q<-Bgamma%*%alphaG+((Dev-c(exp(Bgamma%*%alphaG)))/exp(Bgamma%*%alphaG))}
if(link=="NB"){
W<-diag(as.vector(0.5*(exp(Bgamma%*%alphaG)/(1+exp(Bgamma%*%alphaG)))^2))
q<-Bgamma%*%alphaG+((Dev-1-exp(Bgamma%*%alphaG))/exp(Bgamma%*%alphaG))}
alphaG1<-ginv(t(Bgamma)%*%W%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%W%*%q
alphaGeps<-sum(((alphaG1-alphaG)/alphaG1)^2)
alphaG<-alphaG1 #alphaGeps<-10^{-j*3}
if (alphaGeps <= tol) break
if (j == maxit){
if(j != 1) warning("Iteration limit reached without convergence")
}}
res$hdisp<-Bgamma%*%ginv(t(Bgamma)%*%W%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%W
if(link=="exp"){res$gammaest<-gammaest<-exp(Bgamma%*%alphaG);W<-diag(x=.5,length(y))}
if(link=="NB") {res$gammaest<-gammaest<-(1+exp(Bgamma%*%alphaG))
W<-diag(as.vector(0.5*(exp(Bgamma%*%alphaG)/(1+exp(Bgamma%*%alphaG)))^2))}
res$alphaG<-as.vector(alphaG)
res$numberOfIteration<-nit
res$deviance<-Dev;res$Method<-type
df<-sum(diag(ginv(t(Bgamma)%*%W%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%W%*%Bgamma))
res$df<-df
res$AIC<-sum(log(gammaest)+Dev/gammaest)+(2*df)
res$GCV<-n*sum(log(gammaest)+Dev/gammaest-log(Dev)-1)/((n-df)^2)
res
}







PoissonDispAIC<-function(lambda,y,Bgamma,truemean,alphaG=0,maxit=80,penordG=2,tol=10^-8,type="quasi",link="exp"){
### given a lambda computes the AIC, function needed for numerical optimization
PoissonDisp(y=y,Bgamma=Bgamma,truemean=truemean,alphaG=alphaG,lambdaG=lambda,maxit=maxit,penordG=penordG,tol=tol,type=type,link=link)$AIC
}
PoissonDispGCV<-function(lambda,y,Bgamma,truemean,alphaG=0,maxit=80,penordG=2,tol=10^-8,type="quasi",link="exp"){
### given a lambda computes the GCV, function needed for numerical optimization
PoissonDisp(y=y,Bgamma=Bgamma,truemean=truemean,alphaG=alphaG,lambdaG=lambda,maxit=maxit,penordG=penordG,tol=tol,type=type,link=link)$GCV
}

PoissonDispGCVgrid<-function(lambdagrid,y,Bgamma,truemean,alphaG=0,maxit=80,penordG=2,tol=10^-8,type="quasi",link="exp"){
### computes the GCV for a certain grid of lambda values, function needed for grid search
pp<-list();GCV<-NULL
for(i in 1:length(lambdagrid)){
GCV[i]<-PoissonDisp(y=y,Bgamma=Bgamma,truemean=truemean,alphaG=alphaG,lambdaG=lambdagrid[i],maxit=maxit,penordG=penordG,tol=tol,type=type,link=link)$GCV}
pp$GCV<-GCV;pp$lambdaseq<-lambdagrid;pp
}

PoissonDispAICgrid<-function(lambdagrid,y,Bgamma,truemean,alphaG=0,maxit=80,penordG=2,tol=10^-8,type="quasi",link="exp"){
### computes the AIC for a certain grid of lambda values, function needed for grid search
pp<-list();AIC<-NULL
for(i in 1:length(lambdagrid)){
AIC[i]<-PoissonDisp(y=y,Bgamma=Bgamma,truemean=truemean,alphaG=alphaG,lambdaG=lambdagrid[i],maxit=maxit,penordG=penordG,tol=tol,type=type,link=link)$AIC}
pp$AIC<-AIC;pp$lambdaseq<-lambdagrid;pp
}


poisson2sGCVnum<-function(y,Bmu,Bgamma,lambdaM,lambdaG,maxitM=75,maxitG=85,maxitT=8,penordM=2,penordG=2,tol=10^-6,type="quasi",alphaM=0,alphaG=0,minlM=0.005,minlG=0.005,maxlM=8000,maxlG=8000,critM="GCV",critG="GCV",link="exp"){
##### this chooses 'optimal' value for the smoothing parametersin two steps and fits the model 
##### the smoothing parameter choice is done via numerical minimization of GCV or AIC
if(critM !="AIC" & critM !="GCV") {warning("what is the criterion to choose lambdaM, I use GCV");critM<-"GCV"}
if(critG !="AIC" & critG !="GCV") {warning("what is the criterion to choose lambdaG, I use GCV");critM<-"GCV"}
if(link == "nb") link<-"NB"
if(link != "exp" & link != "NB"){ warning("which link? I use exp");link<-"exp"}
n<-length(y);res<-list();res$conv<-"false";k<-ncol(Bmu);kt<-ncol(Bgamma)
lM0<-lambdaM;lG0<-lambdaG;nit<-0
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bgamma))
}
if(link=="exp"){gammaest<-exp(Bgamma%*%alphaG)}
if(link=="NB") {gammaest<-(1+exp(Bgamma%*%alphaG))}
meanest<-exp(Bmu%*%alphaM)
for(i in 1:maxitT){
nit<-nit+1;alphaMinit<-alphaM;alphaGinit<-alphaG
if(critM=="GCV"){
sl<-optim(p=lambdaM,PoissonMeanGCV,y=y,Bmu=Bmu,thetatrue=gammaest,alphaM=alphaM,maxit=maxitM,penordM=penordM,tol=tol,method="L-BFGS-B",upper=maxlM,lower=minlM)
if(sl$conv > 0.5){
warning(paste("at mean iteration",i,"optim did not converge, I restart from initial lambda"))
sl<-optim(p=lM0,PoissonMeanGCV,y=y,Bmu=Bmu,thetatrue=gammaest,alphaM=alphaM,maxit=maxitM,penordM=penordM,tol=tol,method="L-BFGS-B",upper=maxlM,lower=minlM)}}
if(critM=="AIC"){
sl<-optim(p=lambdaM,PoissonMeanAIC,y=y,Bmu=Bmu,thetatrue=gammaest,alphaM=alphaM,maxit=maxitM,penordM=penordM,tol=tol,method="L-BFGS-B",upper=maxlM,lower=minlM)
if(sl$conv > 0.5){
warning(paste("at mean iteration",i,"optim did not converge, I restart from initial lambda"))
sl<-optim(p=lM0,PoissonMeanAIC,y=y,Bmu=Bmu,thetatrue=gammaest,alphaM=alphaM,maxit=maxitM,penordM=penordM,tol=tol,method="L-BFGS-B",upper=maxlM,lower=minlM)}}
lambdaM<-sl$par
# print(c(sl$conv,sl$par,lambdaM))
pm<-PoissonMean(y=y,Bmu=Bmu,thetatrue=gammaest,lambdaM=lambdaM,alphaM=alphaM,maxit=maxitM,penordM=penordM,tol=tol)
meanest<-exp(Bmu%*%pm$alphaM);alphaM<-pm$alphaM
if(critG=="GCV"){
sl<-optim(p=lambdaG,PoissonDispGCV,y=y,Bgamma=Bgamma,truemean=meanest,alphaG=alphaG,maxit=maxitG,penordG=penordG,tol=tol,type=type,method="L-BFGS-B",upper=maxlG,lower=minlG,link=link)
if(sl$conv > 0.5){
warning(paste("at disp iteration",i,"optim did not converge for dispersion, I restart from initial lambda"))
sl<-optim(p=lG0,PoissonDispGCV,y=y,Bgamma=Bgamma,truemean=meanest,alphaG=alphaG,maxit=maxitG,penordG=penordG,tol=tol,type=type,method="L-BFGS-B",upper=maxlG,lower=minlG,link=link)}}
if(critG=="AIC"){
sl<-optim(p=lambdaG,PoissonDispAIC,y=y,Bgamma=Bgamma,truemean=meanest,alphaG=alphaG,maxit=maxitG,penordG=penordG,tol=tol,type=type,method="L-BFGS-B",upper=maxlG,lower=minlG,link=link)
if(sl$conv > 0.5){
warning(paste("at disp iteration",i,"optim did not converge for dispersion, I restart from initial lambda"))
sl<-optim(p=lG0,PoissonDispAIC,y=y,Bgamma=Bgamma,truemean=meanest,alphaG=alphaG,maxit=maxitG,penordG=penordG,tol=tol,type=type,method="L-BFGS-B",upper=maxlG,lower=minlG,link=link)}}
lambdaG<-sl$par
pg<-PoissonDisp(y=y,Bgamma=Bgamma,truemean=meanest,lambdaG=lambdaG,alphaG=alphaG,maxit=maxitG,penordG=penordG,tol=tol,type=type,link=link)
alphaG<-pg$alphaG;gammaest<-pg$gammaest
# print(c(sl$conv,sl$par,lambdaG))
if (sum((alphaGinit-alphaG)^2)/sum(alphaG^2)<=tol & sum((alphaMinit-alphaM)^2)/sum(alphaM^2)<=tol){
res$conv<-"true"
break}
if(nit == maxitT) print('number of maxitT iteration reached with no convergence')
}
res$numberOfIteration<-nit
res$alphaM<-as.vector(alphaM);res$alphaG<-as.vector(alphaG)
res$lM<-res$lambdaM<-lambdaM;res$lG<-res$lambdaG<-lambdaG;res$vals<-c(lambdaM,lambdaG)
# res$hmean<-pm$hmean;res$hdisp<-pg$hdisp
res$meanest<-meanest;res$gammaest<-gammaest
res$deviance<-pg$deviance
res$dfM<-pm$df;res$dfG<-pg$df
res
}










poisson2sGCV<-function(y,Bmu,Bgamma,seqM,seqG,lenM,lenG,maxitM=75,maxitG=75,maxitT=125,penordM=2,penordG=2,tol=10^-6,tolSTOP=10^-4,type="quasi",alphaM=0,alphaG=0,minlM=0.005,minlG=0.005,maxlM=8000,maxlG=8000,critM="GCV",critG="GCV",link="exp",plotM=FALSE,plotG=FALSE,uppars=TRUE){
##### this chooses 'optimal' value for the smoothing parametersin two steps and fits the model 
##### the smoothing parameter choice is done via a grid search (the same as in the normal model)
if(critM !="AIC" & critM !="GCV") {warning("what is the criterion to choose lambdaM, I use GCV");critM<-"GCV"}
if(critG !="AIC" & critG !="GCV") {warning("what is the criterion to choose lambdaG, I use GCV");critM<-"GCV"}
if(link == "nb") link<-"NB"
if(link != "exp" & link != "NB"){ warning("which link? I use exp");link<-"exp"}
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bgamma))}
n<-length(y);res<-list();res$conv<-"false";k<-ncol(Bmu);kt<-ncol(Bgamma)
llM<-seq(min(seqM),max(seqM),length=lenM)
llG<-seq(min(seqG),max(seqG),length=lenG)
# diffsMean<-NULL;diffsDisp<-NULL
n<-length(y);dispest<-as.vector(rep(1,n))
meanest<-as.vector(rep(mean(y),n))
lambdaM<-10^-5;lambdaG<-10^-5
extremesM<-NULL;extremesG<-NULL
ntot<-0;numsM<-NULL;numsG<-NULL
meanest<-rep(mean(y),l=n)
for(j in 2:6){
ntot<-ntot+1
checkM<-NULL;checkM[1]<-10^-5
checkG<-NULL;checkG[1]<-10^-5
for(i in 1:5){
alphaMinit<-as.vector(alphaM);alphaGinit<-as.vector(alphaG)
if(critM=="GCV"){
gcv<-PoissonMeanGCVgrid(lambdagrid=llM,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit)}
if(critM=="AIC"){
gcv<-PoissonMeanAICgrid(lambdagrid=llM,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit);gcv$GCV<-gcv$AIC}
gcvM<-min(gcv$GCV)
if (plotM==TRUE){
plot(llM,gcv$GCV,type="l")
}
ss<-min(seq(1,lenM)[gcv$GCV==min(gcv$GCV)])
lM<-llM[ss]
checkM[i+1]<-lM#min(gcv$GCV)
if(checkM[i]/checkM[i+1]<1.05 & checkM[i]/checkM[i+1]>0.95){
llM<-sort(seq(lM*0.5,lM*2,length=lenM))
extremesM<-c(extremesM,as.character(ss))
break}
else {
paceM<-llM[2]-llM[1];lM<-llM[ss]
if(ss==lenM){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL;diffGCVFin<-NULL
eps<-min(gcv$GCV)*10^-3;sum5<-NULL
for( k in 1:(lenM-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenM])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))}
sum5<-diffsTOT
ns<-min(ss,seq(1,length(llM))[sum5>5.5])
#ns<-min(ss,seq(1,length(llM))[diffGCV<eps])
rm(sum5)
if(ns != ss){
lM<-llM[ns];extremesM<-c(extremesM,"max - cutted")
llM<-sort(seq(from=max(minlM,min(llM[1]*0.8,lM-6*paceM)),by=1.6*paceM,l=lenM))
break}
if(ns==ss){
lM<-llM[ss]
llM<- sort(seq(from=max(minlM,lM-10*paceM),by=1.2*paceM,l=lenM))
#llM<-sort(seq(from=min(llM[1]*.9,lM-7*paceM),by=1.2*paceM,l=lenM))
extremesM<-c(extremesM,"max")
}}
if(ss==1){
if(lM==minlM){
llM<-sort(seq(minlM,by=1.2*paceM,l=lenM))}
else{
llM<-sort(seq(max(minlM,lM-(0.8*lenM*paceM)),by=1.2*paceM,l=lenM))}
extremesM<-c(extremesM,"min")
}
if(ss!=1 & ss!=lenM){
llM<-sort(seq(from=max(minlM,lM-6*paceM),to=max(minlM+10*paceM,lM+6*paceM),l=lenM))
extremesM<-c(extremesM,as.character(ss))
}}}
extremesM<-c(extremesM," ")
if(lM>maxlM){
lM<-maxlM
warning("lM higher then maxlM forced to be maxlM")
}
numsM<-c(numsM,i)
mm<-PoissonMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lM,maxit=maxitM,penordM=penordM,tol=tol)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2)
alphaM<-mm$alphaM;meanest<-as.vector(exp(Bmu%*%mm$alphaM))
rm(gcv,ss)
for(i in 1:5){
if(critG=="GCV"){
gcv<-PoissonDispGCVgrid(lambdagrid=llG,y=y,Bgamma=Bgamma,truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,link=link)}
if(critG=="AIC"){
gcv<-PoissonDispAICgrid(lambdagrid=llG,y=y,Bgamma=Bgamma,truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,link=link);gcv$GCV<-gcv$AIC}
gcvG<-min(gcv$GCV)
if(plotG==TRUE){
plot(llG,gcv$GCV,type="l")}
ss<-min(seq(1,lenG)[gcv$GCV==min(gcv$GCV)]);lG<-llG[ss]
checkG[i+1]<-lG#min(gcv$GCV)
if(checkG[i]/checkG[i+1]<1.05 & checkG[i]/checkG[i+1]>.95){
llG<-sort(seq(0.5*lG,lG*2,length=lenG))
extremesG<-c(extremesG,as.character(ss))
break}
else{
paceG<-llG[2]-llG[1]
if(ss==lenG){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL;diffGCVFin<-NULL
eps<-min(gcv$GCV)*10^-2
for( k in 1:(lenG-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenG])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))
}
sum6<-diffsTOT
#sum(diffGCV1<eps,diffGCV2<eps,diffGCV3<eps,diffGCV4<eps,diffGCV5<eps)
ns<-min(ss,seq(1,length(llG))[sum6>5.5])
rm(sum6)
if(ns != ss){
lG<-llG[ns];extremesG<-c(extremesG,"max - cutted")
llG<-sort(seq(from=max(minlG,min(llG[1]*0.8,lG-6*paceG)),by=1.6*paceG,l=lenG))
break}
if(ns==ss){
lG<-llG[ss];llG<- sort(seq(from=max(minlG,lG-10*paceG),by=1.2*paceG,l=lenG))
extremesG<-c(extremesG,"max")
}}
if(ss==1){
if(lG==0){
llG<-sort(seq(max(minlG,lG-(0.8*lenG*paceG)),by=1.2*paceG,l=lenG))}
else{
llG<-sort(seq(from=max(minlG,lG-(0.8*lenG*paceG)),by=1.2*paceG,l=lenG))}
extremesG<-c(extremesG,"min")
}
if(ss!=1 & ss!=lenG){
llG<-sort(seq(from=max(minlG,lG-6*paceG),to=max(minlG+10*paceG,lG+6*paceG),l=lenG))
extremesG<-c(extremesG,as.character(ss))
}}
}
extremesG<-c(extremesG," ")
if(lG>maxlG){
lG<-maxlG
warning("lG higher then maxlG forced to be maxlG")}
numsG<-c(numsG,i)
zz<-PoissonDisp(y=y,Bgamma=Bgamma,truemean=meanest,penordG=penordG,alphaG=alphaGinit,lambdaG=lG,type=type,link=link,maxit=maxitG)
alphaGdiff<-sum(((zz$alphaG-alphaG)/(alphaG))^2)
dispest<-zz$gammaest;alphaG<-zz$alphaG
rm(gcv,ss)
lambdaM<-c(lambdaM,lM);lambdaG<-c(lambdaG,lG)
if (uppars == TRUE){
alphaGinit<-alphaG;alphaMinit<-alphaM}
if(max(alphaGdiff,alphaMdiff)<tolSTOP) break
}
res$vals<-c(lM,lG)
res$lM<-res$lambdaM<-lM
res$lG<-res$lambdaG<-lG
res$ntri<-j-1
res$sequencesM<-extremesM
res$sequencesG<-extremesG
res$meanest<-as.vector(meanest)
res$gammaest<-as.vector(dispest)
res$lambdasM<-lambdaM[1:length(lambdaM)]
res$lambdasG<-lambdaG[1:length(lambdaG)]
res$alphaM<-as.vector(mm$alphaM)
res$alphaG<-as.vector(zz$alphaG)
res$dfM<-mm$df
res$dfG<-zz$df
res$vecNumM<-numsM;res$vecNumG<-numsG
res$Dev<-zz$deviance
res
}



poissonMeanOnlyGCV<-function(y,Bmu,seqM,lenM,maxitM=75,penordM=2,tol=10^-6,tolSTOP=10^-4,alphaM=0,minlM=0.005,maxlM=8000,critM="GCV",plotM=FALSE,uppars=TRUE){
##### this function only estimates the mean function keeping the phi=1
##### the smoothing parameter choice is done via a grid search (as above)
if(critM !="AIC" & critM !="GCV") {warning("what is the criterion to choose lambdaM, I use GCV");critM<-"GCV"}
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
n<-length(y);res<-list();res$conv<-"false";k<-ncol(Bmu)
llM<-seq(min(seqM),max(seqM),length=lenM)
n<-length(y);dispest<-as.vector(rep(1,n))
meanest<-as.vector(rep(mean(y),n))
lambdaM<-10^-5;extremesM<-NULL;ntot<-0;numsM<-NULL
for(j in 2:5){
ntot<-ntot+1
checkM<-NULL;checkM[1]<-10^-5
for(i in 1:5){
alphaMinit<-as.vector(alphaM)
if(critM=="GCV"){
gcv<-PoissonMeanGCVgrid(lambdagrid=llM,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit)}
if(critM=="AIC"){
gcv<-PoissonMeanAICgrid(lambdagrid=llM,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit);gcv$GCV<-gcv$AIC}
gcvM<-min(gcv$GCV)
if (plotM==TRUE){
plot(llM,gcv$GCV,type="l")
}
ss<-min(seq(1,lenM)[gcv$GCV==min(gcv$GCV)])
lM<-llM[ss]
checkM[i+1]<-lM#min(gcv$GCV)
if(checkM[i]/checkM[i+1]<1.05 & checkM[i]/checkM[i+1]>0.95){
llM<-sort(seq(lM*0.5,lM*2,length=lenM))
extremesM<-c(extremesM,as.character(ss))
break}
else {
paceM<-llM[2]-llM[1];lM<-llM[ss]
if(ss==lenM){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL;diffGCVFin<-NULL
eps<-min(gcv$GCV)*10^-3;sum5<-NULL
for( k in 1:(lenM-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenM])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))}
sum5<-diffsTOT
ns<-min(ss,seq(1,length(llM))[sum5>5.5])
rm(sum5)
if(ns != ss){
lM<-llM[ns];extremesM<-c(extremesM,"max - cutted")
llM<-sort(seq(from=max(minlM,min(llM[1]*0.8,lM-6*paceM)),by=1.6*paceM,l=lenM))
break}
if(ns==ss){
lM<-llM[ss]
llM<- sort(seq(from=max(minlM,lM-10*paceM),by=1.2*paceM,l=lenM))
extremesM<-c(extremesM,"max")
}}
if(ss==1){
if(lM==minlM){
llM<-sort(seq(minlM,by=1.2*paceM,l=lenM))}
else{
llM<-sort(seq(max(minlM,lM-(0.8*lenM*paceM)),by=1.2*paceM,l=lenM))}
extremesM<-c(extremesM,"min")
}
if(ss!=1 & ss!=lenM){
llM<-sort(seq(from=max(minlM,lM-6*paceM),to=max(minlM+10*paceM,lM+6*paceM),l=lenM))
extremesM<-c(extremesM,as.character(ss))
}}}
extremesM<-c(extremesM," ")
if(lM>maxlM){
lM<-maxlM
warning("lM higher then maxlM forced to be maxlM")
}
numsM<-c(numsM,i)
mm<-PoissonMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lM,maxit=maxitM,penordM=penordM,tol=tol)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2)
alphaM<-mm$alphaM;meanest<-as.vector(exp(Bmu%*%mm$alphaM))
rm(gcv,ss)
lambdaM<-c(lambdaM,lM)
if (uppars == TRUE){
alphaMinit<-alphaM}
if(max(alphaMdiff)<tolSTOP) break
}
res$lM<-res$lambdaM<-lM;res$ntri<-j-1
res$sequencesM<-extremesM;res$meanest<-as.vector(meanest)
res$lambdasM<-lambdaM[1:length(lambdaM)]
res$alphaM<-as.vector(mm$alphaM);res$dfM<-mm$df
res$vecNumM<-numsM
res$Dev<-DevPois(y,meanest)
res
}



#### in this one, we also estimate a constant dispersion parameter...


poissonMeanOnlyGCV2<-function(y,Bmu,seqM,lenM,maxitM=75,penordM=2,tol=10^-6,tolSTOP=10^-4,alphaM=0,minlM=0.005,maxlM=8000,critM="GCV",plotM=FALSE,uppars=TRUE){
  ##### the function estimates the mean function as a smooth function 
  ##### it estimates the dispersion parameter phi as a constant
  ##### the smoothing parameter choice is done via a grid search (as above)
if(critM !="AIC" & critM !="GCV") {warning("what is the criterion to choose lambdaM, I use GCV");critM<-"GCV"}
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
n<-length(y);res<-list();res$conv<-"false";k<-ncol(Bmu)
llM<-seq(min(seqM),max(seqM),length=lenM)
n<-length(y);dispest<-as.vector(rep(1,n))
meanest<-as.vector(rep(mean(y),n))
lambdaM<-10^-5;extremesM<-NULL;ntot<-0;numsM<-NULL
for(j in 2:5){
ntot<-ntot+1
checkM<-NULL;checkM[1]<-10^-5
for(i in 1:5){
alphaMinit<-as.vector(alphaM)
if(critM=="GCV"){
gcv<-PoissonMeanGCVgrid(lambdagrid=llM,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit)}
if(critM=="AIC"){
gcv<-PoissonMeanAICgrid(lambdagrid=llM,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit);gcv$GCV<-gcv$AIC}
gcvM<-min(gcv$GCV)
if (plotM==TRUE){
plot(llM,gcv$GCV,type="l")
}
ss<-min(seq(1,lenM)[gcv$GCV==min(gcv$GCV)])
lM<-llM[ss]
checkM[i+1]<-lM#min(gcv$GCV)
if(checkM[i]/checkM[i+1]<1.05 & checkM[i]/checkM[i+1]>0.95){
llM<-sort(seq(lM*0.5,lM*2,length=lenM))
extremesM<-c(extremesM,as.character(ss))
break}
else {
paceM<-llM[2]-llM[1];lM<-llM[ss]
if(ss==lenM){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL;diffGCVFin<-NULL
eps<-min(gcv$GCV)*10^-3;sum5<-NULL
for( k in 1:(lenM-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenM])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))}
sum5<-diffsTOT
ns<-min(ss,seq(1,length(llM))[sum5>5.5])
rm(sum5)
if(ns != ss){
lM<-llM[ns];extremesM<-c(extremesM,"max - cutted")
llM<-sort(seq(from=max(minlM,min(llM[1]*0.8,lM-6*paceM)),by=1.6*paceM,l=lenM))
break}
if(ns==ss){
lM<-llM[ss]
llM<- sort(seq(from=max(minlM,lM-10*paceM),by=1.2*paceM,l=lenM))
extremesM<-c(extremesM,"max")
}}
if(ss==1){
if(lM==minlM){
llM<-sort(seq(minlM,by=1.2*paceM,l=lenM))}
else{
llM<-sort(seq(max(minlM,lM-(0.8*lenM*paceM)),by=1.2*paceM,l=lenM))}
extremesM<-c(extremesM,"min")
}
if(ss!=1 & ss!=lenM){
llM<-sort(seq(from=max(minlM,lM-6*paceM),to=max(minlM+10*paceM,lM+6*paceM),l=lenM))
extremesM<-c(extremesM,as.character(ss))
}}}
extremesM<-c(extremesM," ")
if(lM>maxlM){
lM<-maxlM
warning("lM higher then maxlM forced to be maxlM")
}
numsM<-c(numsM,i)
mm<-PoissonMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lM,maxit=maxitM,penordM=penordM,tol=tol)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2)
alphaM<-mm$alphaM;meanest<-as.vector(exp(Bmu%*%mm$alphaM))
rm(gcv,ss)
dispest<-rep(sum(mm$dev)/(length(y)-mm$df),l=length(y)) #### dispersion estimate
lambdaM<-c(lambdaM,lM)
if (uppars == TRUE){
alphaMinit<-alphaM}
if(max(alphaMdiff)<tolSTOP) break
}
res$lM<-res$lambdaM<-lM;res$ntri<-j-1
res$sequencesM<-extremesM;res$meanest<-as.vector(meanest)
res$lambdasM<-lambdaM[1:length(lambdaM)]
res$alphaM<-as.vector(mm$alphaM);res$dfM<-mm$df
res$vecNumM<-numsM
res$Dev<-DevPois(y,meanest)
res$dispest<-dispest
res
}





# # The Fabric data as in the paper
# # 
# # 
# # library(gamlss)
# # data(fabric)
# # fabric<-fabric[sort.list(fabric$x),]
# # 
# # xx<-fabric$x;yy<-fabric$y
# # dx<-(max(xx)-min(xx))/7;kn<-seq(min(xx)-3*dx,max(xx)+3*dx,by=dx)
# # basem<-splineDesign(knots=kn,x=xx,ord=(3+1),0*xx,outer=TRUE)
# # dx<-(max(xx)-min(xx))/5;kn<-seq(min(xx)-3*dx,max(xx)+3*dx,by=dx)
# # based<-splineDesign(knots=kn,x=xx,ord=(3+1),0*xx,outer=TRUE)
# # 
# # tt<-poisson2sGCV(y=yy,Bmu=basem,Bgamma=based,seqM=c(2.5,25),seqG=c(405,6825),lenM=51,lenG=91,maxlM=8000,maxlG=10000,minlM=10^-6,minlG=10^-8,critM="GCV",critG="GCV",link="exp",plotG=TRUE,plotM=TRUE)
# # 
# # n<-length(yy)
# # Dev<-2*as.vector(yy*log(yy)-yy-yy*log(tt$meanest)+tt$meanest)
# # ss<-seq(1,n)[yy==0];Dev[ss]<-2*tt$meanest[ss];rm(n)
# # 
# # par(mfrow=c(1,2))
# # plot(xx,yy,xlab="log(roll length)",ylab="Number of faults",bty="l",pch=19,cex=0.7,main="(a)")
# # points(xx,tt$meanest,col=4,type="l",lwd=2)
# # plot(xx,Dev,xlab="log(roll length)",ylab="Deviance",bty="l",pch=19,cex=0.7,main="(b)")
# # points(xx,tt$gammaest,col=4,type="l",lwd=2)







############################################################################
################                Binomial                     ################
############################################################################
############################################################################

### ranBin generates samples from the double binomial

ranBin<-function(size,n,prob,theta){
prob<-rep(prob,l=size);res<-list()
theta<-rep(theta,l=size)
n<-rep(n,l=size);p<-NULL
ss<-NULL
for( i in 1:size){
y<-seq(0,n[i]);p<-NULL
for(j in 1:length(y)){
p[j]<-sqrt(1/theta[i])*choose(n[i],y[j])*((((prob[i])^y[j])*((1-prob[i])^(n[i]-y[j])))^(1/theta[i]))  *((((y[j]/n[i])^y[j])*((1-y[j]/n[i])^(n[i]-y[j])))^(1-1/theta[i]))}
#p<-p/sum(p)
samp<-sample(y,size=1,prob=p)#,replace=TRUE
ss<-c(ss,samp)
}
res$sample<-ss
res$n<-n
#res$numberOFsample<-ncamp
#res$p<-p
res
}



#### deviance of a poisson for a given eta ....
DevBinE<-function(y,eta){
    n<-length(y)
    Dev<-2*(y*(log(y/(1-y)))-y*(eta)-log(1/(1-y))+log(1+exp(eta)))
    seq0<-seq(1,n)[y==0]
    seq1<-seq(1,n)[y==1]
    Dev[seq0]<-2*log(1+exp(eta[seq0]))
    Dev[seq1]<-2*(log(1+exp(eta[seq1]))-(eta[seq1]))
    Dev
}

#### .... and for a given probability
DevBinP<-function(y,p){ eta<-log(p/(1-p));n<-length(y);Dev<-2*(y*(log(y/(1-y)))-y*(eta)-log(1/(1-y))+log(1+exp(eta)));seq0<-seq(1,n)[y==0];seq1<-seq(1,n)[y==1];Dev[seq0]<-2*log(1+exp(eta[seq0]));Dev[seq1]<-2*(log(1+exp(eta[seq1]))-(eta[seq1]));Dev}




#### notes for most functions: 

#### you need to build the Bmu and Bgamma B-splines matrices before (via splines).....

#### for the binomial y should be the proportion of successes and n should be the number of trials ...

##### type can be either "quasi" or "pseudo"
  ## with quasi the procedure described in the paper is performed: use the deviance residuals as a responce for the dispersion model
  ## with pseudo a pseudo likelihood procedure is performed (Davidian and Carrol 1986): use the pearson residuals as a responce for the dispersion model  
  
##### link can be either "bb" or "exp" and it corresponds to the link function for the dispersion function
  ## with "exp" we have gamma=exp(xi), the default in the paper
  ## with "bb" we have gamma=h(xi) as in a beta binomial (see later)





eexpp<-function(x) exp(x)/(1+exp(x))  #### logitinv the inverse of the canonical link (logit)

#### the h() function  when using a link function like in the betabinomial case gamma=h(xi)=1+(n-1)*(xi/(1+xi))
#### the derivative h' and the ratio h'/h....

hf<-function(x,nn){
fun<-(1+nn*x)/(1+x);fun}

hfp<-function(x,nn){
fun<-(nn-1)/((1+x)^2);fun}

hfratio<-function(x,nn){
fun<-(nn-1)/((1+x)*(1+nn*x));fun}

#### the h() function  when using a link function like in the betabinomial case 
#### making sure the function is positive
#### gamma=h(xi)=1+(n-1)*(exp(xi)/(1+exp(xi)))
#### the derivative h' and the ratio h'/h....

hfex<-function(x,nn){
fun<- (1+nn*exp(x))/(1+exp(x));fun}

hfpex<-function(x,nn){
fun<-(nn-1)*exp(x)/((1+exp(x))^2);fun}

hfratioex<-function(x,nn){
fun<- ((nn-1)*exp(x))/((1+nn*exp(x))*(1+exp(x)));fun}





############################
#### one steps procedures



BinomialBoth<-function(y,num,Bmu,Bgamma,lambdaM,lambdaG,penordM=2,penordG=2,tol=10^-10,type="quasi",alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp",maxitT=250,maxitM=55,maxitG=45){
##### for given lambdaM and lambdaG estimate both functions
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bmu))
}
n<-length(y);meanest<-as.vector(Bmu%*%alphaM)
if(linktype=="exp") {dispest<-exp(Bgamma%*%alphaG)}
if(linktype=="BB")  {dispest<-hfex(Bgamma%*%alphaG,nn=num)}
# dispest<-hfratio(Bgamma%*%alphaG,nn=num)}
res<-NULL;ntrial<-NULL;ntot<-0
alphaMinit<-alphaM;alphaGinit<-alphaG
for(j in 1:maxitT){
mm<-BinomialMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lambdaM,tol=tol,maxit=maxitM)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2);alphaM<-mm$alphaM
meanest<-as.vector(exp(Bmu%*%mm$alphaM)/(1+exp(Bmu%*%mm$alphaM)))
zz<-BinomialDisp(y=y,num=num,Bgamma=Bgamma,lambdaG=lambdaG,truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,linktype=linktype,maxit=maxitG)
if(linktype=="exp") {dispest<-exp(Bgamma%*%zz$alphaG)}
if(linktype=="BB")  {dispest<-hfex(Bgamma%*%zz$alphaG,nn=num)}
# dispest<-hfratio(Bgamma%*%zz$alphaG,nn=num)}
alphaGdiff<-sum(((zz$alphaG-alphaG)/(alphaG))^2);alphaG<-zz$alphaG
alphaGinit<-alphaG;alphaMinit<-alphaM
if(max(alphaGdiff,alphaMdiff)<tolSTOP) break
}
res$vals<-c(lambdaM,lambdaG);res$ntri<-j-1
res$meanest<-as.vector(meanest);gammaest<-res$dispest<-as.vector(dispest)
res$alphaM<-as.vector(mm$alphaM);res$alphaG<-as.vector(zz$alphaG)
dfM<-res$dfM<-mm$df;dfG<-res$dfG<-zz$df
res$Dev<-dev<-zz$deviance;res$linktype<-linktype
res$AIC<-sum(log(gammaest)+dev/gammaest)+2*(dfM+dfG)
res$AIC1<-sum(log(gammaest/dev)+dev/gammaest-1)+sum(dev/gammaest)+2*(dfM+dfG)
res$GCV<-(n*sum(log(gammaest/dev)+dev/gammaest))/((n-(dfM+dfG))^2)
res$GCV1<-(n*n*(sum(log(gammaest/dev)+dev/gammaest-1)+sum(dev/gammaest)))/((2*n-(dfM+dfG))^2)
res
}

fisher_def_binomial<-function(y,num,Bmu,Bgamma,lambdaM,lambdaG,penordM=2,penordG=2,tol=10^-10,type="quasi",alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp",maxitT=250,maxitM=55,maxitG=45){
    ##### for given lambdaM and lambdaG estimate both functions
    if (length(alphaM)==1){
        if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
        else alphaM<-rep(alphaM,l=ncol(Bmu))
    }
    if (length(alphaG)==1){
        if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
        else alphaG<-rep(alphaG,l=ncol(Bmu))
    }
    n<-length(y)
    meanest<-as.vector(Bmu%*%alphaM)
    if(linktype=="exp") {dispest<-exp(Bgamma%*%alphaG)}
    if(linktype=="BB")  {dispest<-hfex(Bgamma%*%alphaG,nn=num)}
    # dispest<-hfratio(Bgamma%*%alphaG,nn=num)}
    res<-NULL
    ntrial<-NULL
    ntot<-0
    alphaMinit<-alphaM
    alphaGinit<-alphaG
    for(j in 1:maxitT){
        mm<-BinomialMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lambdaM,tol=tol,maxit=maxitM)
        alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2);alphaM<-mm$alphaM
        meanest<-as.vector(exp(Bmu%*%mm$alphaM)/(1+exp(Bmu%*%mm$alphaM)))
        zz<-BinomialDisp(y=y,num=num,Bgamma=Bgamma,lambdaG=lambdaG,truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,linktype=linktype,maxit=maxitG)
        if(linktype=="exp") {dispest<-exp(Bgamma%*%zz$alphaG)}
        if(linktype=="BB")  {dispest<-hfex(Bgamma%*%zz$alphaG,nn=num)}
        # dispest<-hfratio(Bgamma%*%zz$alphaG,nn=num)}
        alphaGdiff<-sum(((zz$alphaG-alphaG)/(alphaG))^2);alphaG<-zz$alphaG
        alphaGinit<-alphaG;alphaMinit<-alphaM
        if(max(alphaGdiff,alphaMdiff)<tolSTOP) {
            conv <- TRUE
            break
        }
    }
    
    # fill the empty result list
    res$alphaM<-as.vector(mm$alphaM)
    res$alphaG<-as.vector(zz$alphaG)
    res$meanest <- as.vector(num*exp(Bmu%*%mm$alphaM)/(1+exp(Bmu%*%mm$alphaM)))
    res$dispest <- as.vector(exp(Bgamma%*%zz$alphaG)^{-1})
    res$numberOfIteration<-j-1
    res$conv<-conv
    res
}




b1sGCV<-function(lambdas,y,num,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
BinomialBoth(y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,lambdaM=lambdas[1],lambdaG=lambdas[2],penordM=penordM,penordG=penordG,tol=tol,type=type,alphaM=alphaM,alphaG=alphaG,linktype=linktype,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT)$GCV
}
b1sGCV1<-function(lambdas,y,num,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
BinomialBoth(y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,lambdaM=lambdas[1],lambdaG=lambdas[2],penordM=penordM,penordG=penordG,tol=tol,type=type,alphaM=alphaM,alphaG=alphaG,linktype=linktype,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT)$GCV1
}
b1sAIC<-function(lambdas,y,num,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
BinomialBoth(y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,lambdaM=lambdas[1],lambdaG=lambdas[2],penordM=penordM,penordG=penordG,tol=tol,type=type,alphaM=alphaM,alphaG=alphaG,linktype=linktype,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT)$AIC
}
b1sAIC1<-function(lambdas,y,num,Bmu,Bgamma,alphaM=0,alphaG=0,maxitM=85,maxitG=85,maxitT=120,penordM=2,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
BinomialBoth(y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,lambdaM=lambdas[1],lambdaG=lambdas[2],penordM=penordM,penordG=penordG,tol=tol,type=type,alphaM=alphaM,alphaG=alphaG,linktype=linktype,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT)$AIC1
}



Binomial1s<-function(y,num,Bmu,Bgamma,alphaM=0,alphaG=0,lambdaM,lambdaG,maxitM=85,maxitG=95,maxitT=125,maxitOUT=7,penordM=2,penordG=2,tol=10^-10,type="quasi",crit="GCV1",minlM=0.005,minlG=0.005,maxlM=10000,maxlG=10000,linktype="exp"){
if(crit!="GCV"  & crit!="GCV1" & crit!="AIC" & crit!="AIC1"){
##### this chooses 'optimal' value for the smoothing parameters and fits the model 
##### the optimal values for lambdaM and lambdaG are chosen simultaneously via numerical minimization
warning("unidentified criterion, use GCV instead")
crit<-"GCV"}
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bmu))}
n<-length(y);meanest<-as.vector(eexpp(Bmu%*%alphaM))
if(linktype=="exp") {dispest<-exp(Bgamma%*%alphaG)}
if(linktype=="BB")  {dispest<-hfex(Bgamma%*%alphaG,nn=num)}
n<-length(y);res<-list();k<-ncol(Bmu);kt<-ncol(Bgamma);nit<-0
lambdaMinit<-lambdaM;lambdaGinit<-lambdaG
for(i in 1:maxitOUT){
nit<-nit+1;alphaMinit<-alphaM;alphaGinit<-alphaG
if(crit=="AIC"){
ee<-optim(c(lambdaM,lambdaG),b1sAIC,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),b1sAIC,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
if(crit=="AIC1"){
ee<-optim(c(lambdaM,lambdaG),b1sAIC1,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),b1sAIC1,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
if(crit=="GCV"){
ee<-optim(c(lambdaM,lambdaG),b1sGCV,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),b1sGCV,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
if(crit=="GCV1"){
ee<-optim(c(lambdaM,lambdaG),b1sGCV1,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
if(ee$conv > 0.5){
warning(paste("at iteration",i,"optim did not converge, I restart from initial lambda"))
ee<-optim(c(lambdaMinit,lambdaGinit),b1sGCV1,y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,alphaM=alphaM,alphaG=alphaG,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT,penordM=penordM,penordG=penordG,tol=tol,type=type,linktype=linktype,method="L-BFGS-B",upper=c(maxlM,maxlG),lower=c(minlM,minlG))
}}
lambdaM<-ee$par[1];lambdaG<-ee$par[2]
qq<-BinomialBoth(y=y,num=num,Bmu=Bmu,Bgamma=Bgamma,lambdaM=lambdaM,lambdaG=lambdaG,penordM=penordM,penordG=penordG,tol=tol,type=type,alphaM=alphaMinit,alphaG=alphaGinit,linktype=linktype,maxitM=maxitM,maxitG=maxitG,maxitT=maxitT)
alphaM<-qq$alphaM;alphaG<-qq$alphaG
meanest<-as.vector(exp(Bmu%*%qq$alphaM)/(1+exp(Bmu%*%qq$alphaM)))
if(linktype=="exp") {dispest<-exp(Bgamma%*%qq$alphaG)}
if(linktype=="BB")  {dispest<-hfex(Bgamma%*%qq$alphaG,nn=num)}
if (sum((alphaGinit-alphaG)^2)/sum(alphaG^2)<=tol & sum((alphaMinit-alphaM)^2)/sum(alphaM^2)<=tol  ) break
if (i == maxitOUT) warning("Iteration limit reached without convergence")
}
res$lambdaM<-lambdaM;res$lambdaG<-lambdaG
res$vals<-c(lambdaM,lambdaG);res$ntri<-nit
res$meanest<-as.vector(meanest);res$dispest<-as.vector(dispest)
res$alphaM<-as.vector(qq$alphaM);res$alphaG<-as.vector(qq$alphaG)
res$dfM<-qq$dfM;res$dfG<-qq$dfG
res$Dev<-dev<-qq$Dev;res$linktype<-linktype
res$AIC<-qq$AIC;res$AIC1<-qq$AIC1
res$GCV<-qq$GCV;res$GCV1<-qq$GCV1
res
}






############################
#### two steps procedures


BinomialMean<-function(y,Bmu,thetatrue,alphaM=0,lambdaM,maxit=950,penordM=2,tol=10^-10){
#### fit, for a given lambdaM and a given dispersion function, the model for the mean
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
res<-list();k<-ncol(Bmu)
DistMu<-diag(k)
if(penordM>0) DistMu<-diff(DistMu,diff=penordM)
Pmu<-t(DistMu)%*%DistMu
nit<-0
for(i in 1:maxit){
nit<-nit+1
# W<-diag(c((1/thetatrue)*exp(Bmu%*%alphaM)/((1+exp(Bmu%*%alphaM))^2)))
# z<-Bmu%*%alphaM+(y-(exp(Bmu%*%alphaM)/(1+exp(Bmu%*%alphaM))))
# matr<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
# alphaM1<-matr%*%z
W<-diag(as.vector(1/thetatrue))
alphaM1<-alphaM+ginv(t(Bmu)%*%diag(as.vector((1/thetatrue)*exp(Bmu%*%alphaM)/((1+exp(Bmu%*%alphaM))^2)))%*%Bmu+lambdaM*Pmu)%*%(t(Bmu)%*%(W)%*%(y-(exp(Bmu%*%alphaM)/(1+exp(Bmu%*%alphaM))))-lambdaM*Pmu%*%alphaM)
alphaMeps<-sum((alphaM-alphaM1)^2)/sum(alphaM^2);alphaM<-alphaM1
if (alphaMeps<= tol) {
break }
if (i == maxit) warning("Iteration limit reached without convergence in mean parameter")
}
W<-diag(as.vector(exp(Bmu%*%alphaM)/(thetatrue*(1+exp(Bmu%*%alphaM))^2)))
matr<-ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
# res$hmean<-Bmu%*%ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W
# res$hmean<-Bmu%*%matr
res$df<-sum(diag(matr%*%Bmu))
res$alphaM<-alphaM
res$numberOfIteration<-nit
# res$df<-sum(diag(ginv(t(Bmu)%*%W%*%Bmu+lambdaM*Pmu)%*%t(Bmu)%*%W%*%Bmu))
res
}



BinomialMeanGCV<-function(y,Bmu,lambdasM,thetatrue,penordM=2,tol=10^-10,alphaM=0){
### computes the GCV for a certain grid of lambda values, function needed for grid search
ddf<-NULL;GCV2<-GCV<-NULL;n<-length(y);co<-NULL;lambdaseq<-NULL;rr<-list()
for( i in 1:length(lambdasM)){
lambdaM<-lambdasM[i]
pp<-BinomialMean(y=y,Bmu=Bmu,thetatrue=thetatrue,alphaM=alphaM,lambdaM=lambdaM,penordM=penordM,tol=tol);df<-pp$df
eta<-Bmu%*%pp$alphaM
Dev<-2*(y*(log(y/(1-y)))-y*(eta)-log(1/(1-y))+log(1+exp(eta)));seq0<-seq(1,n)[y==0];seq1<-seq(1,n)[y==1]
Dev[seq0]<-2*log(1+exp(eta[seq0]));Dev[seq1]<-2*(log(1+exp(eta[seq1]))-(eta[seq1]))
GCV[i]<-sum(Dev/thetatrue)/((n-df)^2)
GCV2[i]<-sum(Dev)/((n-df)^2)
ddf<-c(ddf,df);co<-c(co,pp$conv)
lambdaseq<-c(lambdaseq,lambdaM)}
rr$GCV<-GCV;rr$conv<-co;rr$lambdaseq<-lambdaseq;rr$df<-ddf;rr$GCV2<-GCV2
rr
}

BinomialMeanGCVl<-function(lambdaMl,y,Bmu,thetatrue,penordM=2,tol=10^-10,alphaM=0){
### computes the GCV for a certain grid of log(lambda) values
### needed for the numeric minimization
ddf<-NULL;GCV<-NULL;n<-length(y);co<-NULL;lambdaseq<-NULL;rr<-list()
pp<-BinomialMean(y=y,Bmu=Bmu,thetatrue=thetatrue,alphaM=alphaM,lambdaM=exp(lambdaMl),penordM=penordM,tol=tol);df<-pp$df
eta<-Bmu%*%pp$alphaM
Dev<-2*(y*(log(y/(1-y)))-y*(eta)-log(1/(1-y))+log(1+exp(eta)));seq0<-seq(1,n)[y==0];seq1<-seq(1,n)[y==1]
Dev[seq0]<-2*log(1+exp(eta[seq0]));Dev[seq1]<-2*(log(1+exp(eta[seq1]))-(eta[seq1]))
GCV<-sum(Dev/thetatrue)/((n-df)^2)
GCV}




BinomialDisp<-function(y,num,Bgamma,truemean,alphaG=0,lambdaG,maxit=50,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
#### fit, for a given lambdaG and a given mean function, the model for the dispersion
n<-length(y);res<-list();kt<-ncol(Bgamma)
if (length(alphaG)==1){
if(alphaG==0) alphaG<-as.vector(rep(0.01,l=ncol(Bgamma)))
else alphaG<-as.vector(rep(alphaG,l=ncol(Bgamma)))}
# if(linktype=="BB") {lambdaG<-10000000+10*lambdaG}
DistGam<-diag(kt);nit<-0
if(penordG>0) DistGam<-diff(DistGam,diff=penordG)
Pgam<-t(DistGam)%*%DistGam;eta<-log(truemean/(1-truemean));pitrue<-truemean
if (type=="quasi"){Dev<-DevBinE(y,eta)}
if(type=="pseudo"){
Dev<-((y-truemean)^2)/(truemean*(1-truemean))}
if(type != "pseudo" & type != "quasi" ){
stop("the type of likelihhod to use is not recognized")
}
for(j in 1:maxit){
nit<-nit+1
if(linktype=="BB"){
q<-(Bgamma%*%alphaG+(num*Dev-hfex((Bgamma%*%alphaG),nn=num))/(hfpex(Bgamma%*%alphaG,nn=num)))
Wg<-diag(0.5*as.vector((hfratioex((Bgamma%*%alphaG),nn=num))^2))}
if(linktype=="exp"){
q<-(Bgamma%*%alphaG+(num*Dev-exp(Bgamma%*%alphaG))/(exp(Bgamma%*%alphaG)));Wg<-diag(0.5,length(Dev))}
alphaG1<-ginv(t(Bgamma)%*%Wg%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%Wg%*%q
alphaGeps<-sum(((alphaG1-alphaG)/alphaG1)^2)
if (j == 1) alphaGeps<-1
alphaG<-alphaG1 #alphaGeps<-10^{-j*3}
if (alphaGeps <= tol) break}
rm(q,Wg)
if(linktype=="BB"){
Wg<-diag(0.5*(as.vector((hfratioex(Bgamma%*%alphaG,nn=num))^2)))}
if(linktype=="exp"){Wg<-diag(0.5,length(Dev))}
# 
# res$hdisp<-Bgamma%*%ginv(t(Bgamma)%*%Wg%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%Wg
res$alphaG<-as.vector(alphaG)
res$numberOfIteration<-nit;res$deviance<-num*Dev
df<-res$df<-sum(diag(ginv(t(Bgamma)%*%Wg%*%Bgamma+lambdaG*Pgam)%*%t(Bgamma)%*%Wg%*%Bgamma))
res
}




BinomialDispGCV<-function(y,num,Bgamma,lambdasG,truemean,alphaG,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
### computes the GCV for a certain grid of lambda values, function needed for grid search
ddf<-NULL;GCV<-NULL;n<-length(y);co<-NULL;lambdaseq<-NULL;rr<-list()
for( i in 1:length(lambdasG)){
lambdaG<-lambdasG[i];n<-length(y)
pp<-BinomialDisp(y=y,num=num,Bgamma=Bgamma,truemean=truemean,alphaG=alphaG,lambdaG=lambdaG,penordG=penordG,tol=tol,type=type,linktype=linktype);df<-pp$df#;eta<-log(truemean/(1-truemean))
if(linktype=="exp") {gest<-exp(Bgamma%*%pp$alphaG)}
if(linktype=="BB")  {gest<-hfex(Bgamma%*%pp$alphaG,nn=num)}# gest<-hf(Bgamma%*%pp$alphaG,nn=num)
rres<-pp$dev;devDev<- (rres/gest)+log(gest/rres)-1
GCV[i]<-n*sum(devDev)/((n-df)^2)
ddf<-c(ddf,df);co<-c(co,pp$conv);lambdaseq<-c(lambdaseq,lambdaG)}
rr$GCV<-GCV;rr$conv<-co;rr$lambdaseq<-lambdaseq;rr$df<-ddf
rr
}



BinomialDispGCVl<-function(lambdaGl,y,num,Bgamma,truemean,alphaG,penordG=2,tol=10^-10,type="quasi",linktype="exp"){
### computes the GCV for a certain grid of log(lambda) values
### needed for numeric minimization
ddf<-NULL;GCV<-NULL;n<-length(y);co<-NULL;lambdaseq<-NULL;rr<-list()
pp<-BinomialDisp(y=y,num=num,Bgamma=Bgamma,truemean=truemean,alphaG=alphaG,lambdaG=exp(lambdaGl),penordG=penordG,tol=tol,type=type,linktype=linktype);df<-pp$df#;eta<-log(truemean/(1-truemean))
rres<-pp$dev
if(linktype=="exp") {gest<-exp(Bgamma%*%pp$alphaG)}
if(linktype=="BB") {gest<-hfex(Bgamma%*%pp$alphaG,nn=num)}
# gest<-hfratio(Bgamma%*%pp$alphaG,nn=num)}
devDev<- (rres/gest)+log(gest/rres)-1
GCV<-n*sum(devDev)/((n-df)^2)
GCV}






BinomialBothGCV<-function(y,num,Bmu,Bgamma,seqM,seqG,lenM,lenG,penordM=2,penordG=2,tol=10^-10,type="quasi",plotM=FALSE,plotG=FALSE,minlM=0.005,maxlM=8000,minlG=0.005,maxlG=8000,uppars=TRUE,alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp"){
##### this chooses 'optimal' value for the smoothing parametersin two steps and fits the model 
##### the smoothing parameter choice is done via a grid search (the same as in the normal model)
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bmu))}
llM<-seq(min(seqM),max(seqM),length=lenM);llG<-seq(min(seqG),max(seqG),length=lenG)
diffsMean<-NULL;diffsDisp<-NULL
n<-length(y);dispest<-as.vector(rep(1,n));meanest<-as.vector(rep(mean(y),n))
gcvMmins<-NULL;gcvGmins<-NULL;lambdaM<-10^-5;lambdaG<-10^-5
res<-NULL;ntrial<-NULL;ntot<-0;numsM<-NULL;numsG<-NULL
lambdaM<-NULL;lambdaG<-NULL;extremesM<-NULL;extremesG<-NULL
for(j in 2:8){
ntot<-ntot+1
checkM<-NULL;checkM[1]<-10^-5;checkG<-NULL;checkG[1]<-10^-5
for(i in 1:5){
alphaMinit<-as.vector(alphaM);alphaGinit<-as.vector(alphaG)
gcv<-BinomialMeanGCV(y=y,Bmu=Bmu,lambdasM=llM,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit)
gcvM<-min(gcv$GCV)
if (plotM==TRUE){
plot(llM,gcv$GCV,type="l")}
ss<-min(seq(1,lenM)[gcv$GCV==min(gcv$GCV)]);lM<-llM[ss];checkM[i+1]<-lM#min(gcv$GCV)
if(checkM[i]/checkM[i+1]<1.05 & checkM[i]/checkM[i+1]>0.95){
llM<-sort(seq(lM*0.5,lM*2,length=lenM))
extremesM<-c(extremesM,as.character(ss))
break}
else{
paceM<-llM[2]-llM[1];lM<-llM[ss]
if(ss==lenM){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL
diffGCVFin<-NULL;eps<-min(gcv$GCV)*10^-3;sum5<-NULL
for( k in 1:(lenM-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenM])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))}
sum5<-diffsTOT;ns<-min(ss,seq(1,length(llM))[sum5>5.5])
rm(sum5)
if(ns != ss){
lM<-llM[ns]
extremesM<-c(extremesM,"max - cutted")
llM<-sort(seq(from=max(minlM,min(llM[1]*0.8,lM-6*paceM)),by=1.6*paceM,l=lenM))
break}
if(ns==ss){
lM<-llM[ss]
llM<- sort(seq(from=max(minlM,lM-10*paceM),by=1.2*paceM,l=lenM))
#llM<-sort(seq(from=min(llM[1]*.9,lM-7*paceM),by=1.2*paceM,l=lenM))
extremesM<-c(extremesM,"max")}}
if(ss==1){
if(lM==minlM){
llM<-sort(seq(minlM,by=1.2*paceM,l=lenM))}
else{
llM<-sort(seq(max(minlM,lM-(0.8*lenM*paceM)),by=1.2*paceM,l=lenM))}
extremesM<-c(extremesM,"min")}
if(ss!=1 & ss!=lenM){
llM<-sort(seq(from=max(minlM,lM-6*paceM),to=max(minlM+10*paceM,lM+6*paceM),l=lenM))
extremesM<-c(extremesM,as.character(ss))
}}}
extremesM<-c(extremesM," ")
if(lM>maxlM){
lM<-maxlM
warning("lM higher then maxlM forced to be maxlM")
}
numsM<-c(numsM,i)
mm<-BinomialMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lM)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2);alphaM<-mm$alphaM
meanest<-as.vector(exp(Bmu%*%mm$alphaM)/(1+exp(Bmu%*%mm$alphaM)))
rm(gcv,ss)
for(i in 1:5){
gcv<-BinomialDispGCV(y=y,num=num,Bgamma=Bgamma,lambdasG=llG,truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,linktype=linktype);gcvG<-min(gcv$GCV)
if(plotG==TRUE){plot(llG,gcv$GCV,type="l")}
ss<-min(seq(1,lenG)[gcv$GCV==min(gcv$GCV)]);lG<-llG[ss];checkG[i+1]<-lG#min(gcv$GCV)
if(checkG[i]/checkG[i+1]<1.05 & checkG[i]/checkG[i+1]>.95){
llG<-sort(seq(0.5*lG,lG*2,length=lenG))
extremesG<-c(extremesG,as.character(ss))
break}
else{
paceG<-llG[2]-llG[1]
if(ss==lenG){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL
diffGCVFin<-NULL
eps<-min(gcv$GCV)*10^-2
for( k in 1:(lenG-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenG])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))
}
sum6<-diffsTOT#sum(diffGCV1<eps,diffGCV2<eps,diffGCV3<eps,diffGCV4<eps,diffGCV5<eps)
ns<-min(ss,seq(1,length(llG))[sum6>5.5]);rm(sum6)
if(ns != ss){
lG<-llG[ns];extremesG<-c(extremesG,"max - cutted")
llG<-sort(seq(from=max(minlG,min(llG[1]*0.8,lG-6*paceG)),by=1.6*paceG,l=lenG))
break
}
if(ns==ss){
lG<-llG[ss];llG<-sort(seq(from=max(minlG,lG-10*paceG),by=1.2*paceG,l=lenG));extremesG<-c(extremesG,"max")
}}
if(ss==1){
if(lG==0){
llG<-sort(seq(max(minlG,lG-(0.8*lenG*paceG)),by=1.2*paceG,l=lenG))
}
else{
llG<-sort(seq(from=max(minlG,lG-(0.8*lenG*paceG)),by=1.2*paceG,l=lenG))}
extremesG<-c(extremesG,"min")}
if(ss!=1 & ss!=lenG){
llG<-sort(seq(from=max(minlG,lG-6*paceG),to=max(minlG+10*paceG,lG+6*paceG),l=lenG))
extremesG<-c(extremesG,as.character(ss))
}}}
extremesG<-c(extremesG," ")
if(lG>maxlG){
lG<-maxlG;warning("lG higher then maxlG forced to be maxlG")}
numsG<-c(numsG,i)
zz<-BinomialDisp(y=y,num=num,Bgamma=Bgamma,truemean=meanest,penordG=penordG,alphaG=as.vector(alphaGinit),lambdaG=lG,type=type,linktype=linktype)
alphaGdiff<-sum(((zz$alphaG-alphaG)/(alphaG))^2);alphaG<-zz$alphaG
if(linktype=="exp") {dispest<-exp(Bgamma%*%zz$alphaG)}
if(linktype=="BB") {dispest<-hfex(Bgamma%*%zz$alphaG,nn=num)}
# dispest<-hfratio(Bgamma%*%zz$alphaG,nn=num)}
lambdaM<-c(lambdaM,lM);lambdaG<-c(lambdaG,lG)
if (uppars == TRUE){
alphaGinit<-alphaG;alphaMinit<-alphaM}
if(max(alphaGdiff,alphaMdiff)<tolSTOP) break}
lval<-c(lM,lG)
res$vals<-c(lM,lG);res$ntri<-j-1
res$sequencesM<-extremesM;res$sequencesG<-extremesG
res$meanest<-as.vector(meanest);res$dispest<-as.vector(dispest)
res$lambdasM<-lambdaM[1:length(lambdaM)];res$lambdasG<-lambdaG[1:length(lambdaG)]
res$alphaM<-as.vector(mm$alphaM);res$alphaG<-as.vector(zz$alphaG)
res$dfM<-mm$df;res$dfG<-zz$df;res$vecNumM<-numsM;res$vecNumG<-numsG
res$Dev<-zz$deviance;res$linktype<-linktype
res
}








BinomialBothGCVnum<-function(y,num,Bmu,Bgamma,lambdaM,lambdaG,penordM=2,penordG=2,tol=10^-10,type="quasi",minlM=0.005,maxlM=8000,minlG=0.005,maxlG=8000,uppars=TRUE,alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp"){

##### this chooses 'optimal' value for the smoothing parametersin two steps and fits the model 
##### this chooses 'optimal' value for the smoothing parametersin via numerical minimization
##### the smoothing parameter choice is done via a grid search (the same as in the normal model)
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))
}
if (length(alphaG)==1){
if(alphaG==0) alphaG<-rep(0.1,l=ncol(Bgamma))
else alphaG<-rep(alphaG,l=ncol(Bmu))
}
n<-length(y);meanest<-as.vector(Bmu%*%alphaM)
if(linktype=="exp") {dispest<-exp(Bgamma%*%alphaG)}
if(linktype=="BB")  {dispest<-hfex(Bgamma%*%alphaG,nn=num)}
# dispest<-hfratio(Bgamma%*%alphaG,nn=num)}
res<-NULL;ntrial<-NULL;ntot<-0
lambdaM<-log(lambdaM);lambdaG<-log(lambdaG)
alphaMinit<-alphaM;alphaGinit<-alphaG
for(j in 2:8){
lambdaM<-optim(lambdaM,BinomialMeanGCVl,y=y,Bmu=Bmu,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit,upper=maxlM,lower=minlM,method="L-BFGS-B")$par
mm<-BinomialMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=exp(lambdaM),tol=tol)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2);alphaM<-mm$alphaM
meanest<-as.vector(exp(Bmu%*%mm$alphaM)/(1+exp(Bmu%*%mm$alphaM)))
lambdaG<-optim(lambdaG,BinomialDispGCVl,y=y,num=num,Bgamma=Bgamma,truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,linktype=linktype,upper=maxlG,lower=minlG,method="L-BFGS-B")$par
zz<-BinomialDisp(y=y,num=num,Bgamma=Bgamma,lambdaG=exp(lambdaG),truemean=meanest,penordG=penordG,tol=tol,type=type,alphaG=alphaGinit,linktype=linktype)
if(linktype=="exp") {dispest<-exp(Bgamma%*%zz$alphaG)}
if(linktype=="BB")  {dispest<-hfex(Bgamma%*%zz$alphaG,nn=num)}
# dispest<-hfratio(Bgamma%*%zz$alphaG,nn=num)}
alphaGdiff<-sum(((zz$alphaG-alphaG)/(alphaG))^2);alphaG<-zz$alphaG
if (uppars == TRUE){
alphaGinit<-alphaG;alphaMinit<-alphaM}
if(max(alphaGdiff,alphaMdiff)<tolSTOP) break
}
res$vals<-c(lambdaM,lambdaG);res$ntri<-j-1
res$meanest<-as.vector(meanest);res$dispest<-as.vector(dispest)
res$alphaM<-as.vector(mm$alphaM);res$alphaG<-as.vector(zz$alphaG)
res$dfM<-mm$df;res$dfG<-zz$df
res$Dev<-zz$deviance;res$linktype<-linktype
res
}




BinomialMeanOnlyGCV<-function(y,num,Bmu,seqM,lenM,penordM=2,tol=10^-10,type="quasi",plotM=FALSE,minlM=0.005,maxlM=8000,uppars=TRUE,alphaM=0,tolSTOP=10^-2,ESTdisp=TRUE){
if (length(alphaM)==1){
if(alphaM==0) alphaM<-rep(log(mean(y)/(1-mean(y))),l=ncol(Bmu))
else alphaM<-rep(alphaM,l=ncol(Bmu))}
llM<-seq(min(seqM),max(seqM),length=lenM)
diffsMean<-NULL
n<-length(y)
dispest<-as.vector(rep(1,n));meanest<-as.vector(rep(mean(y),n))
gcvMmins<-NULL;lambdaM<-10^-5
res<-NULL;ntrial<-NULL;ntot<-0;numsM<-NULL
lambdaM<-NULL;lambdaG<-NULL;extremesM<-NULL
for(j in 2:8){
ntot<-ntot+1
checkM<-NULL;checkM[1]<-10^-5
for(i in 1:5){
alphaMinit<-as.vector(alphaM)
gcv<-BinomialMeanGCV(y=y,Bmu=Bmu,lambdasM=llM,thetatrue=dispest,penordM=penordM,tol=tol,alphaM=alphaMinit)
gcvM<-min(gcv$GCV)
if (plotM==TRUE){
plot(llM,gcv$GCV,type="l")}
ss<-min(seq(1,lenM)[gcv$GCV==min(gcv$GCV)]);lM<-llM[ss];checkM[i+1]<-lM#min(gcv$GCV)
if(checkM[i]/checkM[i+1]<1.05 & checkM[i]/checkM[i+1]>0.95){
llM<-sort(seq(lM*0.5,lM*2,length=lenM))
extremesM<-c(extremesM,as.character(ss))
break}
else{
paceM<-llM[2]-llM[1];lM<-llM[ss]
if(ss==lenM){
diffGCV1<-NULL;diffGCV2<-NULL;diffGCV3<-NULL
diffGCV4<-NULL;diffGCV5<-NULL;diffsTOT<-NULL
diffGCVFin<-NULL;eps<-min(gcv$GCV)*10^-3;sum5<-NULL
for( k in 1:(lenM-5)){
diffGCV1[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+1])<eps)
diffGCV2[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+2])<eps)
diffGCV3[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+3])<eps)
diffGCV4[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+4])<eps)
diffGCV5[k]<-as.integer((gcv$GCV[k]-gcv$GCV[k+5])<eps)
diffGCVFin[k]<-as.integer((gcv$GCV[k]-gcv$GCV[lenM])<eps)
diffsTOT[k]<-as.integer(sum(diffGCV1[k],diffGCV2[k],diffGCV3[k],diffGCV4[k],diffGCV5[k],diffGCVFin[k]))}
sum5<-diffsTOT;ns<-min(ss,seq(1,length(llM))[sum5>5.5])
rm(sum5)
if(ns != ss){
lM<-llM[ns]
extremesM<-c(extremesM,"max - cutted")
llM<-sort(seq(from=max(minlM,min(llM[1]*0.8,lM-6*paceM)),by=1.6*paceM,l=lenM))
break}
if(ns==ss){
lM<-llM[ss]
llM<- sort(seq(from=max(minlM,lM-10*paceM),by=1.2*paceM,l=lenM))
#llM<-sort(seq(from=min(llM[1]*.9,lM-7*paceM),by=1.2*paceM,l=lenM))
extremesM<-c(extremesM,"max")}}
if(ss==1){
if(lM==minlM){
llM<-sort(seq(minlM,by=1.2*paceM,l=lenM))}
else{
llM<-sort(seq(max(minlM,lM-(0.8*lenM*paceM)),by=1.2*paceM,l=lenM))}
extremesM<-c(extremesM,"min")}
if(ss!=1 & ss!=lenM){
llM<-sort(seq(from=max(minlM,lM-6*paceM),to=max(minlM+10*paceM,lM+6*paceM),l=lenM))
extremesM<-c(extremesM,as.character(ss))
}}}
extremesM<-c(extremesM," ")
if(lM>maxlM){
lM<-maxlM
warning("lM higher then maxlM forced to be maxlM")
}
numsM<-c(numsM,i)
mm<-BinomialMean(y=y,Bmu=Bmu,thetatrue=dispest,alphaM=alphaMinit,lambdaM=lM)
alphaMdiff<-sum(((mm$alphaM-alphaM)/(alphaM))^2);alphaM<-mm$alphaM
meanest<-as.vector(exp(Bmu%*%mm$alphaM)/(1+exp(Bmu%*%mm$alphaM)))
rm(gcv,ss)
if(ESTdisp) {dispest<-rep(sum(num*DevBinP(y,meanest)/(n-mm$df)),l=n)}
if(!ESTdisp) j<-7
lambdaM<-c(lambdaM,lM)
if (uppars == TRUE){alphaMinit<-alphaM}
if(max(alphaMdiff)<tolSTOP) break}
lval<-c(lM)
res$vals<-c(lM);res$ntri<-j-1
res$sequencesM<-extremesM
res$meanest<-as.vector(meanest);res$dispest<-as.vector(dispest)
res$lambdasM<-lambdaM[1:length(lambdaM)]
res$alphaM<-as.vector(mm$alphaM)
res$dfM<-mm$df;res$vecNumM<-numsM
res$Dev<-DevBinP(y,meanest)
res
}


# # ### a simulated example....
# # set.seed(16)
# # x<-sort(runif(250))
# # eta<-sin(5*x)
# # pp<-exp(eta)/(1+exp(eta))
# # xi<-0.7+2.5*x^2*cos(5*x)*sin(4*x)
# # gg<-exp(xi)
# # set.seed(17) 
# # yy<-ranBin(250,p=pp,n=280,theta=gg)
# # prop<-yy$sample/yy$n;ntrial<-yy$n
# # 
# # 
# # #### Build the matrices
# # dx<-(max(x)-min(x))/37;kn<-seq(min(x)-3*dx,max(x)+3*dx,by=dx)
# # basem<-splineDesign(knots=kn,x=x,ord=(3+1),0*x,outer=TRUE)
# # dx<-(max(x)-min(x))/31;kn<-seq(min(x)-3*dx,max(x)+3*dx,by=dx)
# # based<-splineDesign(knots=kn,x=x,ord=(3+1),0*x,outer=TRUE)
# # 
# # 
# # 
# # bb<-BinomialBothGCV(y=prop,num=ntrial,Bmu=basem,Bgamma=based,seqM=c(0.1,150),seqG=c(1,1500),lenM=45,lenG=45,penordM=2,penordG=2,tol=10^-10,type="quasi",plotM=TRUE,plotG=TRUE,minlM=0.005,maxlM=8000,minlG=0.005,maxlG=8000,uppars=TRUE,alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp")
# # bb2<-BinomialBothGCV(y=prop,num=ntrial,Bmu=basem,Bgamma=based,seqM=c(0.1,150),seqG=c(1,1500),lenM=45,lenG=45,penordM=2,penordG=2,tol=10^-10,type="quasi",plotM=TRUE,plotG=TRUE,minlM=0.005,maxlM=8000,minlG=0.005,maxlG=8000,uppars=TRUE,alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="BB")
# # bb3<-BinomialBoth(y=prop,num=ntrial,Bmu=basem,Bgamma=based,lambdaM=1.8,lambdaG=5)
# # 
# # par(mfrow=c(1,2))
# # plot(x,prop)
# # lines(x,pp)
# # lines(x,eexpp(basem%*%bb$alphaM),col=2)
# # lines(x,eexpp(basem%*%bb2$alphaM),col=4)
# # lines(x,eexpp(basem%*%bb3$alphaM),col=3)
# # 
# # plot(x,bb$Dev)
# # lines(x,gg)
# # lines(x,bb$dispest,col=2)
# # lines(x,bb2$dispest,col=4)
# # lines(x,bb3$dispest,col=3)

