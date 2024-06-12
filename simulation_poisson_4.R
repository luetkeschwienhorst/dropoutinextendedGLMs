
# libraries
library(splines)
library(MASS)
library(lattice)
library(gridBase)
library(grid)
library(caret)
library(parallel)
library(mgcv)
library(viridis)
library(scales)
library(pryr)
library(gamlss)
library(inline)
library(foreach)
library(doParallel)

# functions
source("functions.R")
source("GijbelsProsdocimiClaeskens.R")

# simulation setup
cores <- detectCores()
registerDoParallel(cores)
data_distr <- "poisson"
phi <- 1
interval_p_mu <- c(0,1)
interval_p_gamma <- c(0,1)
interval_sd_mu <- c(0,3)
interval_sd_gamma <- c(0,6)
m_mu <- 30
m_gamma <- 20
cv_draws <- 500
samplesizes <- c(250,500,1000)
replicates <- 101

X <- list()
set.seed(12345)
for(i in 1:length(samplesizes)) {
    x <- matrix(data=sort(runif(samplesizes[i])),nrow=samplesizes[i])
    for (j in 2:replicates) {
        x <- cbind(x,sort(runif(samplesizes[i])))
    }
    X[[i]] <- x
}

# simulation
i <- 4
for (j in 1:3) {
    print(paste("Model: model_",sep="",i,j,"_poisson"))
    model <- paste("model_",sep="",i,j,"_poisson")
    assign(model,simulation_def(data_distribution=data_distr,phi=phi,nu=1,N=NULL,
                                b=b_poisson,Db=Db_poisson,Db_inv=Db_inv_poisson,D2b=D2b_poisson,X=X,
                                truemean=testfunctions_mean[[i]],truemean_index=i,truedisp=testfunctions_disp[[j]],replicates=replicates,
                                samplesizes=samplesizes,interval_p_mu=interval_p_mu,interval_p_gamma=interval_p_gamma,
                                interval_sd_mu=interval_sd_mu,interval_sd_gamma=interval_sd_gamma,
                                cv_draws=cv_draws,cv_folds=5,learning_rate=list(method="ADADELTA",rho=0.95),
                                batch_size=30,m_mu=m_mu,m_gamma=m_gamma,cores=cores,progress=TRUE))
    image <- paste("model_",sep="",i,1,"_to_",i,j,"_poisson.RData")
    save.image(file=image)
    cat("\n\n")
}

