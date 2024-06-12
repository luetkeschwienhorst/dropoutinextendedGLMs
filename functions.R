

##### TESTFUNCTIONS #####

pdf_normal <- function(x,mean,sd) {
    dnorm(x, mean = mean, sd = sd)
}

cdf_normal <- function(x,mean,sd) {
    pnorm(x, mean = mean, sd = sd)
}

# f_2 <- function(x){3+0.5*pdf_normal(x,0.2,0.04)+0.6*pdf_normal(x,0.35,0.04)}

pdf_skew_normal <- function(x,mean,sd,alpha) {
    (2/sd)*pdf_normal((x-mean)/sd,0,1)*cdf_normal(alpha*(x-mean)/sd,0,1)
}

testfunctions_mean <- list(
    f_1 <- function(x){2*sin(4*pi*x)*pdf_normal(x,0.5,0.05)},
    f_2 <- function(x){3+0.5*pdf_normal(x,0.2,0.04)+0.6*pdf_normal(x,0.35,0.04)},
    f_3 <- function(x){40+10*sin(4*pi*x)*pdf_normal(x,0.5,0.05)},
    f_4 <- function(x){10+1*pdf_normal(x,0.2,0.03)+1.2*pdf_normal(x,0.35,0.03)}
)

testfunctions_disp <- list(
    g_1 <- function(x){exp(-0.1*pdf_normal(x,0.2,0.03)+0.3*sin(4*pi*(x+0.5))*pdf_normal(x,0.5,0.04))},
    g_2 <- function(x){exp(-0.08*pdf_normal(x,0.4,0.03)+0.08*pdf_normal(x,0.7,0.03))},
    g_3 <- function(x){exp(-0.2*pdf_skew_normal(x,0.8,0.15,-4))},
    g_4 <- function(x){exp(-0.08*pdf_normal(x,0.4,0.04)-0.08*pdf_normal(x,0.7,0.02))}
)

##### CGFs, their derivatives and inverses #####

b_normal <- function (x) {(1/2)*x^2}
Db_normal <- function (x) {x}
Db_inv_normal <- function (x) {x}
D2b_normal <- function (x) {1}

b_poisson <- function (x) {exp(x)}
Db_poisson <- function (x) {exp(x)}
Db_inv_poisson <- function (x) {log(x)}
D2b_poisson <- function (x) {exp(x)}

b_binomial <- function(x,n) {n*log(1+exp(x))}
Db_binomial <- function(x,n) {n*exp(x)/(1+exp(x))}
Db_inv_binomial <- function(x,n) {log(x/(n-x))}
D2b_binomial <- function (x,n) {n*exp(x)/(1+exp(x))^2}

##### ESTIMATION #####

trunc_poly_basis <- function(X,knots,degree){
    B_1 <- outer(X,c(0:degree),"^")
    B_2 <- outer(X,knots,">")*outer(X,knots,"-")^degree
    B <- cbind(B_1,B_2)
    B
}

loglikelihood_def <- function(y,B_mu,B_gamma,beta,alpha,phi,nu,b,Db_inv) {
    loglikelihood <- sum((1/2)*B_gamma%*%alpha
                         +(1/phi)*exp(B_gamma%*%alpha)*(y*B_mu%*%beta-b(B_mu%*%beta))*nu
                         +(1/phi)*(1-exp(B_gamma%*%alpha))*(y*Db_inv(y)-b(Db_inv(y)))*nu)
    loglikelihood
}

loglikelihood_def_gijbels <- function(y,B_mu,B_gamma,beta,alpha,phi,nu,b,Db_inv,N=NULL,correct_for_zeros=TRUE,correct_for_N=TRUE) {
    dev <- y*Db_inv(y)-b(Db_inv(y))
    if (correct_for_zeros) {
        dev[which(y==0)] <- 0
    }
    if (correct_for_N) {
        dev[which(y==N)] <- 0
    }
    loglikelihood <- sum((1/2)*B_gamma%*%alpha
                         +(1/phi)*exp(B_gamma%*%alpha)*(y*B_mu%*%beta-b(B_mu%*%beta))*nu
                         +(1/phi)*(1-exp(B_gamma%*%alpha))*(dev)*nu)
    loglikelihood
}

gradient_def <- function(y,B_mu,B_gamma,beta,alpha,phi,nu,b,Db,Db_inv,
                         dropout_matrix_mu=1,dropout_matrix_gamma=1) {
    B_mu_dropout <- B_mu*dropout_matrix_mu
    B_gamma_dropout <- B_gamma*dropout_matrix_gamma
    z <- (1/phi)*nu*exp(B_gamma_dropout%*%alpha)*(y-Db(B_mu_dropout%*%beta))
    gradient_beta <- t(B_mu_dropout)%*%z
    q <- rep(1/2,length(y))+(1/phi)*nu*exp(B_gamma_dropout%*%alpha)*(y*(B_mu_dropout%*%beta-Db_inv(y))-b(B_mu_dropout%*%beta)+b(Db_inv(y))) 
    gradient_alpha <- t(B_gamma_dropout)%*%q
    gradient <- list()
    gradient$gradient_beta <- gradient_beta
    gradient$gradient_alpha <- gradient_alpha
    gradient
}

gradient_def_gijbels <- function(y,B_mu,B_gamma,beta,alpha,phi,nu,b,Db,Db_inv,
                         dropout_matrix_mu=1,dropout_matrix_gamma=1,N=NULL,
                         correct_for_zeros=TRUE,correct_for_N=TRUE) {
    B_mu_dropout <- B_mu*dropout_matrix_mu
    B_gamma_dropout <- B_gamma*dropout_matrix_gamma
    z <- (1/phi)*nu*exp(B_gamma_dropout%*%alpha)*(y-Db(B_mu_dropout%*%beta))
    gradient_beta <- t(B_mu_dropout)%*%z
    devi <- y*(B_mu_dropout%*%beta-Db_inv(y))-b(B_mu_dropout%*%beta)+b(Db_inv(y))
    if (correct_for_zeros) {
        devi[which(y==0)] <- -b(B_mu_dropout%*%beta)[which(y==0)]
    }
    if (correct_for_N) {
        devi[which(y==N)] <- (y*(B_mu_dropout%*%beta)-b(B_mu_dropout%*%beta))[which(y==N)]
    }
    q <- rep(1/2,length(y))+(1/phi)*nu*exp(B_gamma_dropout%*%alpha)*devi
    gradient_alpha <- t(B_gamma_dropout)%*%q
    gradient <- list()
    gradient$gradient_beta <- gradient_beta
    gradient$gradient_alpha <- gradient_alpha
    gradient
}

dropout_matrices <- function(noise_distribution,dropout_mu,dropout_gamma,
                             batch_size,ncol_B_mu,ncol_B_gamma) {
    if (noise_distribution=="Normal") {
        if (isTRUE(is.finite(dropout_mu))) {
            dropout_matrix_mu <- matrix(data = rnorm(batch_size*ncol_B_mu,1,dropout_mu),nrow = batch_size,ncol = ncol_B_mu)
        } else {
            dropout_matrix_mu <- matrix(data=rep(1,batch_size*ncol_B_mu),nrow=batch_size,ncol=ncol_B_mu)
            warning("Dropout parameter for the mean is invalid, such that no dropout is performed in the mean!")
        }
        if (isTRUE(is.finite(dropout_gamma))) {
            dropout_matrix_gamma <- matrix(data = rnorm(batch_size*ncol_B_gamma,1,dropout_gamma),nrow = batch_size,ncol = ncol_B_gamma)
        } else {
            dropout_matrix_gamma <- matrix(data = rep(1,batch_size*ncol_B_gamma),nrow = batch_size,ncol = ncol_B_gamma)
            warning("Dropout parameter for the dispersion is invalid, such that no dropout is performed in the dispersion!")
        }
    } else if (noise_distribution=="Bernoulli") {
        if (isTRUE(is.finite(dropout_mu))) {
            dropout_matrix_mu <- matrix(data=rbinom(batch_size*ncol_B_mu,1,1-dropout_mu),nrow = batch_size,ncol = ncol_B_mu)
            if (isTRUE(dropout_mu!=1)) {
                dropout_matrix_mu <- (1-dropout_mu)^-1*dropout_matrix_mu
            }
        } else {
            dropout_matrix_mu <- matrix(data=rep(1,batch_size*ncol_B_mu),nrow=batch_size,ncol=ncol_B_mu)
            warning("Dropout parameter for the mean is invalid, such that no dropout is performed in the mean!")
        }
        if (isTRUE(is.finite(dropout_gamma))) {
            dropout_matrix_gamma <- matrix(data=rbinom(batch_size*ncol_B_gamma,1,1-dropout_gamma),nrow = batch_size,ncol = ncol_B_gamma)
            if (isTRUE(dropout_gamma!=1)) {
                dropout_matrix_gamma <- (1-dropout_gamma)^-1*dropout_matrix_gamma
            }
        } else {
            dropout_matrix_gamma <- matrix(data = rep(1,batch_size*ncol_B_gamma),nrow = batch_size,ncol = ncol_B_gamma)
            warning("Dropout parameter for the dispersion is invalid, such that no dropout is performed in the dispersion!")
        }
    }
    dropout_matrices <- list()
    dropout_matrices$dropout_matrix_mu <- dropout_matrix_mu
    dropout_matrices$dropout_matrix_gamma <- dropout_matrix_gamma
    dropout_matrices
}

sga_stopping_criterion <- function(epsilon,values_loglikelihood) {
    mean_10 <- mean(values_loglikelihood[0.9*length(values_loglikelihood):length(values_loglikelihood)])
    mean_5 <- mean(values_loglikelihood[0.95*length(values_loglikelihood):length(values_loglikelihood)])
    diff <- abs(mean_10-mean_5)
    fulfilled <- FALSE
    if (isTRUE(is.finite(diff))) {
        if (diff<epsilon) {
            fulfilled <- TRUE
        }
    }
    list(fulfilled,diff)
}

# sga for any def with dropout regularization
sga_def_dropout <- function(y,B_mu,B_gamma,beta=0,alpha=0,phi=1,nu=1,N=NULL,b,Db,Db_inv,D2b,correct_for_zeros=FALSE,correct_for_N=FALSE,
                            noise_distribution="Bernoulli",dropout_mu=0,dropout_gamma=0,learning_rate=list(method="ADADELTA",rho=0.95),
                            batch_size=30,max_iterations=50000,tolerance=10^{-5},progress=TRUE,time=TRUE,
                            save_stopping_crit=TRUE) {
    
    # preparations
    start_time <- Sys.time()
    n <- length(y)
    k <- ncol(B_mu)
    l <- ncol(B_gamma)
    if (length(beta)==1) {
        if (beta==0) {
            #beta <- ginv(t(B_mu)%*%B_mu)%*%t(B_mu)%*%y
            beta <- c(Db_inv(mean(y)),rep(0,k-1))
        } else {
            beta <- rep(beta,l=k)
        }
    }
    if (length(alpha)==1) {
        if(alpha==0) {
            alpha <- rep(0,l=l)
        } else {
            alpha <- rep(alpha,l=l)
        }
    }
    if (length(nu)==1) {
        nu <- rep(1,length(y))
    }
    if (save_stopping_crit) {
        values_stopping_crit <- NULL
    }
    if (learning_rate$method=="ADADELTA") {
        E_gradient_beta_sqrd <- 0
        E_delta_beta_sqrd<- 0
        E_gradient_alpha_sqrd <- 0
        E_delta_alpha_sqrd <- 0
        decay_rate <- learning_rate$rho
        epsilon <- 10^-6
    }
    if (learning_rate$method=="ADAGRAD") {
        diag_G_mu <- 0
        diag_G_gamma <- 0
        base_learning_rate <- learning_rate$rho
    }
    if (learning_rate$method=="VANILLA") {
        stepsize <- learning_rate$rho
    }
    res<-list()
    counter <- 0
    converged <- FALSE
    values_loglikelihood <- NULL
    
    if (correct_for_zeros && isTRUE(is.element(0,y))) {
        y[which(y==0)] <- 1e-100
    }
    
    if (correct_for_N && isTRUE(is.element(N,y))) {
        y[which(y==N)] <- N-(1e-10)
    }
    
    # sga iterations for beta and alpha
    if (progress) pb = txtProgressBar(min = 0,max = max_iterations,initial = 0,style=3)
    
    for (i in 1:max_iterations) {
        
        if (progress) setTxtProgressBar(pb,i)
        
        counter <- counter+1
        
        beta_0 <- beta
        alpha_0 <- alpha
        
        subsample <- sample(c(1:n),batch_size)
        
        dropout_matrices <- dropout_matrices(noise_distribution,dropout_mu,dropout_gamma,batch_size,k,l)
        
        stoch_gradient <- gradient_def(y[subsample],B_mu[subsample,],B_gamma[subsample,],beta_0,alpha_0,
                                       phi,nu[subsample],b,Db,Db_inv,dropout_matrices$dropout_matrix_mu,
                                       dropout_matrices$dropout_matrix_gamma)
        
        gradient_beta <- stoch_gradient$gradient_beta
        gradient_alpha <- stoch_gradient$gradient_alpha
        
        if (learning_rate$method=="ADADELTA") {
            # beta
            E_gradient_beta_sqrd <- decay_rate*E_gradient_beta_sqrd+(1-decay_rate)*gradient_beta^2
            update_beta <- (sqrt(E_delta_beta_sqrd+epsilon)/sqrt(E_gradient_beta_sqrd+epsilon))*gradient_beta
            E_delta_beta_sqrd <- decay_rate*E_delta_beta_sqrd+(1-decay_rate)*update_beta^2
            beta <- beta_0 + update_beta
            # alpha
            E_gradient_alpha_sqrd <- decay_rate*E_gradient_alpha_sqrd+(1-decay_rate)*gradient_alpha^2
            update_alpha <- (sqrt(E_delta_alpha_sqrd+epsilon)/sqrt(E_gradient_alpha_sqrd+epsilon))*gradient_alpha
            E_delta_alpha_sqrd <- decay_rate*E_delta_alpha_sqrd+(1-decay_rate)*update_alpha^2
            alpha <- alpha_0 + update_alpha
        } else if (learning_rate$method=="ADAGRAD") {
            # beta
            diag_G_mu <- diag_G_mu + gradient_beta^2
            beta <- beta_0 + base_learning_rate*(diag_G_mu^(-1/2))*gradient_beta
            # alpha
            diag_G_gamma <- diag_G_gamma + gradient_alpha^2
            alpha <- alpha_0 + base_learning_rate*(diag_G_gamma^(-1/2))*gradient_alpha
        } else if (learning_rate$method=="VANILLA") {
            # beta
            beta <- beta_0 + stepsize*gradient_beta
            # alpha
            alpha <- alpha_0 + stepsize*gradient_alpha
        }
        
        values_loglikelihood[i] <- loglikelihood_def(y,B_mu,B_gamma,beta,alpha,phi,nu,b,Db_inv)

        # check stopping criterion
        if (i>=100) {
            mean_10 <- mean(values_loglikelihood[0.9*length(values_loglikelihood):length(values_loglikelihood)])
            mean_5 <- mean(values_loglikelihood[0.95*length(values_loglikelihood):length(values_loglikelihood)])
            diff <- abs(mean_10-mean_5)
            if (isTRUE(is.finite(diff) & diff < tolerance)) {
                converged <- TRUE
                break
            }
        }
    }
    
    if (progress) close(pb)
    if (time) print(Sys.time()-start_time)
    
    # fill the empty result list
    res$iterations <- counter
    res$converged <- converged
    res$beta <- as.vector(beta)
    res$alpha <- as.vector(alpha)
    res$meanest <- as.vector(Db(B_mu%*%beta))
    res$dispest <- as.vector(exp(B_gamma%*%alpha))
    res$phi <- phi
    res$varest <- as.vector(phi*D2b(B_mu%*%beta)/exp(B_gamma%*%alpha))
    res$values_loglikelihood <- values_loglikelihood
    if (save_stopping_crit) {
        res$values_stopping_crit <- values_stopping_crit
    }
    res
}

##### PLOTTING ESTIMATION #####

plot_mean_dispersion <- function(estimate,y=FALSE,X,truemean,truedisp,truephi,plot_dispersion=TRUE,
                                 rangemean=range(estimate$meanest,truemean(X),y),
                                 rangedisp=range(estimate$dispest,truedisp(X)),rangevar=range(estimate$varest,truephi/truedisp(X))){
    par(mfrow=c(1,2))
    if (length(y)>1) {
        plot(X,y,pch=20,bty="l",cex=0.4,ylab="",xlab="X",ylim=rangemean)
        lines(X,truemean(X),bty="l",col=6,type="l",lwd=3)
    } else {
        plot(X,truemean(X),bty="l",col=6,type="l",lwd=3,ylab="",xlab="X",ylim=rangemean)
    }
    grid(col="grey",lty = "solid",lwd = 0.5)
    lines(X,estimate$meanest,col=4,type="l",lwd=3)
    legend("topright", legend=c("true mean","estimated mean"),col=c(6,4),lty=1,lwd=3, cex=1, box.lty=0,bg="transparent")
    if (plot_dispersion) {
        plot(X,truedisp(X),bty="l",col=6,type="l",lwd=3,ylab="",xlab="X",ylim=rangedisp)
        grid(col="grey",lty = "solid",lwd = 0.5)
        lines(X,estimate$dispest,col=4,type="l",lwd=3)
        legend("topright", legend=c("true dispersion","estimated dispersion"),col=c(6,4),lty=1,lwd=3, cex=1, box.lty=0,bg="transparent")
    } else {
        plot(X,truephi/truedisp(X),bty="l",col=6,type="l",lwd=3,ylab="",xlab="X",ylim=rangevar)
        grid(col="grey",lty = "solid",lwd = 0.5)
        lines(X,estimate$varest,col=4,type="l",lwd=3)
        legend("topright", legend=c("true variance","estimated variance"),col=c(6,4),lty=1,lwd=3, cex=1, box.lty=0,bg="transparent")
    }
    par(mfrow=c(1,1))
}


##### CROSS-VALIDATION #####

# random search k-fold cross-validation
cross_validation_def_dropout <- function(y,B_mu,B_gamma,phi=1,nu=1,N=NULL,b,Db,Db_inv,D2b,correct_for_zeros,correct_for_N,
                                         noise_distribution,interval_dropout_mu,interval_dropout_gamma,draws=100,
                                         learning_rate=list(method="ADADELTA",rho=0.95),batch_size=30,
                                         max_iterations=50000,tolerance=10^-5,folds_cv=10,progress=FALSE,cores=4) {
    
    # check if batch size is small enough
    if(floor(length(y)*((folds_cv-1)/folds_cv))<batch_size){
        message("The batch size of SGA and the fold size of the k-fold cross-validation are incompatible. 
                Adjust accordingly!")
    }
    
    # preparations
    if(nu==1){
        nu <- rep(1,length(y))
    }
    # correct for zeros in y
    # NOTE: no further correction is necessary when calling sga_def_dropout
    if (correct_for_zeros && isTRUE(is.element(0,y))) {
        y[which(y==0)] <- 1e-100
    }
    if (correct_for_N && isTRUE(is.element(N,y))) {
        y[which(y==N)] <- N-(1e-10)
    }
    folds <- createFolds(y,k=folds_cv)
    grid_dropout_mu <- runif(draws,min(interval_dropout_mu),max(interval_dropout_mu))
    grid_dropout_gamma <- runif(draws,min(interval_dropout_gamma),max(interval_dropout_gamma))
    
    # find optimal hyperparameters
    cv_values <- NULL
    if (progress) pb = txtProgressBar(min = 0,max = draws,initial = 0,style=3) 
    for(j in 1:draws){
        if (progress) setTxtProgressBar(pb,j)
        models <- mclapply(folds, function(x) sga_def_dropout(y=y[-x],B_mu=B_mu[-x,],B_gamma=B_gamma[-x,],beta=0,alpha=0,
                                                              phi=phi,nu=nu[-x],N=N,b=b,Db=Db,Db_inv=Db_inv,D2b=D2b,correct_for_zeros=FALSE,correct_for_N=FALSE,
                                                              noise_distribution=noise_distribution,dropout_mu=grid_dropout_mu[j],dropout_gamma=grid_dropout_gamma[j],
                                                              learning_rate=learning_rate,batch_size=batch_size,
                                                              max_iterations=max_iterations,tolerance=tolerance,progress=FALSE,time=FALSE,
                                                              save_stopping_crit=FALSE),mc.cores=cores)
        likelihoods <- mcmapply(function(x,z) loglikelihood_def(y=y[x],B_mu=B_mu[x,],B_gamma=B_gamma[x,],
                                                                beta=z$beta,alpha=z$alpha,phi=phi,nu=nu[x],
                                                                b=b,Db_inv=Db_inv),folds,models,mc.cores=cores)
        cv_values[j] <- -(1/folds_cv)*sum(likelihoods)
    }
    if (progress) close(pb)
    cv_min <- which.min(cv_values)
    opt_dropout_mu <- grid_dropout_mu[cv_min]
    opt_dropout_gamma <- grid_dropout_gamma[cv_min]
    
    # optimal model
    opt_model <- sga_def_dropout(y=y,B_mu=B_mu,B_gamma=B_gamma,beta=0,alpha=0,
                                 phi=phi,nu=nu,N=N,b=b,Db=Db,Db_inv=Db_inv,D2b=D2b,correct_for_zeros=FALSE,correct_for_N=FALSE,
                                 noise_distribution=noise_distribution,dropout_mu=opt_dropout_mu,dropout_gamma=opt_dropout_gamma,
                                 learning_rate=learning_rate,batch_size=batch_size,
                                 max_iterations=max_iterations,tolerance=tolerance,progress=FALSE,time=FALSE,
                                 save_stopping_crit=FALSE)
    
    # results
    res <- list()
    res$opt_dropout_mu <- opt_dropout_mu
    res$opt_dropout_gamma <- opt_dropout_gamma
    res$opt_beta <- opt_model$beta
    res$opt_alpha <- opt_model$alpha
    res$opt_meanest <- opt_model$meanest
    res$opt_dispest <- opt_model$dispest
    res$opt_varest <- opt_model$varest
    res$grid_dropout_mu <- grid_dropout_mu
    res$grid_dropout_gamma <- grid_dropout_gamma
    res$cv_values <- cv_values
    res
}

aise <- function(true,estimate,grid) {
    sum((estimate(grid)-true(grid))^2)/sum(true(grid))
}

get_rmse <- function(true,estimate,grid) {
    sqrt(mean((estimate(grid)-true(grid))^2))
}


##### SIMULATION #####

simulate_double_normal <- function(x,truemean,truedisp,phi) {
    samplesize <- length(x)
    mean <- truemean(x)
    disp <- truedisp(x)
    y <- mean+rnorm(samplesize,0,sqrt(phi/disp))
    y
}

simulate_double_poisson <- function(x,truemean,truedisp) {
    samplesize <- length(x)
    mean <- truemean(x)
    disp <- truedisp(x)
    y <- NULL
    for (i in 1:samplesize) {
        # rDPO requires mu to be the intensity of the Poisson distribution,
        # such that we need to pass mean[i] as lambda
        lambda <- mean[i]
        gamma <- disp[i]^{-1}
        y[i] <- rDPO(1,lambda,gamma)
    }
    y
}

simulate_double_binomial <- function(x,truemean,truedisp,n) {
    samplesize <- length(x)
    mean <- truemean(x)
    disp <- truedisp(x)
    y <- NULL
    for (i in 1:samplesize) {
        # rDBI requires mu to be the success probability, such that we
        # need to pass mean[i]/n as success probability
        mu <- mean[i]/n
        gamma <- disp[i]^{-1}
        y[i] <- rDBI(1,mu,gamma,n)
    }
    y
}

simulate_data <- function(truemean,truedisp,replicates,samplesizes,distribution=list("GAUSSIAN",phi)) {
    
    for (i in 1:length(samplesizes)) {
        X <- matrix(data=sort(runif(samplesizes[i],0,1)),nrow=samplesizes[i])
        for (j in 2:replicates) {
            X <- cbind(X,sort(runif(samplesizes[i],0,1)))
        }
        Z <- matrix(data=sort(runif(samplesizes[i],0,1)),nrow=samplesizes[i])
        for (j in 2:replicates) {
            Z <- cbind(Z,sort(runif(samplesizes[i],0,1)))
        }
        y <- matrix(data=truemean(X[,1])+rnorm(samplesizes[i],0,sqrt(phi/exp(truedisp(Z[,1])))),nrow=samplesizes[i])
        for (j in 2:replicates) {
            y <- cbind(y,truemean(X[,j])+rnorm(samplesizes[i],0,sqrt(phi/exp(truedisp(Z[,j])))))
        }
    }
    data <- list(y,X,Z)
    data
}

simulation_def <- function (data_distribution,phi=1,nu=1,N=NULL,b,Db,Db_inv,D2b,X,truemean,truemean_index,truedisp,
                            replicates,samplesizes,interval_p_mu,interval_p_gamma,interval_sd_mu,interval_sd_gamma,
                            cv_draws,cv_folds,learning_rate=list(method="ADADELTA",rho=0.95),
                            batch_size,m_mu,m_gamma,cores,progress=FALSE) {
    
    # to be filled with simulation results
    results <<- list()
    correct_for_zeros <- FALSE
    correct_for_N <- FALSE
    if (data_distribution=="poisson") {
        correct_for_zeros <- TRUE
    }
    if (data_distribution=="binomial") {
        correct_for_zeros <- TRUE
        correct_for_N <- TRUE
    }
    
    # equidistant knots
    knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
    knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
    
    for (i in 1:length(samplesizes)) {
        
        if (progress) print(paste("Case ",sep="",i,": samplesize of ",samplesizes[i]))
        
        # simulate replicate datasets for all samplesizes
        if (data_distribution=="normal") {
            y <- matrix(data=simulate_double_normal(X[[i]][,1],truemean,truedisp,phi),nrow=samplesizes[i])
            for (j in 2:replicates) {
                y <- cbind(y,simulate_double_normal(X[[i]][,j],truemean,truedisp,phi))
            }
        } else if (data_distribution=="poisson") {
            y <- matrix(data=simulate_double_poisson(X[[i]][,1],truemean,truedisp),nrow=samplesizes[i])
            for (j in 2:replicates) {
                y <- cbind(y,simulate_double_poisson(X[[i]][,j],truemean,truedisp))
            }
        } else if (data_distribution=="binomial") {
            y <- matrix(data=simulate_double_binomial(X[[i]][,1],truemean,truedisp,n=N),nrow=samplesizes[i])
            for (j in 2:replicates) {
                y <- cbind(y,simulate_double_binomial(X[[i]][,j],truemean,truedisp,n=N))
            }
        }

        # cross-validation for dropout noise
        B_mu <- cbind(rep(1,samplesizes[i]),ns(x=X[[i]][,1],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
        B_gamma <- ns(x=X[[i]][,1],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
        
        if (progress) print("Cross-validation for Bernoulli dropout on first replicate:")
        cv_bernoulli <- cross_validation_def_dropout(y=y[,1],B_mu=B_mu,B_gamma=B_gamma,phi=phi,nu=nu,N=N,b=b,Db=Db,Db_inv=Db_inv,D2b=D2b,
                                                     correct_for_zeros=correct_for_zeros,correct_for_N=correct_for_N,noise_distribution="Bernoulli",
                                                     interval_dropout_mu=interval_p_mu,interval_dropout_gamma=interval_p_gamma,
                                                     draws=cv_draws,learning_rate=learning_rate,batch_size=batch_size,
                                                     max_iterations=50000,tolerance=10^-5,folds_cv=cv_folds,
                                                     progress=progress,cores=cores)
        
        if (progress) print("Cross-validation for Normal dropout on first replicate:")
        cv_normal <- cross_validation_def_dropout(y=y[,1],B_mu=B_mu,B_gamma=B_gamma,phi=phi,nu=nu,N=N,b=b,Db=Db,Db_inv=Db_inv,D2b=D2b,
                                                  correct_for_zeros=correct_for_zeros,correct_for_N=correct_for_N,noise_distribution="Normal",
                                                  interval_dropout_mu=interval_sd_mu,interval_dropout_gamma=interval_sd_gamma,
                                                  draws=cv_draws,learning_rate=learning_rate,batch_size=batch_size,
                                                  max_iterations=50000,tolerance=10^-5,folds_cv=cv_folds,
                                                  progress=progress,cores=cores)
        # cross-validation for PMLE
        if (progress) print("Cross-validation for PMLE on first replicate...")
        seq_mu <- c(0,5000)
        seq_gamma <- c(0,5000)
        if (data_distribution=="normal") {
            cv_pmle <- NormalBothGCV(y=y[,1],Bmu=B_mu,Bgamma=B_gamma,seqM=seq_mu,seqG=seq_gamma,
                                          lenM=50,lenG=50,minlM=10^-5,maxlM=15000,minlG=10^-5,maxlG=15000,
                                          plotM=FALSE,plotG=FALSE)
        } else if (data_distribution=="poisson") {
            cv_pmle <- poisson2sGCV(y=y[,1],Bmu=B_mu,Bgamma=B_gamma,seqM=seq_mu,seqG=seq_gamma,
                                         lenM=50,lenG=50,minlM=10^-5,maxlM=15000,minlG=10^-5,maxlG=15000,
                                         plotM=FALSE,plotG=FALSE)
        } else if (data_distribution=="binomial") {
            prob <- y[,1]/N
            cv_pmle <- BinomialBothGCV(y=prob,num=N,Bmu=B_mu,Bgamma=B_gamma,seqM=seq_mu,seqG=seq_gamma,
                                            lenM=50,lenG=50,minlM=10^-5,maxlM=15000,minlG=10^-5,maxlG=15000,
                                            plotM=FALSE,plotG=FALSE)
        }
        
        # get estimates
        if (progress) print("Estimation with Bernoulli dropout for all replicates...")
        estimates_bernoulli <- foreach (k=2:replicates) %dopar% {
            B_mu <- cbind(rep(1,samplesizes[i]),ns(x=X[[i]][,k],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
            B_gamma <- ns(x=X[[i]][,k],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
            sga_def_dropout(y=y[,k],B_mu=B_mu,B_gamma=B_gamma,beta=0,alpha=0,phi=phi,nu=nu,N=N,b=b,Db=Db,Db_inv=Db_inv,D2b=D2b,correct_for_zeros=correct_for_zeros,correct_for_N=correct_for_N,
                                     noise_distribution="Bernoulli",dropout_mu=cv_bernoulli$opt_dropout_mu,dropout_gamma=cv_bernoulli$opt_dropout_gamma,
                                     learning_rate=learning_rate,batch_size=batch_size,max_iterations=50000,tolerance=10^-5,progress=FALSE,time=FALSE,
                                     save_stopping_crit=FALSE)
        }
        
        if (progress) print("Estimation with Normal dropout for all replicates...")
        estimates_normal <- foreach (k=2:replicates) %dopar% {
            B_mu <- cbind(rep(1,samplesizes[i]),ns(x=X[[i]][,k],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
            B_gamma <- ns(x=X[[i]][,k],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
            sga_def_dropout(y=y[,k],B_mu=B_mu,B_gamma=B_gamma,beta=0,alpha=0,phi=phi,nu=nu,N=N,b=b,Db=Db,Db_inv=Db_inv,D2b=D2b,correct_for_zeros=correct_for_zeros,correct_for_N=correct_for_N,
                            noise_distribution="Normal",dropout_mu=cv_normal$opt_dropout_mu,dropout_gamma=cv_normal$opt_dropout_gamma,
                            learning_rate=learning_rate,batch_size=batch_size,max_iterations=50000,tolerance=10^-5,progress=FALSE,time=FALSE,
                            save_stopping_crit=FALSE)
        }
        
        if (progress) print("Estimation with PMLE for all replicates...")
        estimates_pmle <- foreach (k=2:replicates) %dopar% {
            B_mu <- cbind(rep(1,samplesizes[i]),ns(x=X[[i]][,k],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
            B_gamma <- ns(x=X[[i]][,k],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
            if (data_distribution=="normal") {
                fisher_def_normal(y=y[,k],Bmu=B_mu,Bgamma=B_gamma,alphaM=0,alphaG=0,phi=phi,lambdaM=cv_pmle$vals[1],
                                  lambdaG=cv_pmle$vals[2],maxitG=35,maxitT=2000,penordM=2,penordG=2,tol=10^-8)
            } else if (data_distribution=="poisson") {
                fisher_def_poisson(y=y[,k],Bmu=B_mu,Bgamma=B_gamma,alphaM=0,alphaG=0,lambdaM=cv_pmle$vals[1],
                                   lambdaG=cv_pmle$vals[2],maxitG=100,maxitT=2000,penordM=2,penordG=2,tol=10^-8)
            } else if (data_distribution=="binomial") {
                prob <- y[,k]/N
                fisher_def_binomial(y=prob,num=N,Bmu=B_mu,Bgamma=B_gamma,lambdaM=cv_pmle$vals[1],
                                    lambdaG=cv_pmle$vals[2],penordM=2,penordG=2,tol=10^-10,type="quasi",
                                    alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp",maxitT=2000,maxitM=60,maxitG=45)
            }
        }
        
        if (progress) print("Computing RMSE for Bernoulli, Normal and PMLE estimates...")
        grid <- seq(from=0,to=1,length.out=501)
        
        cv_bernoulli$rmse_mean <- get_rmse(truemean,get_est(x,knots_mu,cv_bernoulli$opt_beta,Db,TRUE),grid)
        rmse_mean_bernoulli <- foreach(k=1:(replicates-1), .combine=c) %dopar% {
            mean_est_bernoulli <- get_est(x,knots_mu,estimates_bernoulli[[k]]$beta,Db,TRUE)
            get_rmse(truemean,mean_est_bernoulli,grid)
        }
        
        cv_bernoulli$rmse_disp <- get_rmse(truedisp,get_est(x,knots_gamma,cv_bernoulli$opt_alpha,function(x){exp(x)}),grid)
        rmse_disp_bernoulli <- foreach(k=1:(replicates-1), .combine=c) %dopar% {
            disp_est_bernoulli <- get_est(x,knots_gamma,estimates_bernoulli[[k]]$alpha,function(x){exp(x)})
            get_rmse(truedisp,disp_est_bernoulli,grid)
        }
        
        cv_normal$rmse_mean <- get_rmse(truemean,get_est(x,knots_mu,cv_normal$opt_beta,Db,TRUE),grid)
        rmse_mean_normal <- foreach(k=1:(replicates-1), .combine=c) %dopar% {
            mean_est_normal <- get_est(x,knots_mu,estimates_normal[[k]]$beta,Db,TRUE)
            get_rmse(truemean,mean_est_normal,grid)
        }
        
        cv_normal$rmse_disp <- get_rmse(truedisp,get_est(x,knots_gamma,cv_normal$opt_alpha,function(x){exp(x)}),grid)
        rmse_disp_normal <- foreach(k=1:(replicates-1), .combine=c) %dopar% {
            disp_est_normal <- get_est(x,knots_gamma,estimates_normal[[k]]$alpha,function(x){exp(x)})
            get_rmse(truedisp,disp_est_normal,grid)
        }
        
        cv_pmle$rmse_mean <- get_rmse(truemean,get_est(x,knots_mu,cv_pmle$alphaM,Db,TRUE),grid)
        rmse_mean_pmle <- foreach(k=1:(replicates-1), .combine=c) %dopar% {
            mean_est_pmle <- get_est(x,knots_mu,estimates_pmle[[k]]$alphaM,Db,TRUE)
            get_rmse(truemean,mean_est_pmle,grid)
        }
        
        cv_pmle$rmse_disp <- get_rmse(truedisp,get_est(x,knots_gamma,cv_pmle$alphaG,function(x){phi*exp(x)^{-1}}),grid)
        rmse_disp_pmle <- foreach(k=1:(replicates-1), .combine=c) %dopar% {
            disp_est_pmle <- get_est(x,knots_gamma,estimates_pmle[[k]]$alphaG,function(x){phi*exp(x)^{-1}})
            get_rmse(truedisp,disp_est_pmle,grid)
        }
        
        # save results
        results[[i]] <<- list()
        results[[i]]$replicates <<- replicates
        results[[i]]$samplesize <<- samplesizes[i]
        results[[i]]$data <<- list(y,X[[i]])
        results[[i]]$cv_bernoulli <<- cv_bernoulli
        results[[i]]$cv_normal <<- cv_normal
        results[[i]]$cv_pmle <<- cv_pmle
        results[[i]]$estimates_bernoulli <<- estimates_bernoulli
        results[[i]]$estimates_normal <<- estimates_normal
        results[[i]]$estimates_pmle <<- estimates_pmle
        results[[i]]$rmse_mean_bernoulli <<- rmse_mean_bernoulli
        results[[i]]$rmse_disp_bernoulli <<- rmse_disp_bernoulli
        results[[i]]$rmse_mean_normal <<- rmse_mean_normal
        results[[i]]$rmse_disp_normal <<- rmse_disp_normal
        results[[i]]$rmse_mean_pmle <<- rmse_mean_pmle
        results[[i]]$rmse_disp_pmle <<- rmse_disp_pmle
        
        # save intermediate results into environment
        intermediate_image <- paste("intermediate_models_",sep="",truemean_index,"_",data_distribution,".RData")
        save.image(file=intermediate_image)
    }
    results
}

# returns the estimated function for mu or lambda
get_est <- function(x,knots,parameter,link_inv,intercept=FALSE) {
    if (intercept) {
        function(x){link_inv(cbind(rep(1,length(x)),ns(x=x,knots=knots,intercept = FALSE,Boundary.knots=c(0,1)))%*%parameter)}
    } else {
        function(x){link_inv(ns(x=x,knots=knots,intercept = FALSE,Boundary.knots=c(0,1))%*%parameter)}
    }
}

plot_results <- function () {
    break
} 

# location should be given as a path, e.g. "/home/.../filename.png"

plot_simulation <- function (simulation,location,data_distribution,truemean,mean_index,Db,truedisp,disp_index,
                             m_mu,m_gamma,phi,resolution_scatter,resolution_functions) {
    
    # iterate over different samplesizes
    for (i in 1:length(simulation)) {
        
        # current subsimulation for certain samplesize
        sim <- simulation[[i]]
        replicates <- sim$replicates-1
        
        # scatter plot bernoulli
        scatter_bernoulli <- data.frame(p_mu=sim$cv_bernoulli$grid_dropout_mu,p_gamma=sim$cv_bernoulli$grid_dropout_gamma,cv=sim$cv_bernoulli$cv_values)
        optimum <- data.frame(p_mu=sim$cv_bernoulli$opt_dropout_mu,p_gamma=sim$cv_bernoulli$opt_dropout_gamma)
        title <- paste(location,sep="",data_distribution,"_scatter_",mean_index,disp_index,"_bernoulli_",sim$samplesize,".png")
        png(title,width=500,height=500,res=resolution_scatter,pointsize=10)
        scatter <- ggplot(scatter_bernoulli, aes(x = p_mu, y = p_gamma, color = cv)) + geom_point(size=3,alpha=0.8) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=p_mu,y=p_gamma), color="black", size=4) + ggtitle(paste("n = ",sep="",sim$samplesize))
        print(scatter)
        dev.off()
        
        # scatter plot normal
        scatter_normal <- data.frame(sd_mu=sim$cv_normal$grid_dropout_mu,sd_gamma=sim$cv_normal$grid_dropout_gamma,cv=sim$cv_normal$cv_values)
        optimum <- data.frame(sd_mu=sim$cv_normal$opt_dropout_mu,sd_gamma=sim$cv_normal$opt_dropout_gamma)
        title <- paste(location,sep="",data_distribution,"_scatter_",mean_index,disp_index,"_normal_",sim$samplesize,".png")
        png(title,width=500,height=500,res=resolution_scatter,pointsize=10)
        scatter <- ggplot(scatter_normal, aes(x = sd_mu, y = sd_gamma, color = cv)) + geom_point(size=3,alpha=0.8) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=sd_mu,y=sd_gamma), color="black", size=4) + ggtitle(paste("n = ",sep="",sim$samplesize))
        print(scatter)
        dev.off()
        
        # joint plots
        X <- seq(from=0,to=1,length.out=sim$samplesize)
        knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
        knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
        B_mu <- cbind(rep(1,sim$samplesize),ns(x=X,knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
        B_gamma <- ns(x=X,knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
        
        # results for bernoulli dropout
        range_mean <- range(truemean(X))
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_bernoulli[[k]]$meanest))
            range_disp <- range(range_disp,range(sim$estimates_bernoulli[[k]]$dispest))
        }
        mean_of_means <- 0
        for (k in 1:replicates) {
            mean_of_means <- mean_of_means + (1/replicates)*Db(B_mu%*%sim$estimates_bernoulli[[k]]$beta)
        }
        mean_of_disps <- 0
        for (k in 1:replicates) {
            mean_of_disps <- mean_of_disps + (1/replicates)*exp(B_gamma%*%sim$estimates_bernoulli[[k]]$alpha)
        }
        mean_5percentile <- sim$estimates_bernoulli[order(sim$rmse_mean_bernoulli)][[ceiling(0.05*replicates)]]$meanest
        mean_50percentile <- sim$estimates_bernoulli[order(sim$rmse_mean_bernoulli)][[ceiling(0.5*replicates)]]$meanest
        mean_95percentile <- sim$estimates_bernoulli[order(sim$rmse_mean_bernoulli)][[ceiling(0.95*replicates)]]$meanest
        disp_5percentile <- sim$estimates_bernoulli[order(sim$rmse_disp_bernoulli)][[ceiling(0.05*replicates)]]$dispest
        disp_50percentile <- sim$estimates_bernoulli[order(sim$rmse_disp_bernoulli)][[ceiling(0.5*replicates)]]$dispest
        disp_95percentile <- sim$estimates_bernoulli[order(sim$rmse_disp_bernoulli)][[ceiling(0.95*replicates)]]$dispest
        
        title <- paste(location,sep="",data_distribution,"_",mean_index,disp_index,"_bernoulli_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truemean(X),mean_of_means,mean_5percentile,mean_50percentile,mean_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_means,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_bernoulli)][,ceiling(0.05*replicates)],mean_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_bernoulli)][,ceiling(0.5*replicates)],mean_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_bernoulli)][,ceiling(0.95*replicates)],mean_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true mean","mean of means","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truedisp(X),mean_of_disps,disp_5percentile,disp_50percentile,disp_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_disps,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_bernoulli)][,ceiling(0.05*replicates)],disp_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_bernoulli)][,ceiling(0.5*replicates)],disp_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_bernoulli)][,ceiling(0.95*replicates)],disp_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true dispersion","mean of disps","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        par(mfrow=c(1,1))
        dev.off()
        
        # results for normal dropout
        range_mean <- range(truemean(X))
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_normal[[k]]$meanest))
            range_disp <- range(range_disp,range(sim$estimates_normal[[k]]$dispest))
        }
        mean_of_means <- 0
        for (k in 1:replicates) {
            mean_of_means <- mean_of_means + (1/replicates)*Db(B_mu%*%sim$estimates_normal[[k]]$beta)
        }
        mean_of_disps <- 0
        for (k in 1:replicates) {
            mean_of_disps <- mean_of_disps + (1/replicates)*exp(B_gamma%*%sim$estimates_normal[[k]]$alpha)
        }
        mean_5percentile <- sim$estimates_normal[order(sim$rmse_mean_normal)][[ceiling(0.05*replicates)]]$meanest
        mean_50percentile <- sim$estimates_normal[order(sim$rmse_mean_normal)][[ceiling(0.5*replicates)]]$meanest
        mean_95percentile <- sim$estimates_normal[order(sim$rmse_mean_normal)][[ceiling(0.95*replicates)]]$meanest
        disp_5percentile <- sim$estimates_normal[order(sim$rmse_disp_normal)][[ceiling(0.05*replicates)]]$dispest
        disp_50percentile <- sim$estimates_normal[order(sim$rmse_disp_normal)][[ceiling(0.5*replicates)]]$dispest
        disp_95percentile <- sim$estimates_normal[order(sim$rmse_disp_normal)][[ceiling(0.95*replicates)]]$dispest
        
        title <- paste(location,sep="",data_distribution,"_",mean_index,disp_index,"_normal_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_normal[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_normal[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truemean(X),mean_of_means,mean_5percentile,mean_50percentile,mean_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_means,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_normal)][,ceiling(0.05*replicates)],mean_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_normal)][,ceiling(0.5*replicates)],mean_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_normal)][,ceiling(0.95*replicates)],mean_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true mean","mean of means","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truedisp(X),mean_of_disps,disp_5percentile,disp_50percentile,disp_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_disps,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_normal)][,ceiling(0.05*replicates)],disp_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_normal)][,ceiling(0.5*replicates)],disp_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_normal)][,ceiling(0.95*replicates)],disp_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true dispersion","mean of disps","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        par(mfrow=c(1,1))
        dev.off()
        
        # results for pmle
        range_mean <- range(truemean(X))
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_pmle[[k]]$meanest))
            range_disp <- range(range_disp,range(sim$estimates_pmle[[k]]$dispest))
        }
        mean_of_means <- 0
        for (k in 1:replicates) {
            mean_of_means <- mean_of_means + (1/replicates)*Db(B_mu%*%sim$estimates_pmle[[k]]$alphaM)
        }
        mean_of_disps <- 0
        for (k in 1:replicates) {
            mean_of_disps <- mean_of_disps + (1/replicates)*phi*exp(B_gamma%*%sim$estimates_pmle[[k]]$alphaG)^-1
        }
        mean_5percentile <- sim$estimates_pmle[order(sim$rmse_mean_pmle)][[ceiling(0.05*replicates)]]$meanest
        mean_50percentile <- sim$estimates_pmle[order(sim$rmse_mean_pmle)][[ceiling(0.5*replicates)]]$meanest
        mean_95percentile <- sim$estimates_pmle[order(sim$rmse_mean_pmle)][[ceiling(0.95*replicates)]]$meanest
        disp_5percentile <- sim$estimates_pmle[order(sim$rmse_disp_pmle)][[ceiling(0.05*replicates)]]$dispest
        disp_50percentile <- sim$estimates_pmle[order(sim$rmse_disp_pmle)][[ceiling(0.5*replicates)]]$dispest
        disp_95percentile <- sim$estimates_pmle[order(sim$rmse_disp_pmle)][[ceiling(0.95*replicates)]]$dispest
        
        title <- paste(location,sep="",data_distribution,"_",mean_index,disp_index,"_pmle_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)^
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_pmle[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_pmle[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truemean(X),mean_of_means,mean_5percentile,mean_50percentile,mean_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_means,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_pmle)][,ceiling(0.05*replicates)],mean_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_pmle)][,ceiling(0.5*replicates)],mean_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_mean_pmle)][,ceiling(0.95*replicates)],mean_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true mean","mean of means","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truedisp(X),mean_of_disps,disp_5percentile,disp_50percentile,disp_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_disps,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_pmle)][,ceiling(0.05*replicates)],disp_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_pmle)][,ceiling(0.5*replicates)],disp_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,-1][,order(sim$rmse_disp_pmle)][,ceiling(0.95*replicates)],disp_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true dispersion","mean of disps","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        par(mfrow=c(1,1))
        dev.off()
    }
    
    # boxplots
    RMSE_mean <- NULL
    RMSE_disp <- NULL
    samplesizes <- NULL
    for (i in 1:length(simulation)) {
        RMSE_mean <- c(RMSE_mean,simulation[[i]]$rmse_mean_bernoulli,simulation[[i]]$rmse_mean_normal,simulation[[i]]$rmse_mean_pmle)
        RMSE_disp <- c(RMSE_disp,simulation[[i]]$rmse_disp_bernoulli,simulation[[i]]$rmse_disp_normal,simulation[[i]]$rmse_disp_pmle)
        samplesizes <- c(samplesizes,paste(simulation[[i]]$samplesize))
    }
    
    boxplot_data <- data.frame(RMSE_mean=RMSE_mean,
                               RMSE_disp=RMSE_disp,
                               samplesize=rep(samplesizes,each=replicates*3),
                               Method=rep(rep(c("Bernoulli dropout","Normal dropout","PMLE"),each=replicates),length(simulation)))
    
    title <- paste(location,sep="",data_distribution,"_boxplot_mean_",mean_index,disp_index,".png")
    png(title,width=700,height=400,res=120,pointsize=10)
    boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("samplesize") + ylab("RMSE")
    print(boxplot_mean)
    dev.off()
    
    title <- paste(location,sep="",data_distribution,"_boxplot_disp_",mean_index,disp_index,".png")
    png(title,width=900,height=400,res=120,pointsize=10)
    boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("")
    print(boxplot_disp)
    dev.off()
}

plot_simulation_comparison <- function (list_of_simulations,location,data_distribution,truemean,mean_index,Db,truedisps,disp_indices,
                                        m_mu,m_gamma,phi,resolution_scatter,resolution_functions) {
    n_simulations <- length(list_of_simulations)
    n_samplesizes <- length(list_of_simulations[[1]])
    replicates <- list_of_simulations[[1]][[1]]$replicates-1
    for (i in 1:n_samplesizes) {
        
        samplesize <- list_of_simulations[[1]][[i]]$samplesize
        X <- seq(from=0,to=1,length.out=samplesize)
        knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
        knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
        B_mu <- cbind(rep(1,samplesize),ns(x=X,knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
        B_gamma <- ns(x=X,knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
        
        title <- paste(location,sep="","comparison_",data_distribution,"_",mean_index,"_",samplesize,".png")
        png(title,width=1600,height=800,res=resolution_functions,pointsize=10)
        par(mfrow=c(2,length(list_of_simulations)))
        for (j in 1:length(n_simulations)) {
            range_mean <- range(truemean(X))
            for (k in 1:replicates) {
                range_mean <- range(range_mean,range(sim$estimates_bernoulli[[k]]$meanest))
            }
            
        }
        for (j in 1:length(truedisps)) {
            range_disp <- range(truedisps[[j]](X))
            for (k in 1:replicates) {
                range_disp <- range(range_disp,range(sim$estimates_bernoulli[[k]]$dispest))
            }
        }

        
        title <- paste(location,sep="","comparison_",data_distribution,"_",mean_index,"_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        
    }
    
}

plot_simulation_old <- function (simulation,location,data_distribution,truemean,mean_index,Db,truedisp,disp_index,
                             m_mu,m_gamma,phi,resolution_scatter,resolution_functions) {
    
    # iterate over different samplesizes
    for (i in 1:length(simulation)) {
        
        # current subsimulation for certain samplesize
        sim <- simulation[[i]]
        replicates <- sim$replicates
        
        # scatter plot bernoulli
        scatter_bernoulli <- data.frame(p_mu=sim$cv_bernoulli$grid_sd_mu,p_gamma=sim$cv_bernoulli$grid_sd_gamma,cv=sim$cv_bernoulli$cv_values)
        optimum <- data.frame(p_mu=sim$cv_bernoulli$opt_sd_mu,p_gamma=sim$cv_bernoulli$opt_sd_gamma)
        title <- paste(location,sep="",data_distribution,"_scatter_",mean_index,disp_index,"_bernoulli_",sim$samplesize,".png")
        png(title,width=500,height=500,res=resolution_scatter,pointsize=10)
        scatter <- ggplot(scatter_bernoulli, aes(x = p_mu, y = p_gamma, color = cv)) + geom_point(size=3,alpha=0.8) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=p_mu,y=p_gamma), color="black", size=4) + ggtitle(paste("n = ",sep="",sim$samplesize))
        print(scatter)
        dev.off()
        
        # scatter plot normal
        scatter_normal <- data.frame(sd_mu=sim$cv_normal$grid_sd_mu,sd_gamma=sim$cv_normal$grid_sd_gamma,cv=sim$cv_normal$cv_values)
        optimum <- data.frame(sd_mu=sim$cv_normal$opt_sd_mu,sd_gamma=sim$cv_normal$opt_sd_gamma)
        title <- paste(location,sep="",data_distribution,"_scatter_",mean_index,disp_index,"_normal_",sim$samplesize,".png")
        png(title,width=500,height=500,res=resolution_scatter,pointsize=10)
        scatter <- ggplot(scatter_normal, aes(x = sd_mu, y = sd_gamma, color = cv)) + geom_point(size=3,alpha=0.8) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=sd_mu,y=sd_gamma), color="black", size=4) + ggtitle(paste("n = ",sep="",sim$samplesize))
        print(scatter)
        dev.off()
        
        # joint plots
        X <- seq(from=0,to=1,length.out=sim$samplesize)
        knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
        knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
        B_mu <- cbind(rep(1,sim$samplesize),ns(x=X,knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
        B_gamma <- ns(x=X,knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
        
        # results for bernoulli dropout
        range_mean <- range(truemean(X))
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_bernoulli[[k]]$meanest))
            range_disp <- range(range_disp,range(sim$estimates_bernoulli[[k]]$dispest))
        }
        mean_of_means <- 0
        for (k in 1:replicates) {
            mean_of_means <- mean_of_means + (1/replicates)*Db(B_mu%*%sim$estimates_bernoulli[[k]]$beta)
        }
        mean_of_disps <- 0
        for (k in 1:replicates) {
            mean_of_disps <- mean_of_disps + (1/replicates)*exp(B_gamma%*%sim$estimates_bernoulli[[k]]$alpha)
        }
        mean_5percentile <- sim$estimates_bernoulli[order(sim$rmse_mean_est_bernoulli)][[ceiling(0.05*replicates)]]$meanest
        mean_50percentile <- sim$estimates_bernoulli[order(sim$rmse_mean_est_bernoulli)][[ceiling(0.5*replicates)]]$meanest
        mean_95percentile <- sim$estimates_bernoulli[order(sim$rmse_mean_est_bernoulli)][[ceiling(0.95*replicates)]]$meanest
        disp_5percentile <- sim$estimates_bernoulli[order(sim$rmse_disp_est_bernoulli)][[ceiling(0.05*replicates)]]$dispest
        disp_50percentile <- sim$estimates_bernoulli[order(sim$rmse_disp_est_bernoulli)][[ceiling(0.5*replicates)]]$dispest
        disp_95percentile <- sim$estimates_bernoulli[order(sim$rmse_disp_est_bernoulli)][[ceiling(0.95*replicates)]]$dispest
        
        title <- paste(location,sep="",data_distribution,"_",mean_index,disp_index,"_bernoulli_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k],sim$estimates_bernoulli[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k],sim$estimates_bernoulli[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truemean(X),mean_of_means,mean_5percentile,mean_50percentile,mean_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_means,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_est_bernoulli)][,ceiling(0.05*replicates)],mean_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_est_bernoulli)][,ceiling(0.5*replicates)],mean_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_est_bernoulli)][,ceiling(0.95*replicates)],mean_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true mean","mean of means","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truedisp(X),mean_of_disps,disp_5percentile,disp_50percentile,disp_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_disps,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_est_bernoulli)][,ceiling(0.05*replicates)],disp_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_est_bernoulli)][,ceiling(0.5*replicates)],disp_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_est_bernoulli)][,ceiling(0.95*replicates)],disp_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true dispersion","mean of disps","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        par(mfrow=c(1,1))
        dev.off()
        
        # results for normal dropout
        range_mean <- range(truemean(X))
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_normal[[k]]$meanest))
            range_disp <- range(range_disp,range(sim$estimates_normal[[k]]$dispest))
        }
        mean_of_means <- 0
        for (k in 1:replicates) {
            mean_of_means <- mean_of_means + (1/replicates)*Db(B_mu%*%sim$estimates_normal[[k]]$beta)
        }
        mean_of_disps <- 0
        for (k in 1:replicates) {
            mean_of_disps <- mean_of_disps + (1/replicates)*exp(B_gamma%*%sim$estimates_normal[[k]]$alpha)
        }
        mean_5percentile <- sim$estimates_normal[order(sim$rmse_mean_est_normal)][[ceiling(0.05*replicates)]]$meanest
        mean_50percentile <- sim$estimates_normal[order(sim$rmse_mean_est_normal)][[ceiling(0.5*replicates)]]$meanest
        mean_95percentile <- sim$estimates_normal[order(sim$rmse_mean_est_normal)][[ceiling(0.95*replicates)]]$meanest
        disp_5percentile <- sim$estimates_normal[order(sim$rmse_disp_est_normal)][[ceiling(0.05*replicates)]]$dispest
        disp_50percentile <- sim$estimates_normal[order(sim$rmse_disp_est_normal)][[ceiling(0.5*replicates)]]$dispest
        disp_95percentile <- sim$estimates_normal[order(sim$rmse_disp_est_normal)][[ceiling(0.95*replicates)]]$dispest
        
        title <- paste(location,sep="",data_distribution,"_",mean_index,disp_index,"_normal_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k],sim$estimates_normal[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k],sim$estimates_normal[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truemean(X),mean_of_means,mean_5percentile,mean_50percentile,mean_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_means,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_est_normal)][,ceiling(0.05*replicates)],mean_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_est_normal)][,ceiling(0.5*replicates)],mean_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_est_normal)][,ceiling(0.95*replicates)],mean_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true mean","mean of means","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truedisp(X),mean_of_disps,disp_5percentile,disp_50percentile,disp_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_disps,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_est_normal)][,ceiling(0.05*replicates)],disp_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_est_normal)][,ceiling(0.5*replicates)],disp_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_est_normal)][,ceiling(0.95*replicates)],disp_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true dispersion","mean of disps","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        par(mfrow=c(1,1))
        dev.off()
        
        # results for pmle
        range_mean <- range(truemean(X))
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$benchmarks[[k]]$meanest))
            range_disp <- range(range_disp,range(sim$benchmarks[[k]]$dispest))
        }
        mean_of_means <- 0
        for (k in 1:replicates) {
            mean_of_means <- mean_of_means + (1/replicates)*Db(B_mu%*%sim$benchmarks[[k]]$alphaM)
        }
        mean_of_disps <- 0
        for (k in 1:replicates) {
            mean_of_disps <- mean_of_disps + (1/replicates)*phi*exp(B_gamma%*%sim$benchmarks[[k]]$alphaG)^-1
        }
        mean_5percentile <- sim$benchmarks[order(sim$rmse_mean_bench)][[ceiling(0.05*replicates)]]$meanest
        mean_50percentile <- sim$benchmarks[order(sim$rmse_mean_bench)][[ceiling(0.5*replicates)]]$meanest
        mean_95percentile <- sim$benchmarks[order(sim$rmse_mean_bench)][[ceiling(0.95*replicates)]]$meanest
        disp_5percentile <- sim$benchmarks[order(sim$rmse_disp_bench)][[ceiling(0.05*replicates)]]$dispest
        disp_50percentile <- sim$benchmarks[order(sim$rmse_disp_bench)][[ceiling(0.5*replicates)]]$dispest
        disp_95percentile <- sim$benchmarks[order(sim$rmse_disp_bench)][[ceiling(0.95*replicates)]]$dispest
        
        title <- paste(location,sep="",data_distribution,"_",mean_index,disp_index,"_pmle_",sim$samplesize,".png")
        png(title,width=1600,height=500,res=resolution_functions,pointsize=10)
        par(mfrow=c(1,4))
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_mean)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k],sim$benchmarks[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range_disp)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
            lines(sim$data[[2]][,k],sim$benchmarks[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
        }
        legend("topright", legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truemean(X),mean_of_means,mean_5percentile,mean_50percentile,mean_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_means,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_bench)][,ceiling(0.05*replicates)],mean_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_bench)][,ceiling(0.5*replicates)],mean_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_mean_bench)][,ceiling(0.95*replicates)],mean_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true mean","mean of means","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="X",ylim=range(truedisp(X),mean_of_disps,disp_5percentile,disp_50percentile,disp_95percentile))
        grid(col="grey",lty = "solid",lwd = 0.3)
        lines(X,mean_of_disps,col="darkgreen",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_bench)][,ceiling(0.05*replicates)],disp_5percentile,col="red",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_bench)][,ceiling(0.5*replicates)],disp_50percentile,col="blue",type="l",lwd=1)
        lines(sim$data[[2]][,order(sim$rmse_disp_bench)][,ceiling(0.95*replicates)],disp_95percentile,col="gold",type="l",lwd=1)
        legend("topright", legend=c("true dispersion","mean of disps","5th RMSE percentile","50th RMSE percentile","95th RMSE percentile"),col=c("black","darkgreen","red","blue","gold"),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        par(mfrow=c(1,1))
        dev.off()
    }
    
    # boxplots
    
    boxplot_data <- data.frame(RMSE_mean=c(simulation[[1]]$rmse_mean_est_bernoulli,simulation[[1]]$rmse_mean_est_normal,simulation[[1]]$rmse_mean_bench,
                                           simulation[[2]]$rmse_mean_est_bernoulli,simulation[[2]]$rmse_mean_est_normal,simulation[[2]]$rmse_mean_bench,
                                           simulation[[3]]$rmse_mean_est_bernoulli,simulation[[3]]$rmse_mean_est_normal,simulation[[3]]$rmse_mean_bench),
                               RMSE_disp=c(simulation[[1]]$rmse_disp_est_bernoulli,simulation[[1]]$rmse_disp_est_normal,simulation[[1]]$rmse_disp_bench,
                                           simulation[[2]]$rmse_disp_est_bernoulli,simulation[[2]]$rmse_disp_est_normal,simulation[[2]]$rmse_disp_bench,
                                           simulation[[3]]$rmse_disp_est_bernoulli,simulation[[3]]$rmse_disp_est_normal,simulation[[3]]$rmse_disp_bench),
                               samplesize=rep(c("250","500","1000"),each=300),
                               Method=rep(rep(c("Bernoulli dropout","Normal dropout","PMLE"),each=100),3))
    
    title <- paste(location,sep="",data_distribution,"_boxplot_mean_",mean_index,disp_index,".png")
    png(title,width=700,height=400,res=120,pointsize=10)
    boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=c("250","500","1000")), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("samplesize") + ylab("RMSE")
    print(boxplot_mean)
    dev.off()
    
    title <- paste(location,sep="",data_distribution,"_boxplot_disp_",mean_index,disp_index,".png")
    png(title,width=900,height=400,res=120,pointsize=10)
    boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=c("250","500","1000")), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("")
    print(boxplot_disp)
    dev.off()
}

correct_simulation <- function (model_bernoulli,model_normal,data_distribution,N=NULL,b,Db,Db_inv,phi,truemean,truedisp,m_mu,m_gamma) {
    # corrections:
    # 0: 
    # 1: redo cv_benchmark on first sample
    # 2: recompute benchmarks
    # 3: add cv_normal and estimates_normal
    # 4: recompute rmses for all functions
    
    start_time <- Sys.time()
    
    # model to correct and join
    model <- model_bernoulli
    model_normal <- model_normal
    data_distribution <- data_distribution
    Db_current <- Db
    phi <- phi
    truemean <- truemean
    truedisp <- truedisp
    
    knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
    knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
    samplesizes <- NULL
    for (i in 1:length(model)) {
        samplesizes[i] <- model[[i]]$samplesize
    }
    replicates <- model[[1]]$replicates
    
    cores <- detectCores()
    registerDoParallel(cores)
    
    for(i in 1:length(samplesizes)) {
        
        # renaming
        names(model[[i]])[names(model[[i]])=="cv_results"] <- "cv_bernoulli"
        names(model[[i]])[names(model[[i]])=="estimates"] <- "estimates_bernoulli"
        
        # deleting
        model[[i]]$cv_benchmark <- NULL
        model[[i]]$benchmarks <- NULL
        model[[i]]$rmse_mean_est <- NULL
        model[[i]]$rmse_disp_est <- NULL
        model[[i]]$rmse_mean_bench <- NULL
        model[[i]]$rmse_disp_bench <- NULL
        
        # joining
        model[[i]] <- append(model[[i]],list(cv_normal=model_normal[[i]]$cv_results),which(names(model[[i]])=="cv_bernoulli"))
        model[[i]] <- append(model[[i]],list(estimates_normal=model_normal[[i]]$estimates),which(names(model[[i]])=="estimates_bernoulli"))
        
        # correct estimates where zeros exist
        print("Recompute estimates whith zeros in dependent variable...")
        if (is.element(0,model[[i]]$data[[1]])) {
            for (k in 1:replicates) {
                if (is.element(0,model[[i]]$data[[1]][,k])) {
                    B_mu <- cbind(rep(1,samplesizes[i]),ns(x=model[[i]]$data[[2]][,k],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
                    B_gamma <- ns(x=model[[i]]$data[[2]][,k],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
                    model[[i]]$estimates_bernoulli[[k]] <- sga_def_dropout(y=model[[i]]$data[[1]][,k],B_mu=B_mu,B_gamma=B_gamma,beta=0,alpha=0,phi=phi,nu=1,N=N,b=b,Db=Db,Db_inv=Db_inv,correct_for_zeros=TRUE,
                                                                           noise_distribution="Bernoulli",dropout_mu=model[[i]]$cv_bernoulli$opt_sd_mu,dropout_gamma=model[[i]]$cv_bernoulli$opt_sd_gamma,
                                                                           learning_rate=list(method="ADADELTA",rho=0.95),batch_size=30,max_iterations=50000,tolerance=10^-5,progress=TRUE,time=FALSE,
                                                                           save_stopping_crit=FALSE)
                    model[[i]]$estimates_normal[[k]] <- sga_def_dropout(y=model[[i]]$data[[1]][,k],B_mu=B_mu,B_gamma=B_gamma,beta=0,alpha=0,phi=phi,nu=1,N=N,b=b,Db=Db,Db_inv=Db_inv,correct_for_zeros=TRUE,
                                                                           noise_distribution="Normal",dropout_mu=model[[i]]$cv_normal$opt_sd_mu,dropout_gamma=model[[i]]$cv_normal$opt_sd_gamma,
                                                                           learning_rate=list(method="ADADELTA",rho=0.95),batch_size=30,max_iterations=50000,tolerance=10^-5,progress=TRUE,time=FALSE,
                                                                           save_stopping_crit=FALSE)
                }
            }
        }
        
        # new cross-validation
        print("Redo cross-validation for benchmark...")
        seq_mu <- c(0,5000)
        seq_gamma <- c(0,5000)
        B_mu <- cbind(rep(1,samplesizes[i]),ns(x=model[[i]]$data[[2]][,1],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
        B_gamma <- ns(x=model[[i]]$data[[2]][,1],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
        if (data_distribution=="Normal") {
            cv_benchmark <- NormalBothGCV(y=model[[i]]$data[[1]][,1],Bmu=B_mu,Bgamma=B_gamma,seqM=seq_mu,seqG=seq_gamma,
                                          lenM=50,lenG=50,minlM=10^-5,maxlM=15000,minlG=10^-5,maxlG=15000,
                                          plotM=FALSE,plotG=FALSE)
        } else if (data_distribution=="Poisson") {
            cv_benchmark <- poisson2sGCV(y=model[[i]]$data[[1]][,1],Bmu=B_mu,Bgamma=B_gamma,seqM=seq_mu,seqG=seq_gamma,
                                         lenM=50,lenG=50,minlM=10^-5,maxlM=15000,minlG=10^-5,maxlG=15000,
                                         plotM=FALSE,plotG=FALSE)
        } else if (data_distribution=="Binomial") {
            prob <- model[[i]]$data[[1]][,1]/N
            cv_benchmark <- BinomialBothGCV(y=prob,num=N,Bmu=B_mu,Bgamma=B_gamma,seqM=seq_mu,seqG=seq_gamma,
                                            lenM=50,lenG=50,minlM=10^-5,maxlM=15000,minlG=10^-5,maxlG=15000,
                                            plotM=FALSE,plotG=FALSE)
        }
        model[[i]] <- append(model[[i]],list(cv_benchmark=cv_benchmark),which(names(model[[i]])=="cv_normal"))
        
        # new benchmarks
        print("Recompute the benchmarks...")
        benchmark_estimates <- foreach (k=1:replicates) %dopar% {
            B_mu <- cbind(rep(1,samplesizes[i]),ns(x=model[[i]]$data[[2]][,k],knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
            B_gamma <- ns(x=model[[i]]$data[[2]][,k],knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
            if (data_distribution=="Normal") {
                fisher_def_normal(y=model[[i]]$data[[1]][,k],Bmu=B_mu,Bgamma=B_gamma,alphaM=0,alphaG=0,phi=phi,lambdaM=model[[i]]$cv_benchmark$vals[1],
                                  lambdaG=model[[i]]$cv_benchmark$vals[2],maxitG=35,maxitT=2000,penordM=2,penordG=2,tol=10^-8)
            } else if (data_distribution=="Poisson") {
                fisher_def_poisson(y=model[[i]]$data[[1]][,k],Bmu=B_mu,Bgamma=B_gamma,alphaM=0,alphaG=0,lambdaM=model[[i]]$cv_benchmark$vals[1],
                                   lambdaG=model[[i]]$cv_benchmark$vals[2],maxitG=100,maxitT=2000,penordM=2,penordG=2,tol=10^-8)
            } else if (data_distribution=="Binomial") {
                prob <- model[[i]]$data[[1]][,k]/N
                fisher_def_binomial(y=prob,num=N,Bmu=B_mu,Bgamma=B_gamma,lambdaM=model[[i]]$cv_benchmark$vals[1],
                                    lambdaG=model[[i]]$cv_benchmark$vals[2],penordM=2,penordG=2,tol=10^-10,type="quasi",
                                    alphaM=0,alphaG=0,tolSTOP=10^-2,linktype="exp",maxitT=2000,maxitM=60,maxitG=45)
            }
        }
        model[[i]]$benchmarks <- benchmark_estimates
        
        # new RMSEs
        print("Recompute all RMSEs...")
        grid <- seq(from=0,to=1,length.out=501)
        
        rmse_mean_est_bernoulli <- foreach(k=1:replicates, .combine=c) %dopar% {
            mean_est_bernoulli <- get_est(x,knots_mu,model[[i]]$estimates_bernoulli[[k]]$beta,Db_current,TRUE)
            get_rmse(truemean,mean_est_bernoulli,grid)
        }
        
        rmse_disp_est_bernoulli <- foreach(k=1:replicates, .combine=c) %dopar% {
            disp_est_bernoulli <- get_est(x,knots_gamma,model[[i]]$estimates_bernoulli[[k]]$alpha,function(x){exp(x)})
            get_rmse(truedisp,disp_est_bernoulli,grid)
        }
        
        rmse_mean_est_normal <- foreach(k=1:replicates, .combine=c) %dopar% {
            mean_est_normal <- get_est(x,knots_mu,model[[i]]$estimates_normal[[k]]$beta,Db_current,TRUE)
            get_rmse(truemean,mean_est_normal,grid)
        }
        
        rmse_disp_est_normal <- foreach(k=1:replicates, .combine=c) %dopar% {
            disp_est_normal <- get_est(x,knots_gamma,model[[i]]$estimates_normal[[k]]$alpha,function(x){exp(x)})
            get_rmse(truedisp,disp_est_normal,grid)
        }
        
        rmse_mean_bench <- foreach(k=1:replicates, .combine=c) %dopar% {
            mean_bench <- get_est(x,knots_mu,model[[i]]$benchmarks[[k]]$alphaM,Db_current,TRUE)
            get_rmse(truemean,mean_bench,grid)
        }
        
        rmse_disp_bench <- foreach(k=1:replicates, .combine=c) %dopar% {
            disp_bench <- get_est(x,knots_gamma,model[[i]]$benchmarks[[k]]$alphaG,function(x){phi*exp(x)^{-1}})
            get_rmse(truedisp,disp_bench,grid)
        }
        
        model[[i]]$rmse_mean_est_bernoulli <- rmse_mean_est_bernoulli
        model[[i]]$rmse_disp_est_bernoulli <- rmse_disp_est_bernoulli
        model[[i]]$rmse_mean_est_normal <- rmse_mean_est_normal
        model[[i]]$rmse_disp_est_normal <- rmse_disp_est_normal
        model[[i]]$rmse_mean_bench <- rmse_mean_bench
        model[[i]]$rmse_disp_bench <- rmse_disp_bench
    }
    
    print(Sys.time()-start_time)
    
    model
}
