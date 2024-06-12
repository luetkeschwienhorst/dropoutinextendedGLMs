plot_spline_matrix <- function(B,x) {
  plot(B[,1]~x, ylim=c(0,max(B)), type='l', lwd=2, col=1, xlab=" B-spline basis", ylab="")
  for (j in 2:ncol(B)) lines(B[,j]~x, lwd=2, col=j)
}

def_deviance <- function(y,N=NULL,meanest,phi,b,Db_inv,correct_for_zeros=FALSE,correct_for_N=FALSE){
  if (correct_for_zeros && isTRUE(is.element(0,y))) {
    y[which(y==0)] <- 1e-100
  }
  if (correct_for_N) {
    y[which(y==N)] <- y[which(y==N)]-(1e-10)
  }
  if (correct_for_N) {
    theta <- Db_inv(meanest,N)
    theta_saturated <- Db_inv(y,N)
    deviance <- ((y*theta_saturated-b(theta_saturated,N))/phi-(y*theta-b(theta,N))/phi)
  } else {
    theta <- Db_inv(meanest)
    theta_saturated <- Db_inv(y)
    deviance <- ((y*theta_saturated-b(theta_saturated))/phi-(y*theta-b(theta))/phi)
  }
  deviance
}

def_loglikelihood <- function(y,N,meanest,dispest,phi,nu,b,Db_inv,
                              correct_for_zeros=FALSE,correct_for_N=FALSE) {
  if (correct_for_zeros && isTRUE(is.element(0,y))) {
    y[which(y==0)] <- 1e-100
  }
  if (correct_for_N) {
    y[which(y==N)] <- y[which(y==N)]-(1e-10)
    loglikelihood <- sum((1/2)*log(dispest)
                         +(1/phi)*dispest*(y*Db_inv(meanest,N)-b(Db_inv(meanest,N),N))*nu
                         +(1/phi)*(1-dispest)*(y*Db_inv(y,N)-b(Db_inv(y,N),N))*nu)
  } else {
    loglikelihood <- sum((1/2)*log(dispest)
                         +(1/phi)*dispest*(y*Db_inv(meanest)-b(Db_inv(meanest)))*nu
                         +(1/phi)*(1-dispest)*(y*Db_inv(y)-b(Db_inv(y)))*nu)
  }
  loglikelihood
}

def_log_score <- function(y,N,meanest,dispest,phi,nu,b,Db_inv,
                          correct_for_zeros=FALSE,correct_for_N=FALSE) {
  loglikelihood <- def_loglikelihood(y,N,meanest,dispest,phi,nu,b,Db_inv,
                                     correct_for_zeros,correct_for_N)
  log_score <- -mean(loglikelihood)
  log_score
}

LS <- function(y,N,meanest,dispest,phi,nu,b,Db_inv,
               correct_for_zeros=FALSE,correct_for_N=FALSE) {
  loglikelihood <- def_loglikelihood(y,N,meanest,dispest,phi,nu,b,Db_inv,
                                     correct_for_zeros,correct_for_N)
  log_score <- -loglikelihood
  log_score
}

def_quadratic_score <- function(y_test,N_test,meanest_test,dispest_test,y,N,meanest,dispest,
                                phi,nu,b,Db_inv,correct_for_zeros=FALSE,correct_for_N=FALSE) {
  loglikelihood_test <- def_loglikelihood(y_test,N_test,meanest_test,dispest_test,phi,nu,b,Db_inv,
                                          correct_for_zeros,correct_for_N)
  densities_test <- exp(loglikelihood_test)
  loglikelihood <- def_loglikelihood(y,N,meanest,dispest,phi,nu,b,Db_inv,
                                     correct_for_zeros,correct_for_N)
  densities <- exp(loglikelihood)
  quadratic_score <- 2*mean(densities_test)-sum(densities^2)
  quadratic_score
}

def_spherical_score <- function(y_test,N_test,meanest_test,dispest_test,y,N,meanest,dispest,
                                phi,nu,b,Db_inv,correct_for_zeros=FALSE,correct_for_N=FALSE) {
  loglikelihood_test <- def_loglikelihood(y_test,N_test,meanest_test,dispest_test,phi,nu,b,Db_inv,
                                          correct_for_zeros,correct_for_N)
  densities_test <- exp(loglikelihood_test)
  loglikelihood <- def_loglikelihood(y,N,meanest,dispest,phi,nu,b,Db_inv,
                                     correct_for_zeros,correct_for_N)
  densities <- exp(loglikelihood)
  spherical_score <- mean(densities_test)*(mean(densities^2))^-0.5
  spherical_score
}

def_dawid_sebastiani_score <- function(y,N,meanest,dispest,phi,nu,b,Db_inv,D2b,
                                       correct_for_zeros=FALSE,correct_for_N=FALSE) {
  if (correct_for_zeros && isTRUE(is.element(0,y))) {
    y[which(y==0)] <- 1e-100
  }
  if (correct_for_N) {
    y[which(y==N)] <- y[which(y==N)]-(1e-10)
    varest <- (phi*D2b(Db_inv(meanest,N),N))/(dispest*nu)
  } else {
    varest <- (phi*D2b(Db_inv(meanest)))/(dispest*nu)
  }
  log_varest <- log(varest)
  dawid_sebastiani_score <- mean(log_varest+(y-meanest)^2/varest)
  dawid_sebastiani_score
}

DSS <- function(y,N,meanest,dispest,phi,nu,b,Db_inv,D2b,
                                       correct_for_zeros=FALSE,correct_for_N=FALSE) {
  if (correct_for_zeros && isTRUE(is.element(0,y))) {
    y[which(y==0)] <- 1e-100
  }
  if (correct_for_N) {
    y[which(y==N)] <- y[which(y==N)]-(1e-10)
    varest <- (phi*D2b(Db_inv(meanest,N),N))/(dispest*nu)
  } else {
    varest <- (phi*D2b(Db_inv(meanest)))/(dispest*nu)
  }
  log_varest <- log(varest)
  dawid_sebastiani_score <- log_varest+(y-meanest)^2/varest
  dawid_sebastiani_score
}
