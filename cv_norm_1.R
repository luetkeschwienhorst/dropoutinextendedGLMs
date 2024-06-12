#############################
# cross-validate traffic data
#############################

# clear the environment
rm(list = ls(all = TRUE))
#rm(list=ls()[!grepl("^model",ls())])

# load dropout estimators and some helpers
source("cluster_functions.R")
source("helpers.R")

# required packages
packages <- c("splines","MASS","ggplot2","SparseM","lattice","gridBase","grid",
              "caret","parallel","mgcv","viridis","ggthemes","scales","pryr","future",
              "rstudioapi","gamlss","inline","foreach","doParallel","glmnet",
              "scales","haven", "Rmpi")

# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# load data
load("traffic_detection_data.RData")

# candidates
spandauer_damm_ost_2vR <- 100101010044425
spandauer_damm_west_2vR <- 100101010044223

mariendorfer_damm_nord_2vR <- 100101010002793
mariendorfer_damm_sued_2vR <- 100101010079080

maerkische_allee_sued_3vR <- 100101010040179
maerkische_allee_nord_3vR <- 100101010040482

roederallee_sued_2vR <- 100101010004211
roederallee_nord_2vR <- 100101010004009

ids <- c(spandauer_damm_ost_2vR,spandauer_damm_west_2vR,
         mariendorfer_damm_nord_2vR,mariendorfer_damm_sued_2vR,
         maerkische_allee_sued_3vR,maerkische_allee_nord_3vR,
         roederallee_sued_2vR,roederallee_nord_2vR)

# fix one candidate
id <- ids[1]
data <- data.frame(list(x = sub_trafficdata$hour[sub_trafficdata$detid_15 == id], y = sub_trafficdata$q_pkw[sub_trafficdata$detid_15 == id]))

# parameters
N <- NULL
correct_for_zeros <- TRUE
correct_for_N <- FALSE
learning_rate <- list(method="ADADELTA",rho=0.95)
max_iterations <- 50000
tolerance <- 10^-5
batch_size <- 30 

# cross-validation
n_draws <- 5000
n_folds <- 5

# knots
m_mu <- 15
m_gamma <- 10
knots_mu <- get_knots_by_quantiles(data$x,m_mu)
knots_gamma <- get_knots_by_quantiles(data$x,m_gamma)

# distribution
data_distr <- "poisson"
nu <- 1
phi <- 1
b <- b_poisson
Db <- Db_poisson
Db_inv <- Db_inv_poisson
D2b <- D2b_poisson

# sort data
data <- list_sorter(data,data$x)

# design matrices
y <- data$y
B_mu <- cbind(rep(1,length(data$x)),cSplineDes(x=data$x,knots=knots_mu, ord=4))
B_gamma <- cSplineDes(x=data$x,knots=knots_gamma, ord=4)

# cross-validation
noise_distribution <- "Normal"
interval_dropout_mu <- c(0,2)
interval_dropout_gamma <- c(0,2)
print("Cross-validation for Normal dropout:")
cpus = detectCores()
print(paste("Number of threads: ", cpus))
start.time <- Sys.time()

################################################################################

# check if batch size is small enough
if(floor(length(y)*((n_folds-1)/n_folds))<batch_size){
  message("The batch size of SGA and the fold size of the k-fold cross-validation are incompatible. 
                Adjust accordingly!")
}

# preparations
if(nu==1){
  nu <- rep(1,length(y))
}
if (correct_for_zeros && isTRUE(is.element(0,y))) {
  y[which(y==0)] <- 1e-100
}
if (correct_for_N && isTRUE(is.element(N,y))) {
  y[which(y==N)] <- N-(1e-10)
}
folds <- createFolds(y,k=n_folds)
grid_dropout_mu <- runif(n_draws,min(interval_dropout_mu),max(interval_dropout_mu))
grid_dropout_gamma <- runif(n_draws,min(interval_dropout_gamma),max(interval_dropout_gamma))
drawn_hyperparameters <- mapply(c, grid_dropout_mu, grid_dropout_gamma, SIMPLIFY=FALSE)

# set up cluster
environment <- as.list(ls())
cluster <- makeCluster(cpus)
clusterExport(cluster, environment)

# get cv values for all hyperparameters
cv_values <- parLapply(cl = cluster,
                       X = drawn_hyperparameters,
                       fun = function(x) evaluate_draw(y,B_mu,B_gamma,phi,nu,N,b,Db,Db_inv,D2b,
                                                       correct_for_zeros,correct_for_N,noise_distribution,
                                                       x,learning_rate,batch_size,max_iterations,tolerance,folds))
stopCluster(cluster)

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

################################################################################

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# save data
save.image(file=paste("cv_norm_",id,sep="",".RData"))

