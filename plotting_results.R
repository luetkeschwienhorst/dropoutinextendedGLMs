
# removing what is not important
rm(list = ls(all = TRUE))

# load stuff
source("functions.R")
source("helpers.R")

# required packages
packages <- c("splines","MASS","ggplot2","SparseM","lattice","gridBase","grid",
              "caret","parallel","mgcv","viridis","ggthemes","scales","pryr","future",
              "rstudioapi","gamlss","inline","foreach","doParallel","glmnet",
              "scales","haven")

# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# candidates
# order: west, south, east, north, first inbound then outbound
spandauer_damm_ost_2vR <- 100101010044425
spandauer_damm_west_2vR <- 100101010044223

mariendorfer_damm_nord_2vR <- 100101010002793
mariendorfer_damm_sued_2vR <- 100101010079080

maerkische_allee_sued_3vR <- 100101010040179
maerkische_allee_nord_3vR <- 100101010040482

roederallee_sued_2vR <- 100101010004211
roederallee_nord_2vR <- 100101010004009

# get subset of trafficdata
ids <- c(spandauer_damm_ost_2vR,spandauer_damm_west_2vR,
         mariendorfer_damm_nord_2vR,mariendorfer_damm_sued_2vR,
         maerkische_allee_sued_3vR,maerkische_allee_nord_3vR,
         roederallee_sued_2vR,roederallee_nord_2vR)

cols <- c("red","magenta","blue","green")
positions <- c("west","south","east","north")
directions <- c("out","in")

# to be filled
#B_mus <- list()
#B_gammas <- list()
#results_ber <- list()
#results_norm <- list()
#datasets <- list()

# join all results in the same order as in ids
#B_mus <- c(B_mus, list(B_mu))
#B_gammas <- c(B_gammas, list(B_gamma))

# ONLY ONCE!!!
#results_ber <- c(results_ber, list(res))
#results_norm <- c(results_norm, list(res))
#datasets <- c(datasets, list(data))

# add newest
keep <- list("B_mus", "B_gammas", "results_ber", "results_norm", "datasets")
rm(list=setdiff(ls(),keep))

# turn cv_values into vector
#for (i in 1:8) {
#  results_ber[[i]]$cv_values <- unlist(results_ber[[i]]$cv_values)
#  results_norm[[i]]$cv_values <- unlist(results_norm[[i]]$cv_values)
#}

# remove NANs and upper 10% of cv_values for Bernoulli dropout in order
up_to <- 4800
for(i in 1:8) {
  # remove NaNs
  cv <- results_ber[[i]]$cv_values
  grid_mu <- results_ber[[i]]$grid_dropout_mu[!is.na(cv)]
  grid_gamma <- results_ber[[i]]$grid_dropout_gamma[!is.na(cv)]
  cv <- cv[!is.na(cv)]
  # sort according to cv
  grid_mu <- grid_mu[order(cv)]
  grid_gamma <- grid_gamma[order(cv)]
  cv <- sort(cv)
  # only keep the lowest 4750
  cv <- cv[1:up_to]
  grid_mu <- grid_mu[1:up_to]
  grid_gamma <- grid_gamma[1:up_to]
  # keep changes
  results_ber[[i]]$cv_values <- cv
  results_ber[[i]]$grid_dropout_mu <- grid_mu
  results_ber[[i]]$grid_dropout_gamma <- grid_gamma
}

for (i in 1:8) {
  id <- ids[i]
  data <- datasets[[i]]
  cv_bernoulli <- results_ber[[i]]
  cv_gaussian <- results_norm[[i]]
  pos <- positions[ceiling(i/2)]
  dir <- directions[i%%2+1]
  
  # scatter plot bernoulli
  scatter_bernoulli <- data.frame(p_mu=cv_bernoulli$grid_dropout_mu,p_gamma=cv_bernoulli$grid_dropout_gamma,cv=cv_bernoulli$cv_values)
  optimum <- data.frame(p_mu=cv_bernoulli$opt_dropout_mu,p_gamma=cv_bernoulli$opt_dropout_gamma)
  title <- paste(pos,dir,"scatter_ber.png",sep="_")
  png(title,width=600,height=500,res=130,pointsize=10)
  scatter <- ggplot(scatter_bernoulli, aes(x = p_mu, y = p_gamma, color = cv)) + geom_point(size=2,alpha=0.3) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=p_mu,y=p_gamma), color="black", size=2)
  print(scatter)
  dev.off()
  
  # scatter plot normal
  scatter_normal <- data.frame(sd_mu=cv_gaussian$grid_dropout_mu,sd_gamma=cv_gaussian$grid_dropout_gamma,cv=cv_gaussian$cv_values)
  optimum <- data.frame(sd_mu=cv_gaussian$opt_dropout_mu,sd_gamma=cv_gaussian$opt_dropout_gamma)
  title <- paste(pos,dir,"scatter_norm.png",sep="_")
  png(title,width=600,height=500,res=130,pointsize=10)
  scatter <- ggplot(scatter_normal, aes(x = sd_mu, y = sd_gamma, color = cv)) + geom_point(size=2,alpha=0.3) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=sd_mu,y=sd_gamma), color="black", size=2)
  print(scatter)
  dev.off()
}

# plot results

b_dist <- b_poisson
Db_dist <- Db_poisson
Db_inv_dist <- Db_inv_poisson
D2b_dist <- D2b_poisson

# get mean and disp estimate
m_mu <- 15
m_gamma <- 10
x <- seq(0,23,by=0.01)
knots_mu <- get_knots_by_quantiles(x,m_mu)
knots_gamma <- get_knots_by_quantiles(x,m_gamma)
B_mu <- cbind(rep(1,length(x)), cSplineDes(x=x, knots=knots_mu, ord=4))
B_gamma <- cSplineDes(x=x, knots=knots_gamma, ord=4)

# inbound
meanest_ber_1 <- Db_dist(B_mu%*%results_ber[[1]]$opt_beta)
meanest_ber_3 <- Db_dist(B_mu%*%results_ber[[3]]$opt_beta)
meanest_ber_5 <- Db_dist(B_mu%*%results_ber[[5]]$opt_beta)
meanest_ber_7 <- Db_dist(B_mu%*%results_ber[[7]]$opt_beta)
dispest_ber_1 <- Db_dist(B_gamma%*%results_ber[[1]]$opt_alpha)
dispest_ber_3 <- Db_dist(B_gamma%*%results_ber[[3]]$opt_alpha)
dispest_ber_5 <- Db_dist(B_gamma%*%results_ber[[5]]$opt_alpha)
dispest_ber_7 <- Db_dist(B_gamma%*%results_ber[[7]]$opt_alpha)

meanest_norm_1 <- Db_dist(B_mu%*%results_norm[[1]]$opt_beta)
meanest_norm_3 <- Db_dist(B_mu%*%results_norm[[3]]$opt_beta)
meanest_norm_5 <- Db_dist(B_mu%*%results_norm[[5]]$opt_beta)
meanest_norm_7 <- Db_dist(B_mu%*%results_norm[[7]]$opt_beta)
dispest_norm_1 <- Db_dist(B_gamma%*%results_norm[[1]]$opt_alpha)
dispest_norm_3 <- Db_dist(B_gamma%*%results_norm[[3]]$opt_alpha)
dispest_norm_5 <- Db_dist(B_gamma%*%results_norm[[5]]$opt_alpha)
dispest_norm_7 <- Db_dist(B_gamma%*%results_norm[[7]]$opt_alpha)

# outbound
meanest_ber_2 <- Db_dist(B_mu%*%results_ber[[2]]$opt_beta)
meanest_ber_4 <- Db_dist(B_mu%*%results_ber[[4]]$opt_beta)
meanest_ber_6 <- Db_dist(B_mu%*%results_ber[[6]]$opt_beta)
meanest_ber_8 <- Db_dist(B_mu%*%results_ber[[8]]$opt_beta)
dispest_ber_2 <- Db_dist(B_gamma%*%results_ber[[2]]$opt_alpha)
dispest_ber_4 <- Db_dist(B_gamma%*%results_ber[[4]]$opt_alpha)
dispest_ber_6 <- Db_dist(B_gamma%*%results_ber[[6]]$opt_alpha)
dispest_ber_8 <- Db_dist(B_gamma%*%results_ber[[8]]$opt_alpha)

meanest_norm_2 <- Db_dist(B_mu%*%results_norm[[2]]$opt_beta)
meanest_norm_4 <- Db_dist(B_mu%*%results_norm[[4]]$opt_beta)
meanest_norm_6 <- Db_dist(B_mu%*%results_norm[[6]]$opt_beta)
meanest_norm_8 <- Db_dist(B_mu%*%results_norm[[8]]$opt_beta)
dispest_norm_2 <- Db_dist(B_gamma%*%results_norm[[2]]$opt_alpha)
dispest_norm_4 <- Db_dist(B_gamma%*%results_norm[[4]]$opt_alpha)
dispest_norm_6 <- Db_dist(B_gamma%*%results_norm[[6]]$opt_alpha)
dispest_norm_8 <- Db_dist(B_gamma%*%results_norm[[8]]$opt_alpha)

# inbound
linewidth <- 1.5
title <- paste("cv_inbound",sep="",".png")
png(title,width=1600,height=1000,res=300,pointsize=8)
par(mfrow=c(1,2),oma = c(0, 0, 1, 0))
# mean
plot(x,meanest_ber_1,type="l",lwd=linewidth,col=cols[1],lty="dotted",bty="l",ylab="y (cars)",xlab="x (hour)",
     ylim=range(meanest_ber_1,meanest_ber_3,meanest_ber_5,meanest_ber_7,meanest_norm_1,meanest_norm_3,meanest_norm_5,meanest_norm_7),main="Mean")
grid(col="grey",lty = "solid",lwd = 0.5)
lines(x,meanest_ber_3,col=cols[2],type="l",lty="dotted",lwd=linewidth)
lines(x,meanest_ber_5,col=cols[3],type="l",lty="dotted",lwd=linewidth)
lines(x,meanest_ber_7,col=cols[4],type="l",lty="dotted",lwd=linewidth)
lines(x,meanest_norm_1,col=cols[1],type="l",lty="dashed",lwd=linewidth)
lines(x,meanest_norm_3,col=cols[2],type="l",lty="dashed",lwd=linewidth)
lines(x,meanest_norm_5,col=cols[3],type="l",lty="dashed",lwd=linewidth)
lines(x,meanest_norm_7,col=cols[4],type="l",lty="dashed",lwd=linewidth)
legend("topright",legend=c("bernoulli","gaussian"),lty=c("dotted","dashed"),
       col="black",lwd=linewidth, cex=1, box.lty=0,bg="transparent")
# dispersion
plot(x,dispest_ber_1,type="l",lwd=linewidth,col=cols[1],lty="dotted",bty="l",ylab="",xlab="x",
     ylim=range(dispest_ber_1,dispest_ber_3,dispest_ber_5,dispest_ber_7,dispest_norm_1,dispest_norm_3,dispest_norm_5,dispest_norm_7),main="Dispersion")
grid(col="grey",lty = "solid",lwd = 0.5)
lines(x,dispest_ber_3,col=cols[2],type="l",lty="dotted",lwd=linewidth)
lines(x,dispest_ber_5,col=cols[3],type="l",lty="dotted",lwd=linewidth)
lines(x,dispest_ber_7,col=cols[4],type="l",lty="dotted",lwd=linewidth)
lines(x,dispest_norm_1,col=cols[1],type="l",lty="dashed",lwd=linewidth)
lines(x,dispest_norm_3,col=cols[2],type="l",lty="dashed",lwd=linewidth)
lines(x,dispest_norm_5,col=cols[3],type="l",lty="dashed",lwd=linewidth)
lines(x,dispest_norm_7,col=cols[4],type="l",lty="dashed",lwd=linewidth)
legend("topright",legend=c("bernoulli","gaussian"),lty=c("dotted","dashed"),
       col="black",lwd=linewidth, cex=1, box.lty=0,bg="transparent")
#mtext("traffic detection data", line=-1.5, outer = TRUE, cex = 2)
dev.off()
par(mfrow=c(1,1))

# outbound
linewidth <- 1.5
title <- paste("cv_outbound",sep="",".png")
png(title,width=1600,height=1000,res=300,pointsize=8)
par(mfrow=c(1,2),oma = c(0, 0, 1, 0))
# mean
plot(x,meanest_ber_2,type="l",lwd=linewidth,col=cols[1],lty="dotted",bty="l",ylab="y (cars)",xlab="x (hour)",
     ylim=range(meanest_ber_2,meanest_ber_4,meanest_ber_6,meanest_ber_8,meanest_norm_2,meanest_norm_4,meanest_norm_6,meanest_norm_8),main="Mean")
grid(col="grey",lty = "solid",lwd = 0.5)
lines(x,meanest_ber_4,col=cols[2],type="l",lty="dotted",lwd=linewidth)
lines(x,meanest_ber_6,col=cols[3],type="l",lty="dotted",lwd=linewidth)
lines(x,meanest_ber_8,col=cols[4],type="l",lty="dotted",lwd=linewidth)
lines(x,meanest_norm_2,col=cols[1],type="l",lty="dashed",lwd=linewidth)
lines(x,meanest_norm_4,col=cols[2],type="l",lty="dashed",lwd=linewidth)
lines(x,meanest_norm_6,col=cols[3],type="l",lty="dashed",lwd=linewidth)
lines(x,meanest_norm_8,col=cols[4],type="l",lty="dashed",lwd=linewidth)
legend("topright",legend=c("bernoulli","gaussian"),lty=c("dotted","dashed"),
       col="black",lwd=linewidth, cex=1, box.lty=0,bg="transparent")
# dispersion
plot(x,dispest_ber_2,type="l",lwd=linewidth,col=cols[1],lty="dotted",bty="l",ylab="",xlab="x",
     ylim=range(dispest_ber_2,dispest_ber_4,dispest_ber_6,dispest_ber_8,dispest_norm_2,dispest_norm_4,dispest_norm_6,dispest_norm_8),main="Dispersion")
grid(col="grey",lty = "solid",lwd = 0.5)
lines(x,dispest_ber_4,col=cols[2],type="l",lty="dotted",lwd=linewidth)
lines(x,dispest_ber_6,col=cols[3],type="l",lty="dotted",lwd=linewidth)
lines(x,dispest_ber_8,col=cols[4],type="l",lty="dotted",lwd=linewidth)
lines(x,dispest_norm_2,col=cols[1],type="l",lty="dashed",lwd=linewidth)
lines(x,dispest_norm_4,col=cols[2],type="l",lty="dashed",lwd=linewidth)
lines(x,dispest_norm_6,col=cols[3],type="l",lty="dashed",lwd=linewidth)
lines(x,dispest_norm_8,col=cols[4],type="l",lty="dashed",lwd=linewidth)
legend("topright",legend=c("bernoulli","gaussian"),lty=c("dotted","dashed"),
       col="black",lwd=linewidth, cex=1, box.lty=0,bg="transparent")
#mtext("traffic detection data", line=-1.5, outer = TRUE, cex = 2)
dev.off()
par(mfrow=c(1,1))


