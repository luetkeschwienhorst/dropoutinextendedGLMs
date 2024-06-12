
# required packages
library(splines)
library(MASS)
library(ggplot2)
library(SparseM)
library(lattice)
library(gridBase)
library(grid)
library(caret)
library(parallel)
library(mgcv)
library(viridis)
library(ggthemes)
library(scales)
library(pryr)
library(future)
library(rstudioapi)
library(gamlss)
library(inline)
library(foreach)
library(doParallel)

# source functions
source("functions.R")

# bases

n <- 500
X <-seq(0,1,1/(n-1))
m_mu <- 30
m_gamma <- 20
knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
B_mu <- ns(x=X,knots=knots_mu,intercept=FALSE,Boundary.knots=c(0,1))
B_gamma <- ns(x=X,knots=knots_gamma,intercept=FALSE,Boundary.knots=c(0,1))

# filename needs to be a path such as "/home/.../basis_mean.png"
png(filename, width=1600, height=600, res=150, pointsize=13)
plot(B_mu[,1]~X, ylim=range(min(B_mu),max(B_mu)), type='l', lwd=3, col="red", 
     xlab="X (30 knots)", ylab="Y",bty="l")
for (j in 2:ncol(B_mu)) lines(B_mu[,j]~X, lwd=3, col="red")
grid(col="grey",lty = "solid",lwd = 0.3)
dev.off()

# filename needs to be a path such as "/home/.../basis_disp.png"
png(filename, width=1600, height=600, res=150, pointsize=13)
plot(B_gamma[,1]~X, ylim=range(min(B_gamma),max(B_gamma)), type='l', lwd=3, col="blue", 
     xlab="X (20 knots)", ylab="Y",bty="l")
for (j in 2:ncol(B_gamma)) lines(B_gamma[,j]~X, lwd=3, col="blue")
grid(col="grey",lty = "solid",lwd = 0.3)
dev.off()

# plot all testfunctions
# filename needs to be a path such as "/home/.../testfunctions_all.png"
png(filename,width=1600, height=500,res=220,pointsize=12)
par(mfrow=c(1,3))
par(mar=c(3.5,2.5,1,2.5))
curve(f_1,type='l',bty="l",col="red",lwd=2,ylab="",xlab="",ylim=range(f_1(seq(0,1,by=1/100)),f_2(seq(0,1,by=1/100))),main="Mean (Normal)",panel.first=grid(col="grey",lty = "solid",lwd = 0.5))
title(xlab="x",line=2,cex.lab=1.2)
curve(f_2,col="darkorange",lwd=2,add=TRUE)
legend("topright", legend=c("f_1","f_2"),col=c("red","darkorange"),lwd=2,lty=1, cex=1, box.lty=0,bg="transparent")
curve(f_3,type='l',bty="l",col="deepskyblue",lwd=2,ylab="",xlab="",ylim=range(f_3(seq(0,1,by=1/100)),f_4(seq(0,1,by=1/100))),main="Mean (Poisson, Binomial)",panel.first=grid(col="grey",lty = "solid",lwd = 0.5))
title(xlab="x",line=2,cex.lab=1.2)
curve(f_4,col="cyan",lwd=2,add=TRUE)
legend("topright", legend=c("f_3","f_4"),col=c("deepskyblue","cyan"),lwd=2,lty=1, cex=1, box.lty=0,bg="transparent")
curve(g_1,type='l',bty="l",col="darkseagreen",lwd=2,ylab="",xlab="",
      ylim=range(g_1(seq(0,1,by=1/100)),g_2(seq(0,1,by=1/100)),g_3(seq(0,1,by=1/100)),g_4(seq(0,1,by=1/100))),main="Dispersion",panel.first=grid(col="grey",lty = "solid",lwd = 0.5))
title(xlab="x",line=2,cex.lab=1.2)
curve(g_2,col="darkgreen",lwd=2,add=TRUE)
curve(g_3,col="forestgreen",lwd=2,add=TRUE)
curve(g_4,col="green",lwd=2,add=TRUE)
legend("topleft", legend=c("g_1","g_2","g_3","g_4"),col=c("darkseagreen","darkgreen","forestgreen","green"),lwd=2,lty=1, cex=1, box.lty=0,bg="transparent")
dev.off()
par(mfrow=c(1,1))

# plot camparison across certain models

#main=paste("Dispersion ",disps_indices[j])

plot_simulation_comparison <- function (list_of_simulations,location,data_distribution,truemean,mean_index,Db,truedisps,disps_indices,
                                        m_mu,m_gamma,phi,resolution_scatter=130,resolution_functions=150,pointsize=12,width=1600,height=800,
                                        color_bernoulli="firebrick1",color_normal="goldenrod",color_pmle="green4",legend_mean="topright",
                                        legend_disp="topright",lwd_percentiles=2) {
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
        
        # bernoulli
        title <- paste(location,sep="","comparison_bernoulli_",data_distribution,"_",mean_index,"_",samplesize,".png")
        png(title,width=width,height=height,res=resolution_functions,pointsize=pointsize)
        par(mfrow=c(2,n_simulations))
        for (j in 1:n_simulations) {
            sim <- list_of_simulations[[j]][[i]]
            range_mean <- range(truemean(X))
            for (k in 1:replicates) {
                range_mean <- range(range_mean,range(sim$estimates_bernoulli[[k]]$meanest))
            }
            par(mar=c(2.5,2.5,2.5,2.5))
            plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean,main=paste("M_",sep="",mean_index,disps_indices[j]))
            grid(col="grey",lty = "solid",lwd = 0.3)
            for (k in 1:replicates) {
                lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
            }
            lines(X,truemean(X),col="black",type="l",lwd=2)
            if (j==n_simulations) {
              legend(legend_mean, legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
            }
        }
        for (j in 1:n_simulations) {
            sim <- list_of_simulations[[j]][[i]]
            range_disp <- range(truedisps[[j]](X))
            for (k in 1:replicates) {
                range_disp <- range(range_disp,range(sim$estimates_bernoulli[[k]]$dispest))
            }
            par(mar=c(3.5,2.5,1,2.5))
            plot(X,truedisps[[j]](X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp)
            title(xlab="x",line=2,cex.lab=1.2)
            grid(col="grey",lty = "solid",lwd = 0.3)
            for (k in 1:replicates) {
              lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
            }
            lines(X,truedisps[[j]](X),col="black",type="l",lwd=2)
            if (j==n_simulations) {
              legend(legend_disp, legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
              
            }
        }
        dev.off()
        
        #normal
        title <- paste(location,sep="","comparison_normal_",data_distribution,"_",mean_index,"_",samplesize,".png")
        png(title,width=width,height=height,res=resolution_functions,pointsize=pointsize)
        par(mfrow=c(2,n_simulations))
        for (j in 1:n_simulations) {
          sim <- list_of_simulations[[j]][[i]]
          range_mean <- range(truemean(X))
          for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_normal[[k]]$meanest))
          }
          par(mar=c(2.5,2.5,2.5,2.5))
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean,main=paste("M_",sep="",mean_index,disps_indices[j]))
          grid(col="grey",lty = "solid",lwd = 0.3)
          for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_normal[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
          }
          lines(X,truemean(X),col="black",type="l",lwd=2)
          if (j==n_simulations) {
            legend(legend_mean, legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
          }        }
        for (j in 1:n_simulations) {
          sim <- list_of_simulations[[j]][[i]]
          range_disp <- range(truedisps[[j]](X))
          for (k in 1:replicates) {
            range_disp <- range(range_disp,range(sim$estimates_normal[[k]]$dispest))
          }
          par(mar=c(3.5,2.5,1,2.5))
          plot(X,truedisps[[j]](X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp)
          title(xlab="x",line=2,cex.lab=1.2)
          grid(col="grey",lty = "solid",lwd = 0.3)
          for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_normal[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
          }
          lines(X,truedisps[[j]](X),col="black",type="l",lwd=2)
          if (j==n_simulations) {
            legend(legend_disp, legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
            
          }        }
        dev.off()
        
        #pmle
        title <- paste(location,sep="","comparison_pmle_",data_distribution,"_",mean_index,"_",samplesize,".png")
        png(title,width=width,height=height,res=resolution_functions,pointsize=pointsize)
        par(mfrow=c(2,n_simulations))
        for (j in 1:n_simulations) {
          sim <- list_of_simulations[[j]][[i]]
          range_mean <- range(truemean(X))
          for (k in 1:replicates) {
            range_mean <- range(range_mean,range(sim$estimates_pmle[[k]]$meanest))
          }
          par(mar=c(2.5,2.5,2.5,2.5))
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean,main=paste("M_",sep="",mean_index,disps_indices[j]))
          grid(col="grey",lty = "solid",lwd = 0.3)
          for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_pmle[[k]]$meanest,col=alpha("red",0.7),type="l",lwd=0.5)
          }
          lines(X,truemean(X),col="black",type="l",lwd=2)
          if (j==n_simulations) {
            legend(legend_mean, legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
          }        }
        for (j in 1:n_simulations) {
          sim <- list_of_simulations[[j]][[i]]
          range_disp <- range(truedisps[[j]](X))
          for (k in 1:replicates) {
            range_disp <- range(range_disp,range(sim$estimates_pmle[[k]]$dispest))
          }
          par(mar=c(3.5,2.5,1,2.5))
          range_disp <- c(max(-10,range_disp[1]),min(10,range_disp[2]))
          plot(X,truedisps[[j]](X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp)
          title(xlab="x",line=2,cex.lab=1.2)
          grid(col="grey",lty = "solid",lwd = 0.3)
          for (k in 1:replicates) {
            lines(sim$data[[2]][,k+1],sim$estimates_pmle[[k]]$dispest,col=alpha("forestgreen",0.7),type="l",lwd=0.5)
          }
          lines(X,truedisps[[j]](X),col="black",type="l",lwd=2)
          if (j==n_simulations) {
            legend(legend_disp, legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
          }
        }
        dev.off()
    }
}

plot_heatmaps <- function(list_of_simulations,location,data_distribution,truemean,mean_index,truedisps,disps_indices,resolution_scatter) {
  for (j in 1:length(list_of_simulations)) {
    simulation <- list_of_simulations[[j]]
    for (i in 1:length(simulation)) {
      sim <- simulation[[i]]
      replicates <- sim$replicates-1

      # scatter plot bernoulli
      scatter_bernoulli <- data.frame(p_mu=sim$cv_bernoulli$grid_dropout_mu,p_gamma=sim$cv_bernoulli$grid_dropout_gamma,cv=sim$cv_bernoulli$cv_values)
      optimum <- data.frame(p_mu=sim$cv_bernoulli$opt_dropout_mu,p_gamma=sim$cv_bernoulli$opt_dropout_gamma)
      title <- paste(location,sep="","scatter_",data_distribution,"_",mean_index,disps_indices[j],"_bernoulli_",sim$samplesize,".png")
      png(title,width=500,height=500,res=resolution_scatter,pointsize=10)
      scatter <- ggplot(scatter_bernoulli, aes(x = p_mu, y = p_gamma, color = cv)) + geom_point(size=3,alpha=0.8) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=p_mu,y=p_gamma), color="black", size=4) + ggtitle(paste("n = ",sep="",sim$samplesize)) + theme(plot.title = element_text(hjust = 0.5))
      print(scatter)
      dev.off()
      
      # scatter plot normal
      scatter_normal <- data.frame(sd_mu=sim$cv_normal$grid_dropout_mu,sd_gamma=sim$cv_normal$grid_dropout_gamma,cv=sim$cv_normal$cv_values)
      optimum <- data.frame(sd_mu=sim$cv_normal$opt_dropout_mu,sd_gamma=sim$cv_normal$opt_dropout_gamma)
      title <- paste(location,sep="","scatter_",data_distribution,"_",mean_index,disps_indices[j],"_normal_",sim$samplesize,".png")
      png(title,width=500,height=500,res=resolution_scatter,pointsize=10)
      scatter <- ggplot(scatter_normal, aes(x = sd_mu, y = sd_gamma, color = cv)) + geom_point(size=3,alpha=0.8) + scale_color_gradient(low = "#3399FF", high = "#ff0000") + theme_classic() + geom_point(data=optimum, aes(x=sd_mu,y=sd_gamma), color="black", size=4) + ggtitle(paste("n = ",sep="",sim$samplesize)) + theme(plot.title = element_text(hjust = 0.5))
      print(scatter)
      dev.off()
    }
  }
}

plot_boxplots <- function(list_of_simulations,location,data_distribution,truemean,mean_index,disps_indices,resolution=120,pointsize=12,width=400,height=400) {
  replicates <- 100
  for (j in 1:length(list_of_simulations)) {
    simulation <- list_of_simulations[[j]]
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
                               Method=rep(rep(c("Bernoulli","Gaussian","PMLE"),each=replicates),length(simulation)))
    
    title <- paste(location,sep="","boxplot_mean_",data_distribution,"_",mean_index,disps_indices[j],".png")
    png(title,width=width,height=height,res=resolution,pointsize=pointsize)
    if (j==length(list_of_simulations)) {
      boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position=c(0.85,0.85),legend.title = element_text( size=8), legend.text=element_text(size=8),plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("") + ylab("")
    } else if (j==1) {
      boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("") + ylab("RMSE")
    } else {
      boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("") + ylab("")
    }
    print(boxplot_mean)
    dev.off()
    
    title <- paste(location,sep="","boxplot_disp_",data_distribution,"_",mean_index,disps_indices[j],".png")
    png(title,width=width,height=height,res=resolution,pointsize=pointsize)
    if (j==1) {
      boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("RMSE") + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95)))
    } else {
      boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("") + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95)))
    }
    print(boxplot_disp)
    dev.off()
  }
}

pfad <- "/home/.../folder/"

# normal
# mean 1
models <- list(model_11,model_12,model_13,model_14)
plot_simulation_comparison(list_of_simulations=models,location=pfad,data_distribution="normal",
                           truemean=f_1,mean_index=1,Db=Db_normal,
                           truedisps=testfunctions_disp,disps_indices=c(1,2,3,4),m_mu=30,m_gamma=20,phi=0.8^2,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1600,height=700)
plot_heatmaps(models,pfad,"normal",f_1,1,testfunctions_disp,c(1,2,3,4),160)
plot_boxplots(models,pfad,"normal",f_1,1,c(1,2,3,4),150)
# mean 2
models <- list(model_21,model_22,model_23,model_24)
plot_simulation_comparison(list_of_simulations=models,location=pfad,data_distribution="normal",
                           truemean=f_2,mean_index=2,Db=Db_normal,
                           truedisps=testfunctions_disp,disps_indices=c(1,2,3,4),m_mu=30,m_gamma=20,phi=0.8^2,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1600,height=700)
plot_heatmaps(models,pfad,"normal",f_2,2,testfunctions_disp,c(1,2,3,4),160)
plot_boxplots(models,pfad,"normal",f_2,2,c(1,2,3,4),150)



# poisson
# mean 3
models <- list(model_31_poisson,model_32_poisson,model_33_poisson,model_34_poisson)
plot_simulation_comparison(list_of_simulations=models,location=pfad,data_distribution="poisson",
                           truemean=f_3,mean_index=3,Db=Db_poisson,
                           truedisps=testfunctions_disp,disps_indices=c(1,2,3,4),m_mu=30,m_gamma=20,phi=1,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1600,height=700)
plot_heatmaps(models,pfad,"poisson",f_3,3,testfunctions_disp,c(1,2,3,4),160)
plot_boxplots(models,pfad,"poisson",f_3,3,c(1,2,3,4),150)
# mean 4
models <- list(model_41_poisson,model_42_poisson,model_43_poisson,model_44_poisson)
plot_simulation_comparison(list_of_simulations=models,location=pfad,data_distribution="poisson",
                           truemean=f_4,mean_index=4,Db=Db_poisson,
                           truedisps=testfunctions_disp,disps_indices=c(1,2,3,4),m_mu=30,m_gamma=20,phi=1,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1600,height=700)
plot_heatmaps(models,pfad,"poisson",f_4,4,testfunctions_disp,c(1,2,3,4),160)
plot_boxplots(models,pfad,"poisson",f_4,4,c(1,2,3,4),150)



# binomial
# mean 3
N <- 70
models <- list(model_31_binomial,model_32_binomial,model_33_binomial,model_34_binomial)
plot_simulation_comparison(list_of_simulations=models,location=pfad,data_distribution="binomial",
                           truemean=f_3,mean_index=3,Db=function(x){Db_binomial(x,N)},
                           truedisps=testfunctions_disp,disps_indices=c(1,2,3,4),m_mu=30,m_gamma=20,phi=1,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1600,height=700)
plot_heatmaps(models,pfad,"binomial",f_3,3,testfunctions_disp,c(1,2,3,4),160)
plot_boxplots(models,pfad,"binomial",f_3,3,c(1,2,3,4),150)


# mean 4
N <- 30
models <- list(model_41_binomial,model_42_binomial,model_43_binomial,model_44_binomial)
plot_simulation_comparison(list_of_simulations=models,location=pfad,data_distribution="binomial",
                           truemean=f_4,mean_index=4,Db=function(x){Db_binomial(x,N)},
                           truedisps=testfunctions_disp,disps_indices=c(1,2,3,4),m_mu=30,m_gamma=20,phi=1,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1600,height=700)
plot_heatmaps(models,pfad,"binomial",f_4,4,testfunctions_disp,c(1,2,3,4),160)
plot_boxplots(models,pfad,"binomial",f_4,4,c(1,2,3,4),150)

