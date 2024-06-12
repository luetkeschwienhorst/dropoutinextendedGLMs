
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

source("functions.R")

# location needs to be given as a path

# plot relevant testfunctions
png(filename=paste(location,sep="","testfunctions.png"),width=1600, height=500,res=220,pointsize=12)
par(mfrow=c(1,3))
par(mar=c(3.5,2.5,1,2.5))
curve(f_1,type='l',bty="l",col="red",lwd=2,ylab="",xlab="",ylim=range(f_1(seq(0,1,by=1/100))),main="Mean (Normal)",panel.first=grid(col="grey",lty = "solid",lwd = 0.5))
title(xlab="x",line=2,cex.lab=1.2)
legend("topright", legend=c("f_1"),col=c("red"),lwd=2,lty=1, cex=1, box.lty=0,bg="transparent")
curve(f_3,type='l',bty="l",col="deepskyblue",lwd=2,ylab="",xlab="",ylim=range(f_3(seq(0,1,by=1/100))),main="Mean (Poisson, Binomial)",panel.first=grid(col="grey",lty = "solid",lwd = 0.5))
title(xlab="x",line=2,cex.lab=1.2)
legend("topright", legend=c("f_2"),col=c("deepskyblue"),lwd=2,lty=1, cex=1, box.lty=0,bg="transparent")
curve(g_4,type='l',bty="l",col="orange",lwd=2,ylab="",xlab="",
      ylim=range(g_2(seq(0,1,by=1/100)),g_3(seq(0,1,by=1/100)),g_4(seq(0,1,by=1/100))),main="Dispersion",panel.first=grid(col="grey",lty = "solid",lwd = 0.5))
title(xlab="x",line=2,cex.lab=1.2)
curve(g_2,col="blueviolet",lwd=2,add=TRUE)
curve(g_3,col="green",lwd=2,add=TRUE)
curve(g_4,col="orange",lwd=2,add=TRUE)
legend("topleft", legend=c("g_1","g_2","g_3"),col=c("orange","blueviolet","green"),lwd=2,lty=1, cex=1, box.lty=0,bg="transparent")
dev.off()
par(mfrow=c(1,1))

# colors for methods
color_bernoulli <- "firebrick1"
color_normal <- "forestgreen"
color_pmle <- "dodgerblue"

plot_simulation_by_method <- function (simulation,location,data_distribution,truemean,mean_index,Db,truedisp,disp_index,
                                        m_mu,m_gamma,phi,resolution_scatter=130,resolution_functions=150,pointsize=12,width=1600,height=400,
                                        color_bernoulli="firebrick1",color_normal="forestgreen",color_pmle="dodgerblue",titled=FALSE,x=FALSE,lwd_percentiles=2) {
    n_samplesizes <- length(simulation)
    replicates <- simulation[[1]]$replicates-1
    for (i in 1:n_samplesizes) {
      
        sim <- simulation[[i]]
        samplesize <- simulation[[i]]$samplesize
        X <- seq(from=0,to=1,length.out=samplesize)
        knots_mu <- seq(0,1-1/(m_mu-1),1/(m_mu-1))
        knots_gamma <- seq(0,1-1/(m_gamma-1),1/(m_gamma-1))
        B_mu <- cbind(rep(1,samplesize),ns(x=X,knots=knots_mu,intercept = FALSE,Boundary.knots=c(0,1)))
        B_gamma <- ns(x=X,knots=knots_gamma,intercept = FALSE,Boundary.knots=c(0,1))
        
        # comparison of methods for given samplesize
        title <- paste(location,sep="","mean_comparison_",data_distribution,"_",mean_index,disp_index,"_",samplesize,".png")
        png(title,width=width,height=height,res=resolution_functions,pointsize=pointsize)
        par(mfrow=c(1,3))
        
        # bernoulli mean
        range_mean <- range(truemean(X))
        for (k in 1:replicates) {
          range_mean <- range(range_mean,range(sim$estimates_bernoulli[[k]]$meanest))
        }
        if (x) {
          if (titled) {
            par(mar=c(3.5,3.5,2.5,2.5))
          } else {
            par(mar=c(3.5,3.5,0,2.5))
          }
        } else {
          par(mar=c(2.5,3.5,2.5,2.5))
        }
        if (titled) {
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean,main="Bernoulli dropout")
        } else {
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean)
        }
        if (x) {
          title(xlab="x",line=2,cex.lab=1.2)
        }
        title(ylab="mean",line=2,cex.lab=1.2)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
          lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$meanest,col=alpha(color_bernoulli,0.7),type="l",lwd=0.5)
        }
        lines(X,truemean(X),col="black",type="l",lwd=2)
        
        # gaussian mean
        range_mean <- range(truemean(X))
        for (k in 1:replicates) {
          range_mean <- range(range_mean,range(sim$estimates_normal[[k]]$meanest))
        }
        if (x) {
          par(mar=c(3.5,2.5,2.5,2.5))
        } else {
          par(mar=c(2.5,2.5,2.5,2.5))
        }
        if (titled) {
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean,main="Gaussian dropout")
        } else {
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean)
        }
        if (x) {
          title(xlab="x",line=2,cex.lab=1.2)
        }
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
          lines(sim$data[[2]][,k+1],sim$estimates_normal[[k]]$meanest,col=alpha(color_normal,0.7),type="l",lwd=0.5)
        }
        lines(X,truemean(X),col="black",type="l",lwd=2)
        
        # pmle mean
        range_mean <- range(truemean(X))
        for (k in 1:replicates) {
          range_mean <- range(range_mean,range(sim$estimates_pmle[[k]]$meanest))
        }
        if (x) {
          par(mar=c(3.5,2.5,2.5,2.5))
        } else {
          par(mar=c(2.5,2.5,2.5,2.5))
        }
        if (titled) {
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean,main="PMLE")
        } else {
          plot(X,truemean(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_mean)
        }
        if (x) {
          title(xlab="x",line=2,cex.lab=1.2)
        }
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
          lines(sim$data[[2]][,k+1],sim$estimates_pmle[[k]]$meanest,col=alpha(color_pmle,0.7),type="l",lwd=0.5)
        }
        lines(X,truemean(X),col="black",type="l",lwd=2)
        #legend(legend_mean, legend=c("true mean","estimates"),col=c("black",alpha("red",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        dev.off()
        
        title <- paste(location,sep="","disp_comparison_",data_distribution,"_",mean_index,disp_index,"_",samplesize,".png")
        png(title,width=width,height=height,res=resolution_functions,pointsize=pointsize)
        par(mfrow=c(1,3))
        
        # bernoulli disp
        range_disp <- range(truedisp(X))
        for (k in 1:replicates) {
          range_disp <- range(range_disp,range(sim$estimates_bernoulli[[k]]$dispest),range(sim$estimates_normal[[k]]$dispest),range(sim$estimates_pmle[[k]]$dispest))
        }
        if (x) {
          par(mar=c(3.5,3.5,2.5,2.5))
        } else {
          par(mar=c(2.5,3.5,2.5,2.5))
        }
        if (titled) {
          plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp,main="Bernoulli dropout")
        } else {
          plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp)
        }
        if (x) {
          title(xlab="x",line=2,cex.lab=1.2)
        }
        title(ylab="dispersion",line=2,cex.lab=1.2)
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
          lines(sim$data[[2]][,k+1],sim$estimates_bernoulli[[k]]$dispest,col=alpha(color_bernoulli,0.7),type="l",lwd=0.5)
        }
        lines(X,truedisp(X),col="black",type="l",lwd=2)
        
        # gaussian disp
        if (x) {
          par(mar=c(3.5,2.5,2.5,2.5))
        } else {
          par(mar=c(2.5,2.5,2.5,2.5))
        }
        if (titled) {
          plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp,main="Gaussian dropout")
        } else {
          plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp)
        }
        if (x) {
          title(xlab="x",line=2,cex.lab=1.2)
        }
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
          lines(sim$data[[2]][,k+1],sim$estimates_normal[[k]]$dispest,col=alpha(color_normal,0.7),type="l",lwd=0.5)
        }
        lines(X,truedisp(X),col="black",type="l",lwd=2)
        
        # pmle disp
        if (x) {
          par(mar=c(3.5,2.5,2.5,2.5))
        } else {
          par(mar=c(2.5,2.5,2.5,2.5))
        }
        if (titled) {
          plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp,main="PMLE")
        } else {
          plot(X,truedisp(X),bty="l",col="black",type="l",lwd=2,ylab="",xlab="",ylim=range_disp)
        }
        if (x) {
          title(xlab="x",line=2,cex.lab=1.2)
        }
        grid(col="grey",lty = "solid",lwd = 0.3)
        for (k in 1:replicates) {
          lines(sim$data[[2]][,k+1],sim$estimates_pmle[[k]]$dispest,col=alpha(color_pmle,0.7),type="l",lwd=0.5)
        }
        lines(X,truedisp(X),col="black",type="l",lwd=2)
        #legend(legend_disp, legend=c("true dispersion","estimates"),col=c("black",alpha("forestgreen",0.7)),lty=1,lwd=2, cex=1, box.lty=0,bg="transparent")
        dev.off()
    }
}

plot_boxplots_neurips_0 <- function(simulation,legend=FALSE,RMSE=FALSE,location,data_distribution,truemean,mean_index,disp_index,resolution=120,pointsize=12,width=400,height=400) {
  replicates <- 100
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
  
  title <- paste(location,sep="","boxplot_mean_",data_distribution,"_",mean_index,disp_index,".png")
  png(title,width=width,height=height,res=resolution,pointsize=pointsize)
  if (legend) {
    if (RMSE) {
      todo
    } else {
      todo
    }
    
  }
  boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position=c(0.85,0.85),legend.title = element_text( size=8), legend.text=element_text(size=8),plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("samplesize") + ylab("RMSE")
  print(boxplot_mean)
  dev.off()
  
  title <- paste(location,sep="","boxplot_disp_",data_distribution,"_",mean_index,disp_index,".png")
  png(title,width=width,height=height,res=resolution,pointsize=pointsize)
  boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("") + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95)))
  print(boxplot_disp)
  dev.off()
}

plot_boxplots_neurips_1 <- function(list_of_simulations,location,data_distributions,mean_indices,disps_indices,resolution=120,pointsize=12,width=400,height=400) {
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
    
    title <- paste(location,sep="","boxplot_mean_",data_distributions[[j]],"_",mean_indices[j],disps_indices[j],".png")
    png(title,width=width,height=height,res=resolution,pointsize=pointsize)
    if (j==length(list_of_simulations)) {
      boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position=c(0.85,0.85),legend.title = element_text( size=8), legend.text=element_text(size=8),plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("") + ylab("") + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
    } else if (j==1) {
      boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("") + ylab("RMSE") + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
    } else {
      boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Mean") + xlab("") + ylab("") + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
    }
    print(boxplot_mean)
    dev.off()
    
    title <- paste(location,sep="","boxplot_disp_",data_distributions[[j]],"_",mean_indices[j],disps_indices[j],".png")
    png(title,width=width,height=height,res=resolution,pointsize=pointsize)
    if (j==1) {
      boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("RMSE") + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95))) + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
    } else {
      boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5)) + ggtitle("Dispersion") + xlab("samplesize") + ylab("") + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95))) + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
    }
    print(boxplot_disp)
    dev.off()
  }
}

plot_boxplots_neurips_2 <- function(simulation,location,data_distribution,plot_title,x,y,legend=FALSE,mean_index,disp_index,resolution=120,pointsize=12,width=400,height=400) {
  replicates <- 100
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
  
  title <- paste(location,sep="","boxplot_mean_",data_distribution,"_",mean_index,disp_index,".png")
  if (legend) {
    png(title,width=1.4*width,height=height,res=resolution,pointsize=pointsize)
  } else {
    png(title,width=width,height=height,res=resolution,pointsize=pointsize)
  }
  if (legend) {
    boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.title = element_text(size=8), legend.text=element_text(size=8), plot.title = element_text(face = "bold",hjust = 0.5), plot.margin = margin(0, 0, 0, 0, "cm")) + ggtitle(plot_title) + xlab(x) + ylab(y) + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
  } else {
    boxplot_mean <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_mean, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size = 0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5), plot.margin = margin(0, 0, 0, 0, "cm")) + ggtitle(plot_title) + xlab(x) + ylab(y) + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
  }
  print(boxplot_mean)
  dev.off()
  
  title <- paste(location,sep="","boxplot_disp_",data_distribution,"_",mean_index,disp_index,".png")
  if (legend) {
    png(title,width=1.4*width,height=height,res=resolution,pointsize=pointsize)
  } else {
    png(title,width=width,height=height,res=resolution,pointsize=pointsize)
  }
  if (legend) {
    boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.title = element_text(size=8), legend.text=element_text(size=8), plot.title = element_text(face = "bold",hjust = 0.5), plot.margin = margin(0, 0, 0, 0, "cm")) + ggtitle(plot_title) + xlab(x) + ylab(y) + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95))) + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
  } else {
    boxplot_disp <- ggplot(boxplot_data,aes(x=factor(samplesize,levels=samplesizes), y=RMSE_disp, fill=Method)) + stat_boxplot(geom= 'errorbar' , width = 0.2, position = position_dodge(width = 0.75) ) + geom_boxplot(outlier.colour = alpha("black",alpha=0.7), outlier.shape = 20, outlier.size=0.7) + theme_linedraw() + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = 0.5), plot.margin = margin(0, 0, 0, 0, "cm")) + ggtitle(plot_title) + xlab(x) + ylab(y) + scale_y_continuous(limits = quantile(boxplot_data$RMSE_disp, c(0, 0.95))) + scale_fill_manual(values=c(color_bernoulli,color_normal,color_pmle))
  }
  print(boxplot_disp)
  dev.off()
}

# boxplots
plot_boxplots_neurips_2(model_14,location,"normal","Normal","","RMSE",FALSE,1,4,150)
plot_boxplots_neurips_2(model_12,location,"normal","","","RMSE",FALSE,1,2,150)
plot_boxplots_neurips_2(model_13,location,"normal","","samplesize","RMSE",FALSE,1,3,150)

plot_boxplots_neurips_2(model_34_poisson,location,"poisson","Poisson","","",FALSE,3,4,150)
plot_boxplots_neurips_2(model_32_poisson,location,"poisson","","","",FALSE,3,2,150)
plot_boxplots_neurips_2(model_33_poisson,location,"poisson","","samplesize","",FALSE,3,3,150)

plot_boxplots_neurips_2(model_34_binomial,location,"binomial","Binomial","","",TRUE,3,4,150)
plot_boxplots_neurips_2(model_32_binomial,location,"binomial","","","",FALSE,3,2,150)
plot_boxplots_neurips_2(model_33_binomial,location,"binomial","","samplesize","",FALSE,3,3,150)

# plots of estimates
# normal
plot_simulation_by_method(model_14,location=location,data_distribution="normal",
                           truemean=f_1,mean_index=1,Db=Db_normal,
                           truedisp=g_4,disp_index=4,m_mu=30,m_gamma=20,phi=0.8^2,
                           resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=TRUE,x=FALSE)

plot_simulation_by_method(model_12,location=location,data_distribution="normal",
                          truemean=f_1,mean_index=1,Db=Db_normal,
                          truedisp=g_2,disp_index=2,m_mu=30,m_gamma=20,phi=0.8^2,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=FALSE,x=FALSE)

plot_simulation_by_method(model_13,location=location,data_distribution="normal",
                          truemean=f_1,mean_index=1,Db=Db_normal,
                          truedisp=g_3,disp_index=3,m_mu=30,m_gamma=20,phi=0.8^2,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=FALSE,x=TRUE)


plot_simulation_by_method(model_34_poisson,location=location,data_distribution="poisson",
                          truemean=f_3,mean_index=3,Db=Db_poisson,
                          truedisp=g_4,disp_index=4,m_mu=30,m_gamma=20,phi=1,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=TRUE,x=FALSE)

plot_simulation_by_method(model_32_poisson,location=location,data_distribution="poisson",
                          truemean=f_3,mean_index=3,Db=Db_poisson,
                          truedisp=g_2,disp_index=2,m_mu=30,m_gamma=20,phi=1,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=FALSE,x=FALSE)

plot_simulation_by_method(model_33_poisson,location=location,data_distribution="poisson",
                          truemean=f_3,mean_index=3,Db=Db_poisson,
                          truedisp=g_3,disp_index=3,m_mu=30,m_gamma=20,phi=1,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=FALSE,x=TRUE)


plot_simulation_by_method(model_34_binomial,location=location,data_distribution="binomial",
                          truemean=f_3,mean_index=3,Db=function(x){Db_binomial(x,70)},
                          truedisp=g_4,disp_index=4,m_mu=30,m_gamma=20,phi=1,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=TRUE,x=FALSE)

plot_simulation_by_method(model_32_binomial,location=location,data_distribution="binomial",
                          truemean=f_3,mean_index=3,Db=function(x){Db_binomial(x,70)},
                          truedisp=g_2,disp_index=2,m_mu=30,m_gamma=20,phi=1,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=FALSE,x=FALSE)

plot_simulation_by_method(model_33_binomial,location=location,data_distribution="binomial",
                          truemean=f_3,mean_index=3,Db=function(x){Db_binomial(x,70)},
                          truedisp=g_3,disp_index=3,m_mu=30,m_gamma=20,phi=1,
                          resolution_scatter=130,resolution_functions=150,pointsize=14,width=1200,height=400,titled=FALSE,x=TRUE)
