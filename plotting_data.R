###################
# plot traffic data
###################

# removing stuff
rm(list=setdiff(ls(), "test_data"))

# load dropout and PMLE estimators
source("functions.R")
source("GijbelsProsdocimiClaeskens.R")
source("helpers.R")

# required packages
library("ggplot2")
packages <- c("splines","MASS","ggplot2","SparseM","lattice","gridBase","grid",
              "caret","parallel","mgcv","viridis","ggthemes","scales","pryr","future",
              "rstudioapi","gamlss","inline","foreach","doParallel","glmnet",
              "scales","haven",)

# install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# load trafficdata
trafficdata <- read.csv("verkehrsdetektion.csv", header = TRUE)

# candidates
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
sub_trafficdata <- lapply(trafficdata, function(list){list[(trafficdata$detid_15 == ids[1]
                                                            | trafficdata$detid_15 == ids[2]
                                                            | trafficdata$detid_15 == ids[3]
                                                            | trafficdata$detid_15 == ids[4]
                                                            | trafficdata$detid_15 == ids[5]
                                                            | trafficdata$detid_15 == ids[6]
                                                            | trafficdata$detid_15 == ids[7]
                                                            | trafficdata$detid_15 == ids[8])
                                                           & trafficdata$weekend == 0
                                                           & trafficdata$year == 2019
                                                           & (trafficdata$month == 6 | trafficdata$month == 7 | trafficdata$month == 8)
                                                           & !(trafficdata$month==6 & trafficdata$day==10)]})

test_trafficdata <- lapply(trafficdata, function(list){list[(trafficdata$detid_15 == ids[1]
                                                            | trafficdata$detid_15 == ids[2]
                                                            | trafficdata$detid_15 == ids[3]
                                                            | trafficdata$detid_15 == ids[4]
                                                            | trafficdata$detid_15 == ids[5]
                                                            | trafficdata$detid_15 == ids[6]
                                                            | trafficdata$detid_15 == ids[7]
                                                            | trafficdata$detid_15 == ids[8])
                                                           & trafficdata$weekend == 0
                                                           & trafficdata$year == 2019
                                                           & trafficdata$month == 9]})

test_data <- list()
for (i in 1:8) {
  id <- ids[i]
  data <- data.frame(list(x = test_trafficdata$hour[test_trafficdata$detid_15 == id], y = test_trafficdata$q_pkw[test_trafficdata$detid_15 == id]))
  test_data[[i]] <- data
}
id <- ids[1]
data <- data.frame(list(x = sub_trafficdata$hour[sub_trafficdata$detid_15 == id], y = sub_trafficdata$q_pkw[sub_trafficdata$detid_15 == id]))


# plot
for (i in 1:length(ids)) {
  filename <- paste(positions[ceiling(i/2)],sep="","_",directions[i%%2+1],".png")
  png(filename,width=700, height=500,res=110,pointsize=18)
  par(mar = c(2.5, 2.5, 2.5, 2.5))
  plot(sub_trafficdata$hour[sub_trafficdata$detid_15 == ids[i]], sub_trafficdata$q_pkw[sub_trafficdata$detid_15 == ids[i]],
       pch=19, col=alpha(cols[ceiling(i/2)],0.3), bty="l", xlab = "", ylab = "")
  dev.off()
}
