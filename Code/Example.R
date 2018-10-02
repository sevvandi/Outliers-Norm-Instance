#---- Set working directory to the folder which has the files
#---- 1. Feature_Functions.R and
#---- 2. yeast_C2_P05_V05_CNA.arff
setwd("Folder with the above two files")

source('./Feature_Functions.R', echo=TRUE)

library("RWeka")
library("infotheo")
library("quantmod")
library("FNN")
library("igraph")
library("ks")


dat <- read.arff(paste("yeast_C2_P05_V05_CNA.arff", sep=""))
features <- ComputeMetaFeatures(dat)
