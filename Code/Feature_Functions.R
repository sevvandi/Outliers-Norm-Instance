#####################################################################################
#### All features
#### BEGIN
#####################################################################################
ComputeMetaFeatures <- function(dat.o){
  
  f0 <- ComputeFeaturesNoScaling(dat.o)
  f1 <- ComputeNewSetWithScaling4(dat.o,1)
  f2 <- ComputeNewSetWithScaling4(dat.o,3)
  f3 <- ComputeNewSetWithScaling5(dat.o,1)
  f4 <- ComputeNewSetWithScaling5(dat.o,3)
  f5 <- ComputeNormRelatedFeatures(dat.o)
  
  
  features <- cbind.data.frame(f0,f1,f2,f3,f4, f5) 
  colnames(features) <- c(colnames(f0), colnames(f1), colnames(f2), colnames(f3),colnames(f4),colnames(f5) )
  return(features)
}

#####################################################################################
#### All features
#### END
#####################################################################################


ComputeFeaturesNoScaling <- function(dat.o){
  dat.o <- data.frame(dat.o) 
  
  ## Pre-processing
  dat.2.mid <- dat.o[,colnames(dat.o)[colnames(dat.o)!="id"]]
  dat.1.mid <- dat.o[,colnames(dat.o)[(colnames(dat.o)!="id")& (colnames(dat.o)!="outlier")]]
  if (sum(colnames(dat.o)=="outlier")==0){
    dat.mid <- dat.1.mid[,-1*dim(dat.1.mid)[2]]
  }else{
    dat.mid <- dat.1.mid
  }
  dat.mid <- as.data.frame(dat.mid)
  
  for(i in 1:dim(dat.mid)[2]){
    if(!is.numeric(dat.mid[,i])){
      dat.mid[,i]<- as.numeric(dat.mid[,i])
    }
  }
  
  temp <- apply(dat.mid,2,function(x)ifelse(sd(x)==0,1,0) )
  cols.sd.not.zero <- setdiff(1:dim(dat.mid)[2], which(temp==1))
  dat <- dat.mid[,cols.sd.not.zero]
  cols.sd.not.zero <- setdiff(1:dim(dat.2.mid)[2], which(temp==1))
  dat.2 <- dat.2.mid[,cols.sd.not.zero]
  dat.2 <- as.data.frame(dat.2)
  
  
  n= 24
  feature.vector <- matrix(NA, nrow=1, ncol=n)
  # F1 - Number of observations
  feature.vector[1] <- dim(dat)[1]
  # F2 - Number of attributes
  feature.vector[2] <- dim(dat)[2]
  # F3 - Number of binary attributes
  feature.vector[3] <- NumberOfBinaryAttributes(dat)  
  # F4 - Number of numerical attributes
  feature.vector[4] <- dim(dat)[2] - feature.vector[3]
  # F5 - Ratio of nominal attributes to numerical
  feature.vector[5] <- (feature.vector[3]+1)/(feature.vector[4]+1)
  # F6 - Mean Absolute Correlation of attributes
  feature.vector[6] <- MeanAbsoluteCorrelation(dat)
  # F8 - SD explained by the First PC component
  feature.vector[7] <- PCAVarComp1(dat)
  # F9 - Mean Entropy Of Attributes
  feature.vector[8] <- MeanEntropyOfAttributes(dat)
  # F10 - Entropy Of Total dataset without class and id
  feature.vector[9] <- infotheo::entropy(discretize(dat))
  # F11 - Mean Mutual Information 
  feature.vector[10] <- MeanMutualInformation(dat.2)
  # F12 - Noise to signal ratio
  feature.vector[11] <- SignalToNoiseRatio(dat)
  # F13 - Modality 
  feature.vector[12] <- Modality(dat)[1]
  # F14 - Density disparity 
  feature.vector[13] <-  Modality(dat)[2]
  # Skewness Stats
  feature.vector[14:17] <- Skewness_stats(dat)
  # Kurtosis Stats
  feature.vector[18:21] <- Kurtosis_stats(dat)
  
  iqr.to.sd <- apply(dat,2,IQRtoSDRatio)
  feature.vector[22] <- max(iqr.to.sd)
  feature.vector[23] <- quantile(iqr.to.sd, probs = 0.95)
  feature.vector[24] <- SpreadMeasure3(dat)

  
  colnames(feature.vector)<- c("Num_Observations", "Num_Attributes", "Binary_Attributes",  "Numerical_Attributes", "Ratio_Nomial_Numerical","Mean_Abs_Corr",  "First_PC_SD", "Mean_Entropy_Attr", "Total_Entropy_Dataset", "Mean_Mutual_Info", "SNR", "Modality", "Density_Disparity",  "Skew_Mean", "Skew_Med", "Skew_Max", "Skew_95","Kuto_Mean", "Kuto_Med", "Kuto_Max", "Kuto_95", "MAX_IQR_TO_SD", "IQR_TO_SD_95","Spread_0")
  return(feature.vector)
}




NumberOfBinaryAttributes <- function(dat){
  n.col <- dim(dat)[2]
  z = 0
  for(i in 1:n.col){
    if(length(unique(dat[,i]))<=2){
      z <- z+1
    }
  }
  return(z)
}




MeanAbsoluteCorrelation <- function(dat){
  # This feature computes mean value of correlation between attributes, when taken into account all pairs of attributes
  # z is the number of columns
  z <- dim(dat)[2]
  for(i in 1:z){
    if(!is.numeric(dat[,i])){
      dat[,i]<- as.numeric(dat[,i])
    }
  }
  mean.corr <- (sum(abs(cor(dat)))-z)/(z^2-z)
  return(mean.corr)
}




PCAVarComp1 <- function(dat){
  pca.dat <- prcomp(dat, center = TRUE, scale=TRUE )
  firstpc.sd <- pca.dat$sdev[1]/sum(pca.dat$sdev)
  return(firstpc.sd)
}

MeanEntropyOfAttributes <- function(dat){
  mean.entropy <- mean(apply(discretize(dat),2,infotheo::entropy))
  return(mean.entropy)
}




MeanMutualInformation <- function(dat.2){
  dat <- dat.2[,colnames(dat.2)[colnames(dat.2)!="outlier"]]
  output <- sapply(discretize(dat), mutinformation, Y = as.matrix(dat.2$outlier))
  return(mean(output))
}


SignalToNoiseRatio <- function(dat){
  mean.all <- apply(dat,2, mean)
  sd.all <- apply(dat,2,sd)
  return(mean(mean.all/sd.all))
}






Modality <- function(dat){
  num.col <- dim(dat)[2]
  numeric.list <- NA
  for(i in 1:num.col){
    if(length(unique(dat[,i]))>2){
      numeric.list <- c(numeric.list,i)
    }
  }
  if(length(numeric.list)>2){
    numeric.list <- numeric.list[2:length(numeric.list)]
  }else{
    return(c(0,0))
  }
  
  d <-density(as.matrix(dat[,numeric.list]))
  peak.pos <- quantmod::findPeaks(d$y)
  num.peaks <- length(peak.pos)
  y.vals.1 <- d$y[peak.pos]
  y.vals <- y.vals.1[y.vals.1!=0]
  if(num.peaks>1){
    combos <- combn(y.vals,2)
    max.ratio <- max(combos[1,]/combos[2,], combos[2,]/combos[1,])
  }else{
    max.ratio <- 1
  }
  
  return(c(num.peaks, max.ratio))
}






Skewness_stats <- function(z){
  sk <- moments::skewness(z)
  sk0 <- mean(sk)
  sk1 <- median(sk)
  sk2 <- max(sk)
  sk3 <- quantile(sk,probs=0.95)
  return(c(sk0,sk1, sk2, sk3))
}

IQRtoSDRatio <-  function(z){
  if(sd(z)==0){
    return(0)
  }else{
    return(IQR(z)/sd(z))
  }
}

Kurtosis_stats <- function(z){
  kt <- moments::kurtosis(z)
  kt0 <- mean(kt)
  kt1 <- median(kt)
  kt2 <- max(kt)
  kt3 <- quantile(kt,probs=0.95)
  return(c(kt0,kt1, kt2, kt3))
}




unitize_1 <- function(z) {
  # min-max normalization - 0 -1
  min.z <- min(z)
  max.z <- max(z)
  if ((max.z-min.z)==0)
    return(z)
  (z - min.z )/(max.z-min.z)
}

unitize_2 <- function(z) {
  # Mean and SD normalization
  mean.z <- mean(z)
  sd.z <- sd(z)
  if (sd.z==0)
    return(z)
  (z - mean.z )/sd.z
}

unitize_3 <- function(z) {
  # Median and IQR normalization
  median.z <- median(z)
  iqr.z <- IQR(z)
  if (iqr.z==0)
    return(z)
  (z - median.z )/iqr.z
}


unitize_4 <- function(z) {
  # Median and MAD normalization
  median.z <- median(z)
  mad.z <- 1.4826*median(abs(z-median.z))
  if (mad.z==0)
    return(z)
  (z - median.z )/mad.z
}


SpreadMeasure3 <- function(dat){
  require(dbscan)
  dat <- data.frame(dat)
  nn <- max(floor(dim(dat)[1]/200),10)
  knn.out <- knn.dist(dat,nn)
  epsilon <- quantile(knn.out[,2],probs =0.9)
  dbscan.ex <- dbscan::dbscan(dat, eps=epsilon)
  output <- length(unique(dbscan.ex$cluster))-1
  return(output)
}



#########################################################################################
##########################################################################################
#### FUNCTIONS FOR OUTLIERS SET 1 - DENSITY AND RESIDUAL ANALYSIS                    #####
#### BEGIN                                                                           #####
##########################################################################################

ComputeNewSetWithScaling4 <- function(dat.o,norm_tech){
  require("FNN")
  dat.o <- data.frame(dat.o) 
  
  ## Pre-processing
  dat.2.mid <- dat.o[,colnames(dat.o)[colnames(dat.o)!="id"]]
  dat.1.mid <- dat.o[,colnames(dat.o)[(colnames(dat.o)!="id")& (colnames(dat.o)!="outlier")]]
  dat.outlier <- dat.o[,colnames(dat.o)[colnames(dat.o)=="outlier"]]
  if (sum(colnames(dat.o)=="outlier")==0){
    dat.mid <- dat.1.mid[,-1*dim(dat.1.mid)[2]]
  }else{
    dat.mid <- dat.1.mid
  }
  dat.mid<- as.data.frame(dat.mid)
  
  for(i in 1:dim(dat.mid)[2]){
    if(!is.numeric(dat.mid[,i])){
      dat.mid[,i]<- as.numeric(dat.mid[,i])
    }
  }
  
  temp <- apply(dat.mid,2,function(x)ifelse(sd(x)==0,1,0) )
  cols.sd.not.zero <- setdiff(1:dim(dat.mid)[2], which(temp==1))
  dat <- dat.mid[,cols.sd.not.zero]
  cols.sd.not.zero <- setdiff(1:dim(dat.2.mid)[2], which(temp==1))
  dat.2 <- dat.2.mid[,cols.sd.not.zero]
  dat <- as.data.frame(dat)
  
  if(norm_tech==1){
    ## Min-max normalization
    dat.normed <- apply(dat,2,unitize_1)
  }else if(norm_tech==2){
    ## Mean_sd normalization
    dat.normed <- apply(dat,2,unitize_2)
  }else if(norm_tech==3){
    ## median-IQR normalization
    dat.normed <- apply(dat,2,unitize_3)
  }else{
    ## median-MAD normalization
    dat.normed <- apply(dat,2,unitize_4)
  }
  
  #n= 34
  #feature.vector <- matrix(NA, nrow=1, ncol=n)
  
  out1 <- DBSCANProperties(dat.normed,dat.outlier)
  out2 <-  DensityOfProxiesAndOutliers(dat.normed,dat.outlier)   
  out3 <-  ResidualsOfProxisAndOutliers(dat.normed,dat.outlier)
  feature.vector <- cbind(out1, out2,out3)
  colnames(feature.vector)<- paste(c("OPO_" ), colnames(feature.vector), "_", norm_tech, sep="")
  return(feature.vector)
}





FindNonBinaryVariables <- function(dat){
  un.vals <- apply(dat,2,function(x) length(unique(x)))
  cols <- which(un.vals >2)
  return(cols)
}


DBSCANProperties <- function(dat, oo){
  require("FNN")
  require("dbscan")
  dat <- as.data.frame(dat)
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  z.n <- knn.dist(dat,pref.k)
  z1 <- z.n[,pref.k]
  pot.outliers <- order(z1,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  
  outliers <-which(oo=='yes')
  
  
  k.value <- min(ceiling(dim(dat)[1]/20),50)
  knn.dist <- knn.dist(dat,k.value)
  eps.val <- mean(knn.dist[,k.value])
  if(dim(dat)[2]>10){
    ## Compute PCA and do the first 10 PCs
    pca.obj <- prcomp(dat)
    dbobj <- dbscan::dbscan(pca.obj$x[,1:10],eps=eps.val)
  }else{
    dbobj <- dbscan::dbscan(dat, eps=eps.val)
  }
  
  dbscan.clust.0 <- which(dbobj$cluster==0) 
  set1 <- intersect(outliers, dbscan.clust.0)
  pp <- which(table(dbobj$cluster)<k.value)
  pp.val <- sort(unique(dbobj$cluster))[pp]
  set2 <- intersect(outliers, which(dbobj$cluster  %in% pp.val))
  
  F1 <- (length(set1)+ length(set2))/length(outliers)
  F2 <- length(set1)/length(outliers)
  F3 <- length(set2)/length(outliers)
  
  out <- cbind.data.frame(F1,F2,F3) 
  colnames(out) <- c("DBSCAN_R1","DBSCAN_R2","DBSCAN_R3" )
  return(out)
}

DensityOfProxiesAndOutliers <- function(dat, oo){
  require("FNN")
  dat <- as.data.frame(dat)
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  z.n <- knn.dist(dat,pref.k)
  z1 <- z.n[,pref.k]
  pot.outliers <- order(z1,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  outliers <- which(oo=='yes')
  
  cols1<- which(apply(dat,2,sd)!=0)
  dat1 <- dat[,cols1]
  dat1 <- as.data.frame(dat1)
  pca.obj <- prcomp(dat1)
  if(dim(dat1)[2]>10){
    cols = 1:10
  }else{
    cols = 1:dim(dat1)[2]
  }
  
  
  out <- data.frame(Den_KNOut_Mean=numeric(),Den_KNOut_Median=numeric(),Den_KNOut_SD=numeric(), Den_KNOut_IQR=numeric(), Den_KNOut_Max=numeric(), Den_KNOut_Min=numeric(), Den_KNOut_95P=numeric(), Den_KNOut_05P=numeric(), Den_DenOut_Mean=numeric(),Den_DenOut_Median=numeric(),Den_DenOut_SD=numeric(), Den_DenOut_IQR=numeric(), Den_DenOut_Max=numeric(), Den_DenOut_Min=numeric(), Den_DenOut_95P=numeric(), Den_DenOut_05P=numeric(), DenOut_KNOut_Common=numeric(), Out_KNOut_1=numeric(), Out_KNOut_2=numeric(), Out_DenOut_1=numeric(), Out_DenOut_2=numeric(), Den_Out_Mean=numeric(),Den_Out_Median=numeric(),Den_Out_SD=numeric(), Den_Out_IQR=numeric(), Den_Out_Max=numeric(), Den_Out_Min=numeric(), Den_Out_95P=numeric(), Den_Out_05P=numeric(), DenOut_Out_Mean=numeric(),DenOut_Out_Median=numeric(),DenOut_Out_SD=numeric(), DenOut_Out_IQR=numeric(), DenOut_Out_Max=numeric(), DenOut_Out_Min=numeric(), DenOut_Out_95P=numeric(), DenOut_Out_05P=numeric() )
  
  quant1 <- length(intersect(pot.outliers,outliers))/length(outliers)
  quant2 <- length(intersect(pot.outliers,outliers))/length(pot.outliers)
  
  
  if(length(cols)>1){
    col.pairs <- matrix(c(cols[-length(cols)], cols[-1]), ncol=2)
    # if(dim(col.pairs.1)[1]>50){
    #   col.pairs <- col.pairs.1[1:50,]
    # }else{
    #   col.pairs <- col.pairs.1
    # }
    abmatrix <- matrix(0,nrow=dim(col.pairs)[1],ncol=37)
    do.not.include.rows <- c()
    for(ii in 1:dim(col.pairs)[1]){
      sel.cols <- col.pairs[ii,]
      tryCatch(
        { ## Try part
          temp <- kde(pca.obj$x[,sel.cols],compute.cont=TRUE, eval.points = dat[,sel.cols])$estimate
          abmatrix[ii,1:8] <- ComputeDensityFeatures(temp,pot.outliers)
          pot.outliers.density <- order(temp,decreasing=FALSE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
          abmatrix[ii,9:16] <- ComputeDensityFeatures(temp,pot.outliers.density)
          abmatrix[ii,17] <- length(intersect(pot.outliers,pot.outliers.density))/length(pot.outliers)
          
          ## Measures on outliers
          abmatrix[ii,18] <- quant1
          abmatrix[ii,19] <- quant2
          abmatrix[ii,20] <- length(intersect(pot.outliers.density,outliers))/length(outliers)
          abmatrix[ii,21] <- length(intersect(pot.outliers.density,outliers))/length(pot.outliers.density)
          abmatrix[ii,22:29] <- ComputeDensityFeatures(temp,outliers)
          
          abmatrix[ii,30:37] <- ComparePOWithOut(temp, outliers,1)
          
        },error=function(cond){
          abmatrix[ii,1:17] <- rep(1,17)
          do.not.include.rows <- c(do.not.include.rows,ii)
          message("Here's the original error message:")
          message(cond)
        }, warning=function(cond) {
          message("Here's the original warning message:")
          message(cond)
        }
      )
    }
    if(length(do.not.include.rows)>0){
      out[1,] <- apply(abmatrix[-do.not.include.rows,],2,mean)
    }else{
      out[1,] <- apply(abmatrix,2,mean)
    }
    
  }else{  ### only one column in dat
    temp <- kde(dat[,cols],eval.points = dat[,cols])$estimate
    out1 <- ComputeDensityFeatures(temp,pot.outliers)
    pot.outliers.density <- order(temp,decreasing=FALSE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
    out2 <- ComputeDensityFeatures(temp,pot.outliers.density)
    out3 <- length(intersect(pot.outliers,pot.outliers.density))/length(pot.outliers)
    
    ## Measures on outliers
    quant3 <- length(intersect(pot.outliers.density,outliers))/length(outliers)
    quant4 <- length(intersect(pot.outliers.density,outliers))/length(pot.outliers.density)
    out4 <- ComputeDensityFeatures(temp,outliers)
    out5 <- ComparePOWithOut(temp, outliers,1)
    out[1,] <- c(out1,out2,out3, quant1, quant2,quant3,quant4, out4, out5  )
  }
  ##colnames(out) <- c("Den_KNOut_Mean","Den_KNOut_Median","Den_KNOut_SD", "Den_KNOut_IQR", "Den_KNOut_Max", "Den_KNOut_Min", "Den_KNOut_95P", "Den_KNOut_05P", "Den_DenOut_Mean","Den_DenOut_Median","Den_DenOut_SD", "Den_DenOut_IQR", "Den_DenOut_Max", "Den_DenOut_Min", "Den_DenOut_95P", "Den_DenOut_05P", "DenOut_KNOut_Common" )
  return(out)
}


ComputeDensityFeatures <- function(dx, pot.outliers){
  # dx is the density estimate
  # pot.outliers are the potential outliers 
  
  # F1 - (mean density of non-pot.outliers)/(mean density of pot.outliers) 
  # F2 - (median density of non-pot.outliers)/(median density of pot.outliers) 
  # F3 - (SD density of non-pot.outliers)/(SD density of pot.outliers)
  # F4 - (IQR density of non-pot.outliers)/(IQR percentile density of pot.outliers)
  # F5 - (max density of non-pot.outliers)/(max density of pot.outliers)
  # F6 - (min density of non-pot.outliers)/(min density of pot.outliers)
  # F7 - (95% percentile density of non-pot.outliers)/(95% percentile density of pot.outliers)
  # F8 - (5% percentile density of non-pot.outliers)/(5% percentile density of pot.outliers)
  
  F1.0 <- (mean(dx[-pot.outliers],na.rm = TRUE) +0.01 )/(0.01 + mean(dx[pot.outliers],na.rm = TRUE))
  F2.0 <- (0.01+median(dx[-pot.outliers],na.rm = TRUE))/(0.01+median(dx[pot.outliers],na.rm = TRUE))
  if(length(pot.outliers)==1){
    F3.0 <- (0.01+sd(dx[-pot.outliers],na.rm = TRUE))/(0.01+0)
  }else{
    F3.0 <- (0.01+sd(dx[-pot.outliers],na.rm = TRUE))/(0.01+sd(dx[pot.outliers],na.rm = TRUE))
  }
  
  F4.0 <- (IQR(dx[-pot.outliers],na.rm = TRUE)+0.01)/(0.01+IQR(dx[pot.outliers],na.rm = TRUE))
  F5.0 <- (max(dx[-pot.outliers],na.rm = TRUE)+0.01)/(0.01+max(dx[pot.outliers],na.rm = TRUE))
  F6.0 <- (min(dx[-pot.outliers],na.rm = TRUE)+0.01)/(0.01+min(dx[pot.outliers],na.rm = TRUE))
  F7.0 <- (quantile(dx[-pot.outliers], prob=0.95,na.rm = TRUE)+0.01)/(0.01+quantile(dx[pot.outliers], prob=0.95,na.rm = TRUE))
  F8.0 <- (quantile(dx[-pot.outliers], prob=0.05,na.rm = TRUE)+0.01)/(0.01+quantile(dx[pot.outliers], prob=0.05,na.rm = TRUE))
  
  F1 <- ifelse(mean(dx[pot.outliers],na.rm = TRUE)==0, F1.0, mean(dx[-pot.outliers],na.rm = TRUE)/mean(dx[pot.outliers],na.rm = TRUE))
  F2 <- ifelse(median(dx[pot.outliers],na.rm = TRUE)==0, F2.0, median(dx[-pot.outliers],na.rm = TRUE)/median(dx[pot.outliers],na.rm = TRUE))
  
  if(length(pot.outliers)==1){
    F3 <- F3.0 
  }else{
    F3 <- ifelse(sd(dx[pot.outliers],na.rm = TRUE)==0, F3.0, sd(dx[-pot.outliers],na.rm = TRUE)/sd(dx[pot.outliers],na.rm = TRUE))
  }
  
  
  F4 <- ifelse(IQR(dx[pot.outliers],na.rm = TRUE)==0, F4.0, IQR(dx[-pot.outliers],na.rm = TRUE)/IQR(dx[pot.outliers],na.rm = TRUE))
  F5 <- ifelse(max(dx[pot.outliers],na.rm = TRUE)==0, F5.0, max(dx[-pot.outliers],na.rm = TRUE)/max(dx[pot.outliers],na.rm = TRUE))
  F6 <- ifelse(min(dx[pot.outliers],na.rm = TRUE)==0, F6.0, min(dx[-pot.outliers],na.rm = TRUE)/min(dx[pot.outliers],na.rm = TRUE))
  F7 <- ifelse(quantile(dx[pot.outliers], prob=0.95,na.rm = TRUE)==0, F7.0, quantile(dx[-pot.outliers], prob=0.95)/quantile(dx[pot.outliers], prob=0.95,na.rm = TRUE))
  F8 <- ifelse(quantile(dx[pot.outliers], prob=0.05,na.rm = TRUE)==0, F8.0,quantile(dx[-pot.outliers], prob=0.05,na.rm = TRUE)/quantile(dx[pot.outliers], prob=0.05,na.rm = TRUE))
  
  return(c(F1,F2,F3,F4,F5,F6,F7,F8))
}


ComparePOWithOut<- function(dx, outliers, sw){
  # dx is the density estimate
  # outliers are real outliers 
  # sw is a switch
  # if sw==1 then it is density so order should be increasing
  # if sw==2 then it is residuals so order should be decreasing
  if(sw==1){
    pot.outliers <-order(dx,decreasing=FALSE)[1:length(outliers)]
  }else{
    pot.outliers <-order(dx,decreasing=TRUE)[1:length(outliers)]
  }
  # F1 - (mean density of pot.outliers)/(mean density of outliers) 
  # F2 - (median density of pot.outliers)/(median density of outliers) 
  # F3 - (SD density of pot.outliers)/(SD density of outliers)
  # F4 - (IQR density of pot.outliers)/(IQR percentile density of outliers)
  # F5 - (max density of pot.outliers)/(max density of outliers)
  # F6 - (min density of pot.outliers)/(min density of outliers)
  # F7 - (95% percentile density of pot.outliers)/(95% percentile density of outliers)
  # F8 - (5% percentile density of pot.outliers)/(5% percentile density of outliers)
  
  F1.0 <- (mean(dx[pot.outliers],na.rm = TRUE) +0.01 )/(0.01 + mean(dx[outliers],na.rm = TRUE))
  F2.0 <- (0.01+median(dx[pot.outliers],na.rm = TRUE))/(0.01+median(dx[outliers],na.rm = TRUE))
  if(length(outliers)==1){
    F3.0 <- (0.01+sd(dx[pot.outliers],na.rm = TRUE))/(0.01+0)
  }else{
    F3.0 <- (0.01+sd(dx[pot.outliers],na.rm = TRUE))/(0.01+sd(dx[outliers],na.rm = TRUE))
  }
  
  F4.0 <- (IQR(dx[pot.outliers],na.rm = TRUE)+0.01)/(0.01+IQR(dx[outliers],na.rm = TRUE))
  F5.0 <- (max(dx[pot.outliers],na.rm = TRUE)+0.01)/(0.01+max(dx[outliers],na.rm = TRUE))
  F6.0 <- (min(dx[pot.outliers],na.rm = TRUE)+0.01)/(0.01+min(dx[outliers],na.rm = TRUE))
  F7.0 <- (quantile(dx[pot.outliers], prob=0.95,na.rm = TRUE)+0.01)/(0.01+quantile(dx[outliers], prob=0.95,na.rm = TRUE))
  F8.0 <- (quantile(dx[pot.outliers], prob=0.05,na.rm = TRUE)+0.01)/(0.01+quantile(dx[outliers], prob=0.05,na.rm = TRUE))
  
  F1 <- ifelse(mean(dx[outliers],na.rm = TRUE)==0, F1.0, mean(dx[pot.outliers],na.rm = TRUE)/mean(dx[outliers],na.rm = TRUE))
  F2 <- ifelse(median(dx[outliers],na.rm = TRUE)==0, F2.0, median(dx[pot.outliers],na.rm = TRUE)/median(dx[outliers],na.rm = TRUE))
  
  if(length(outliers)==1){
    F3 <- F3.0 
  }else{
    F3 <- ifelse(sd(dx[outliers],na.rm = TRUE)==0, F3.0, sd(dx[pot.outliers],na.rm = TRUE)/sd(dx[outliers],na.rm = TRUE))
  }
  
  
  F4 <- ifelse(IQR(dx[outliers],na.rm = TRUE)==0, F4.0, IQR(dx[pot.outliers],na.rm = TRUE)/IQR(dx[outliers],na.rm = TRUE))
  F5 <- ifelse(max(dx[outliers],na.rm = TRUE)==0, F5.0, max(dx[pot.outliers],na.rm = TRUE)/max(dx[outliers],na.rm = TRUE))
  F6 <- ifelse(min(dx[outliers],na.rm = TRUE)==0, F6.0, min(dx[pot.outliers],na.rm = TRUE)/min(dx[outliers],na.rm = TRUE))
  F7 <- ifelse(quantile(dx[outliers], prob=0.95,na.rm = TRUE)==0, F7.0, quantile(dx[pot.outliers], prob=0.95)/quantile(dx[outliers], prob=0.95,na.rm = TRUE))
  F8 <- ifelse(quantile(dx[outliers], prob=0.05,na.rm = TRUE)==0, F8.0,quantile(dx[pot.outliers], prob=0.05,na.rm = TRUE)/quantile(dx[outliers], prob=0.05,na.rm = TRUE))
  
  return(c(F1,F2,F3,F4,F5,F6,F7,F8))
  
}



ResidualsOfProxisAndOutliers <- function(dat, oo){
  require("FNN")
  dat <- as.data.frame(dat)
  
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  z.n <- knn.dist(dat,pref.k)
  z1 <- z.n[,pref.k]
  pot.outliers <- order(z1,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  outliers <- which(oo=='yes')
  
  cols <- FindNonBinaryVariables(dat)
  if(length(cols)==0){  ## all columns are binary
    cols = 1:dim(dat)[2]
  }
  
  if(length(cols)==1){  ## only one numeric column in the dataset
    cols = 1:dim(dat)[2]
  }
  
  if(dim(dat)[1]/length(cols)< 10){  ## less than 10 observations per attribute
    num.cols.to.be.chosen <- min(floor(dim(dat)[1]/10), 50)
    vars <- apply(dat[,cols],2,var)    
    cols <- order(vars, decreasing=TRUE)[1:num.cols.to.be.chosen]
  }
  
  if(length(cols)>1){
    set.seed(101)
    cols <- sample(cols,length(cols))
  }
  
  out <- data.frame(Res_KNOut_Mean=numeric(),Res_KNOut_Median=numeric(),Res_KNOut_SD=numeric(), Res_KNOut_IQR=numeric(), Res_KNOut_Max=numeric(), Res_KNOut_Min=numeric(), Res_KNOut_95P=numeric(), Res_KNOut_05P=numeric(), Res_ResOut_Mean=numeric(),Res_ResOut_Median=numeric(),Res_ResOut_SD=numeric(), Res_ResOut_IQR=numeric(), Res_ResOut_Max=numeric(), Res_ResOut_Min=numeric(), Res_ResOut_95P=numeric(), Res_ResOut_05P=numeric(), ResOut_KNOut_Common=numeric(), Out_ResOut_1=numeric(), Out_ResOut_2=numeric(), Res_Out_Mean=numeric(),Res_Out_Median=numeric(),Res_Out_SD=numeric(), Res_Out_IQR=numeric(), Res_Out_Max=numeric(), Res_Out_Min=numeric(), Res_Out_95P=numeric(), Res_Out_05P=numeric(), ResOut_Out_Mean=numeric(),ResOut_Out_Median=numeric(),ResOut_Out_SD=numeric(), ResOut_Out_IQR=numeric(), ResOut_Out_Max=numeric(), ResOut_Out_Min=numeric(), ResOut_Out_95P=numeric(), ResOut_Out_05P=numeric())
  
  
  
  if(length(cols)>1){
    kk <- min(length(cols),50)
    output.mat <- matrix(0, nrow=kk, ncol=35)
    for(ii in 1:kk){
      y <- dat[cols[ii]]
      x <- dat[cols[-ii]]
      model <- lm(unlist(y)~., data=x)
      pot.outliers.residuals <- order(abs(model$residuals),decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
      output.mat[ii,1:8]<- ComputeResidualFeatures(abs(model$residuals), pot.outliers)
      output.mat[ii,9:16]<- ComputeResidualFeatures(abs(model$residuals), pot.outliers.residuals)
      output.mat[ii,17] <- length(intersect(pot.outliers,pot.outliers.residuals))/length(pot.outliers)
      
      ## Measures on outliers
      output.mat[ii,18] <- length(intersect(pot.outliers.residuals,outliers))/length(outliers)
      output.mat[ii,19] <- length(intersect(pot.outliers.residuals,outliers))/length(pot.outliers.residuals)
      output.mat[ii,20:27] <- ComputeResidualFeatures(abs(model$residuals),outliers)
      
      output.mat[ii,28:35] <-  ComparePOWithOut(abs(model$residuals), outliers,2)
    }
    out[1,] <- apply(output.mat,2,mean)
  }else{  ## only one column anyway
    residuals <- dat[,cols]-mean(dat[,cols])
    pot.outliers.residuals <- order(abs(residuals),decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
    out1 <- ComputeResidualFeatures(abs(residuals), pot.outliers)
    out2 <- ComputeResidualFeatures(abs(residuals), pot.outliers.residuals)
    out3 <- length(intersect(pot.outliers,pot.outliers.residuals))/length(pot.outliers)
    
    ## Measures on outliers
    quant3 <- length(intersect(pot.outliers.residuals,outliers))/length(outliers)
    quant4 <- length(intersect(pot.outliers.residuals,outliers))/length(pot.outliers.residuals)
    out4 <- ComputeDensityFeatures(abs(residuals),outliers)
    out5 <-  ComparePOWithOut(abs(residuals), outliers,2)
    out[1,] <- c(out1,out2,out3, quant3,quant4, out4,out5  )
  }
  return(out)
}



ComputeResidualFeatures <- function(res, pot.outliers){
  
  # res is the residuals measure
  # pot.outliers are the potential outliers 
  
  # F1 - (mean density of non-pot.outliers)/(mean density of pot.outliers) 
  # F2 - (median density of non-pot.outliers)/(median density of pot.outliers) 
  # F3 - (SD density of non-pot.outliers)/(SD density of pot.outliers)
  # F4 - (IQR density of non-pot.outliers)/(IQR percentile density of pot.outliers)
  # F5 - (max density of non-pot.outliers)/(max density of pot.outliers)
  # F6 - (min density of non-pot.outliers)/(min density of pot.outliers)
  # F7 - (95% percentile density of non-pot.outliers)/(95% percentile density of pot.outliers)
  # F8 - (5% percentile density of non-pot.outliers)/(5% percentile density of pot.outliers)
  res <- abs(res)
  
  F1.0 <- (0.01 + mean(res[pot.outliers],na.rm = TRUE))/(mean(res[-pot.outliers],na.rm = TRUE) +0.01 )
  F2.0 <- (0.01+median(res[pot.outliers],na.rm = TRUE))/(0.01+median(res[-pot.outliers],na.rm = TRUE))
  
  if(length(pot.outliers)==1){
    F3.0 <- (0.01+0)/(0.01+sd(res[-pot.outliers],na.rm = TRUE))
  }else{
    F3.0 <- (0.01+sd(res[pot.outliers],na.rm = TRUE))/(0.01+sd(res[-pot.outliers],na.rm = TRUE))
  }
  
  
  F4.0 <- (0.01+IQR(res[pot.outliers],na.rm = TRUE))/(IQR(res[-pot.outliers],na.rm = TRUE)+0.01)
  F5.0 <- (0.01+max(res[pot.outliers],na.rm = TRUE))/(max(res[-pot.outliers],na.rm = TRUE)+0.01)
  F6.0 <- (0.01+min(res[pot.outliers],na.rm = TRUE))/(min(res[-pot.outliers],na.rm = TRUE)+0.01)
  F7.0 <- (0.01+quantile(res[pot.outliers], prob=0.95,na.rm = TRUE))/(quantile(res[-pot.outliers], prob=0.95,na.rm = TRUE)+0.01)
  F8.0 <- (0.01+quantile(res[pot.outliers], prob=0.05,na.rm = TRUE))/(quantile(res[-pot.outliers], prob=0.05,na.rm = TRUE)+0.01)
  
  F1 <- ifelse(mean(res[-pot.outliers],na.rm = TRUE)==0, F1.0, mean(res[pot.outliers],na.rm = TRUE)/mean(res[-pot.outliers],na.rm = TRUE))
  F2 <- ifelse(median(res[-pot.outliers],na.rm = TRUE)==0, F2.0, median(res[pot.outliers],na.rm = TRUE)/median(res[-pot.outliers],na.rm = TRUE))
  
  if(length(pot.outliers)==1){
    F3 <- F3.0 
  }else{
    F3 <- ifelse(sd(res[-pot.outliers],na.rm = TRUE)==0, F3.0, sd(res[pot.outliers],na.rm = TRUE)/sd(res[-pot.outliers],na.rm = TRUE))
  }
  
  F4 <- ifelse(IQR(res[-pot.outliers],na.rm = TRUE)==0, F4.0, IQR(res[pot.outliers],na.rm = TRUE)/IQR(res[-pot.outliers],na.rm = TRUE))
  F5 <- ifelse(max(res[-pot.outliers],na.rm = TRUE)==0, F5.0, max(res[pot.outliers],na.rm = TRUE)/max(res[-pot.outliers],na.rm = TRUE))
  F6 <- ifelse(min(res[-pot.outliers],na.rm = TRUE)==0, F6.0, min(res[pot.outliers],na.rm = TRUE)/min(res[-pot.outliers],na.rm = TRUE))
  F7 <- ifelse(quantile(res[-pot.outliers], prob=0.95,na.rm = TRUE)==0, F7.0, quantile(res[pot.outliers], prob=0.95,na.rm = TRUE)/quantile(res[-pot.outliers], prob=0.95,na.rm = TRUE))
  F8 <- ifelse(quantile(res[-pot.outliers], prob=0.05,na.rm = TRUE)==0, F8.0,quantile(res[pot.outliers], prob=0.05,na.rm = TRUE)/quantile(res[-pot.outliers], prob=0.05,na.rm = TRUE))
  
  return(c(F1,F2,F3,F4,F5,F6,F7,F8))
}


#########################################################################################
##########################################################################################
#### FUNCTIONS FOR OUTLIERS SET 1 - DENSITY AND RESIDUAL ANALYSIS                    #####
#### END                                                                             #####
##########################################################################################



##########################################################################################
##########################################################################################
#### FUNCTIONS FOR OUTLIERS SET 2 - LOCAL DENSITY AND IGRAPH FEATURES                #####
#### BEGIN                                                                           #####
##########################################################################################


ComputeNewSetWithScaling5 <- function(dat.o,norm_tech){
  require("FNN")
  dat.o <- data.frame(dat.o) 
  
  ## Pre-processing
  dat.2.mid <- dat.o[,colnames(dat.o)[colnames(dat.o)!="id"]]
  dat.1.mid <- dat.o[,colnames(dat.o)[(colnames(dat.o)!="id")& (colnames(dat.o)!="outlier")]]
  dat.outlier <- dat.o[,colnames(dat.o)[colnames(dat.o)=="outlier"]]
  if (sum(colnames(dat.o)=="outlier")==0){
    dat.mid <- dat.1.mid[,-1*dim(dat.1.mid)[2]]
  }else{
    dat.mid <- dat.1.mid
  }
  dat.mid<- as.data.frame(dat.mid)
  
  for(i in 1:dim(dat.mid)[2]){
    if(!is.numeric(dat.mid[,i])){
      dat.mid[,i]<- as.numeric(dat.mid[,i])
    }
  }
  
  temp <- apply(dat.mid,2,function(x)ifelse(sd(x)==0,1,0) )
  cols.sd.not.zero <- setdiff(1:dim(dat.mid)[2], which(temp==1))
  dat <- dat.mid[,cols.sd.not.zero]
  cols.sd.not.zero <- setdiff(1:dim(dat.2.mid)[2], which(temp==1))
  dat.2 <- dat.2.mid[,cols.sd.not.zero]
  dat <- as.data.frame(dat)
  
  if(norm_tech==1){
    ## Min-max normalization
    dat.normed <- apply(dat,2,unitize_1)
  }else if(norm_tech==2){
    ## Mean_sd normalization
    dat.normed <- apply(dat,2,unitize_2)
  }else if(norm_tech==3){
    ## median-IQR normalization
    dat.normed <- apply(dat,2,unitize_3)
  }else{
    ## median-MAD normalization
    dat.normed <- apply(dat,2,unitize_4)
  }
  
 
  out1 <-  LocDensityOfProxiesAndOutliers(dat.normed,dat.outlier)   
  out2 <-  GraphFeaturesOfProxisAndOutliers(dat.normed,dat.outlier)
  feature.vector <- cbind(out1, out2)
  colnames(feature.vector)<- paste(c("OPO_" ), colnames(feature.vector), "_", norm_tech, sep="")
  return(feature.vector)
}



LocDensityOfProxiesAndOutliers <- function(dat, oo){
  require("FNN")
  dat <- as.data.frame(dat)
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  knn.nbrs <- knn.index(dat,pref.k)
  z.n <- knn.dist(dat,pref.k)
  z1 <- z.n[,pref.k]
  pot.outliers <- order(z1,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  outliers <- which(oo=='yes')
  
  cols1<- which(apply(dat,2,sd)!=0)
  dat1 <- dat[,cols1]
  dat1 <- as.data.frame(dat1)
  pca.obj <- prcomp(dat1)
  if(dim(dat1)[2]>10){
    cols = 1:10
  }else{
    cols = 1:dim(dat1)[2]
  }
  
  
  out <- data.frame(LocDen_KNOut_Mean=numeric(),LocDen_KNOut_Median=numeric(),LocDen_KNOut_SD=numeric(), LocDen_KNOut_IQR=numeric(), LocDen_KNOut_Max=numeric(), LocDen_KNOut_Min=numeric(), LocDen_KNOut_95P=numeric(), LocDen_KNOut_05P=numeric(), LocDen_LocDenOut_Mean=numeric(),LocDen_LocDenOut_Median=numeric(),LocDen_LocDenOut_SD=numeric(), LocDen_LocDenOut_IQR=numeric(), LocDen_LocDenOut_Max=numeric(), LocDen_LocDenOut_Min=numeric(), LocDen_LocDenOut_95P=numeric(), LocDen_LocDenOut_05P=numeric(), LocDenOut_KNOut_Common=numeric(), Out_KNOut_1=numeric(), Out_KNOut_2=numeric(), Out_LocDenOut_1=numeric(), Out_LocDenOut_2=numeric(), LocDen_Out_Mean=numeric(),LocDen_Out_Median=numeric(),LocDen_Out_SD=numeric(), LocDen_Out_IQR=numeric(), LocDen_Out_Max=numeric(), LocDen_Out_Min=numeric(), LocDen_Out_95P=numeric(), LocDen_Out_05P=numeric(), LocDenOut_Out_Mean=numeric(),LocDenOut_Out_Median=numeric(),LocDenOut_Out_SD=numeric(), LocDenOut_Out_IQR=numeric(), LocDenOut_Out_Max=numeric(), LocDenOut_Out_Min=numeric(), LocDenOut_Out_95P=numeric(), LocDenOut_Out_05P=numeric() )
  
  quant1 <- length(intersect(pot.outliers,outliers))/length(outliers)
  quant2 <- length(intersect(pot.outliers,outliers))/length(pot.outliers)
  
  
  if(length(cols)>1){
    col.pairs <- matrix(c(cols[-length(cols)], cols[-1]), ncol=2)
 
    abmatrix <- matrix(0,nrow=dim(col.pairs)[1],ncol=37)
    do.not.include.rows <- c()
    for(ii in 1:dim(col.pairs)[1]){
      sel.cols <- col.pairs[ii,]
      tryCatch(
        { ## Try part
          temp0 <- kde(pca.obj$x[,sel.cols],compute.cont=TRUE, eval.points = dat[,sel.cols])$estimate
          temp <- GetLocalDensity(temp0,knn.nbrs)
          abmatrix[ii,1:8] <- ComputeDensityFeatures(temp,pot.outliers)
          pot.outliers.density <- order(temp,decreasing=FALSE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
          abmatrix[ii,9:16] <- ComputeDensityFeatures(temp,pot.outliers.density)
          abmatrix[ii,17] <- length(intersect(pot.outliers,pot.outliers.density))/length(pot.outliers)
          
          ## Measures on outliers
          abmatrix[ii,18] <- quant1
          abmatrix[ii,19] <- quant2
          abmatrix[ii,20] <- length(intersect(pot.outliers.density,outliers))/length(outliers)
          abmatrix[ii,21] <- length(intersect(pot.outliers.density,outliers))/length(pot.outliers.density)
          abmatrix[ii,22:29] <- ComputeDensityFeatures(temp,outliers)
          
          abmatrix[ii,30:37] <- ComparePOWithOut(temp, outliers,1)
          
        },error=function(cond){
          abmatrix[ii,1:17] <- rep(1,17)
          do.not.include.rows <- c(do.not.include.rows,ii)
          message("Here's the original error message:")
          message(cond)
        }, warning=function(cond) {
          message("Here's the original warning message:")
          message(cond)
        }
      )
    }
    if(length(do.not.include.rows)>0){
      out[1,] <- apply(abmatrix[-do.not.include.rows,],2,mean)
    }else{
      out[1,] <- apply(abmatrix,2,mean)
    }
    
  }else{  ### only one column in dat
    temp0 <- kde(dat[,cols],eval.points = dat[,cols])$estimate
    temp <- GetLocalDensity(temp0,knn.nbrs)
    out1 <- ComputeDensityFeatures(temp,pot.outliers)
    pot.outliers.density <- order(temp,decreasing=FALSE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
    out2 <- ComputeDensityFeatures(temp,pot.outliers.density)
    out3 <- length(intersect(pot.outliers,pot.outliers.density))/length(pot.outliers)
    
    ## Measures on outliers
    quant3 <- length(intersect(pot.outliers.density,outliers))/length(outliers)
    quant4 <- length(intersect(pot.outliers.density,outliers))/length(pot.outliers.density)
    out4 <- ComputeDensityFeatures(temp,outliers)
    out5 <- ComparePOWithOut(temp, outliers,1)
    out[1,] <- c(out1,out2,out3, quant1, quant2,quant3,quant4, out4, out5  )
  }
  return(out)
}


GetLocalDensity <- function(dx,knn.nbrs){
  ## dx has the kde densities
  ## knn.nbrs has the knn neighbours
  out <- rep(0, length(dx))
  for(kk in 1:length(dx)){
    out[kk]<- ifelse(mean(dx[knn.nbrs[kk,]],na.rm = TRUE)==0, 0, dx[kk]/mean(dx[knn.nbrs[kk,]],na.rm = TRUE) )
  }
  return(out)
}

GraphFeaturesOfProxisAndOutliers <- function(dat,oo){
  require("igraph")
  require("FNN")
  dat <- as.data.frame(dat)
  outliers <- which(oo=='yes')
  
  out <- data.frame(GDegOut_Out_1=numeric(), GDegOut_Out_2=numeric(), GDeg_Out_Mean=numeric(), GDeg_Out_Median=numeric(), GDeg_Out_SD=numeric(), GDeg_Out_IQR=numeric(), GDeg_Out_Max=numeric(), GDeg_Out_Min=numeric(), GDeg_Out_Q95=numeric(), GDeg_Out_Q05=numeric(),  GDeg_PO_Mean=numeric(), GDeg_PO_Median=numeric(), GDeg_PO_SD=numeric(), GDeg_PO_IQR=numeric(), GDeg_PO_Max=numeric(), GDeg_PO_Min=numeric(), GDeg_PO_Q95=numeric(), GDeg_PO_Q05=numeric(),  GComp_Out_Mean=numeric(), GComp_Out_Median=numeric(), GComp_Out_SD=numeric(), GComp_Out_IQR=numeric(), GComp_Out_Max=numeric(), GComp_Out_Min=numeric(), GComp_Out_Q95=numeric(), GComp_Out_Q05=numeric(), GComp_PO_Mean=numeric(), GComp_PO_Median=numeric(), GComp_PO_SD=numeric(), GComp_PO_IQR=numeric(), GComp_PO_Max=numeric(), GComp_PO_Min=numeric(), GComp_PO_Q95=numeric(), GComp_PO_Q05=numeric(), GDist_O_Zero =numeric(), GDist_O_One = numeric(), GDist_O_Two = numeric(), GDist_O_Three = numeric(), GDist_O_Four = numeric(), GDist_O_Five = numeric(), GDist_O_AboveFive = numeric())
  
  pref.k <- min(ceiling(dim(dat)[1]/20),200)
  relations <- knn.index(dat,pref.k)
  dist <- knn.dist(dat,pref.k)
  vert <- 1:dim(dat)[1]
  g <- graph_from_data_frame(relations, directed=TRUE, vertices=vert)
  
  deg.vals <- degree(g)
  pot.outliers <- order(deg.vals,decreasing=TRUE)[1:min(ceiling(dim(dat)[1]*3/100),200)]
  quant1 <- length(intersect(pot.outliers,outliers))/length(outliers)
  quant2 <- length(intersect(pot.outliers,outliers))/length(pot.outliers)
  
  out1 <- RatiosBetweenTwoGroups(deg.vals,outliers)
  out2 <- RatiosBetweenTwoGroups(deg.vals,pot.outliers)
  
  
  
  comp <- components(g)
  comp.mem.o <-unique(comp$membership[outliers])
  comp.mem.po <-unique(comp$membership[pot.outliers])
  
  comp.csize <- comp$csize
  
  out3 <- RatiosBetweenTwoGroups(comp.csize,comp.mem.o)
  out4 <- RatiosBetweenTwoGroups(comp.csize,comp.mem.po)
  
  dist.from.out <- distances(g, v=V(g)[outliers], to=V(g)[-outliers], weights=NA)
  dist.from.out <- data.frame(dist.from.out)
  stats.d.f.o <- apply(dist.from.out,1,function(x)sum(is.finite(x)))
  out5 <- GetStatsFromGraphDistances(stats.d.f.o)
  
  out[1,] <- c(quant1, quant2, out1, out2, out3, out4, out5)
  return(out)
}


RatiosBetweenTwoGroups <-  function(dx, grp){
  ## dx is some 1D quantity like, degree, number of components, residuals, density etc
  ## grp is the group that it belongs to like outliers, pot.outliers
  
  F1 <- ifelse(mean(dx[-grp],na.rm = TRUE)==0, 0, mean(dx[grp],na.rm = TRUE)/mean(dx[-grp],na.rm = TRUE))
  F2 <- ifelse(median(dx[-grp],na.rm = TRUE)==0, 0, median(dx[grp],na.rm = TRUE)/median(dx[-grp],na.rm = TRUE))
  F3 <- ifelse(sd(dx[-grp],na.rm = TRUE)==0, 0, sd(dx[grp],na.rm = TRUE)/sd(dx[-grp],na.rm = TRUE))
  F4 <- ifelse(IQR(dx[-grp],na.rm = TRUE)==0, 0, IQR(dx[grp],na.rm = TRUE)/IQR(dx[-grp],na.rm = TRUE))
  F5 <- ifelse(max(dx[-grp],na.rm = TRUE)==0, 0, max(dx[grp],na.rm = TRUE)/max(dx[-grp],na.rm = TRUE))
  F6 <- ifelse(min(dx[-grp],na.rm = TRUE)==0, 0, min(dx[grp],na.rm = TRUE)/min(dx[-grp],na.rm = TRUE))
  F7 <- ifelse(quantile(dx[-grp],prob=0.95,na.rm = TRUE)==0, 0, quantile(dx[grp],prob=0.95,na.rm = TRUE)/quantile(dx[-grp],prob=0.95,na.rm = TRUE))
  F8 <- ifelse(quantile(dx[-grp],prob=0.05,na.rm = TRUE)==0, 0, quantile(dx[grp],prob=0.05,na.rm = TRUE)/quantile(dx[-grp],prob=0.05,na.rm = TRUE))
  out <- c(F1,F2,F3, F4, F5, F6, F7, F8)
  return(out)
}

GetStatsFromGraphDistances <- function(dx){
  ## dx has stats from dist.from.out from outliers to non-outliers
  F0 <- sum(dx==0)/length(dx) 
  F1 <- sum(dx==1)/length(dx) 
  F2 <- sum(dx==2)/length(dx)
  F3 <- sum(dx==3)/length(dx)
  F4 <- sum(dx==4)/length(dx)
  F5 <- sum(dx==5)/length(dx)
  F6 <- sum(dx>5)/length(dx)
  
  out <- c(F0,F1, F2, F3, F4, F5, F6)
  return(out)
}

##########################################################################################
##########################################################################################
#### FUNCTIONS FOR OUTLIERS SET 2 - LOCAL DENSITY AND IGRAPH FEATURES                #####
#### END                                                                             #####
##########################################################################################


#####################################################################################
#### Features that capture normalization
#### BEGIN
#####################################################################################
ComputeNormRelatedFeatures <- function(dat.o){
  dat.o <- data.frame(dat.o) 
  
  ## Pre-processing
  dat.2.mid <- dat.o[,colnames(dat.o)[colnames(dat.o)!="id"]]
  dat.1.mid <- dat.o[,colnames(dat.o)[(colnames(dat.o)!="id")& (colnames(dat.o)!="outlier")]]
  oo <- dat.o[,colnames(dat.o)[colnames(dat.o)=="outlier"]]
  
  dat.mid<- as.data.frame(dat.1.mid)
  
  
  dat <- data.frame(dat.mid)
  vecs <- ComputeWVectors(dat)
  features <- matrix(0, nrow=1, ncol=15)
  
  colnames(features) <- c("Angle_MM_IQR", "MAX_PROD_MM", "MIN_PROD_MM", "MEAN_PROD_MM", "MEDIAN_PROD_MM", "SD_PROD_MM", "IQR_PROD_MM", "PL1_PROD_MM", "MAX_PROD_IQR", "MIN_PROD_IQR", "MEAN_PROD_IQR", "MEDIAN_PROD_IQR", "SD_PROD_IQR", "IQR_PROD_IQR", "PL1_PROD_IQR" )
  if(dim(dat)[2]>1){
    vecs <- ComputeWVectors(dat)
    vec.mm <- vecs[,1]
    vec.iqr <- vecs[,2]   
    angle_1 <- sum(vec.mm*vec.iqr)
    
    out.feat <- ComputeAngleFeatures(dat, oo,vec.mm, vec.iqr, S=FALSE)
    norm.feat <- ComputeAngleFeatures(dat, oo,vec.mm, vec.iqr, S=TRUE)
    
    n <- dim(out.feat)[1]
    norm1 <-sort(norm.feat[,1], decreasing=TRUE)[1:n]
    norm2 <-sort(norm.feat[,2], decreasing=TRUE)[1:n]
    tt1 <- sort(out.feat[,1],decreasing=FALSE)/norm1
    tt2 <- sort(out.feat[,2],decreasing=FALSE)/norm2
    F1 <- GetBasicStats(tt1)
    F2 <- GetBasicStats(tt2)
    
    features[1,] <- c(angle_1,F1,F2 )
  }else{
    outl <- which(oo=="yes")
    tt<- sort(dat[outl,1])/sort(dat[-outl,1],decreasing=TRUE)[1:length(outl)]
    F1 <- GetBasicStats(tt)
    features[1,] <- c(1,F1,F1 )
    
  }
  
  return(features)
}


ComputeWVectors <- function(dat){
  dat <- data.frame(dat)
  temp.1 <- apply(dat,2, function(x)max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
  if(length(which(temp.1==0))>0){
    temp.1[which(temp.1==0)] <- 1
  }
  w1 <- normm(1/(temp.1^2))
  temp.2 <- apply(dat,2, function(x)IQR(x, na.rm=TRUE))
  if(length(which(temp.2==0))>0){
    temp.2[which(temp.2==0)] <- 1
  }
  w2 <- normm(1/(temp.2^2))
  W <- cbind(w1,w2)
  return(W)
}


ComputeAngleFeatures <- function(dat,oo,vec.mm, vec.iqr, S){
  #S = TRUE for normal data S is a switch
  require("FNN")
  if(S){
    recs <- which(oo!='yes')
  }else{
    recs <- which(oo=='yes')
  }
  kk <- min(ceiling(dim(dat)[1]/20), 100)
  kn.ind <-   knn.index(dat,kk)[,kk]
  x.vecs <- (dat[recs,] - dat[kn.ind[recs],])^2
  x.vecs.length <- apply(x.vecs,1,function(x) sum(x^2))
  x.vecs.n <- t(apply(x.vecs,1,normm))
  prod.n.mm <- apply(x.vecs.n,1,function(x)sum(x*vec.mm))
  prod.n.iqr <- apply(x.vecs.n,1,function(x)sum(x*vec.iqr))
  
  prod.1 <- apply(x.vecs,1,function(x)sum(x*vec.mm))
  prod.2 <- apply(x.vecs,1,function(x)sum(x*vec.iqr))
  prods <- cbind(prod.1,prod.2)
  colnames(prods) <- c("Min_Max_Vec", "IQR_Vec")
  return(prods)
}

GetBasicStats <- function(x){
  max.x <- max(x, na.rm=TRUE)
  min.x <- min(x, na.rm=TRUE)
  mean.x <- mean(x, na.rm=TRUE)
  median.x <- median(x, na.rm=TRUE)
  std.x <- sd(x, na.rm=TRUE)
  iqr.x <- IQR(x, na.rm=TRUE)
  pl1 <- sum(x<=1)/length(x)  # percentage less than 1
  return(c(max.x,min.x,mean.x,median.x,std.x,iqr.x,pl1))
}


normm <- function(x) {x / sqrt(sum(x^2))}

#####################################################################################
#### Features that capture normalization
#### END
#####################################################################################




