# This script loads all dynamic graphs....
# It would be very useful to use with competitors...
# Database: RIPPLE
# 

rm(list = ls()) # Remove all the objects we created so far.
# Path
setwd('/.../COMPETITORS/Zhao/RIPPLE')

library(igraph)
library(visNetwork)
library(networkD3)
library(gSeg)
library(ade4)
library(matrixcalc)
library(Rmisc)
library(Matrix)
source('Util.R')

# ----- To load the summary... 
nameSummary = 'RIP5Y90summary.txt' 
tblSummary = read.table(nameSummary, header = FALSE, sep = "", dec = ".")
daysStart <- as.Date(tblSummary[,3])
daysEnd <- as.Date(tblSummary[,4])
NDays <- as.numeric(tblSummary[,5])
nameNetworks <- as.character(tblSummary[,6])
lisSelTokens <- c(1, 2, 3, 4, 5)

# **** To choose a type of graph ****
#pathEdgeListFile = 'RIPPLE_EdgeList_A/A' # Graph based on Adjacency (only 0 and 1)
#pathEdgeListFile = 'RIPPLE_EdgeList_W/W' # Graph based on Weighted Adjacency (values between 0 and 1)
#pathEdgeListFile = 'RIPPLE_EdgeList_A_GD/GDA' # Graph based on the Geodesic Distance - Graph based on Adjacency
pathEdgeListFile = 'RIPPLE_EdgeList_W_GD/GDW' # Graph based on the Geodesic Distance - Graph based on Weight 

maxNetworks <- 5; # Number of Tokens 
nameFolderINDEX = 'RIPPLE_INDEX/INDEXnetwork' 

# ----- To Open every network 
for(kNet in lisSelTokens){
  print(paste0('********* Working on: ', nameNetworks[kNet]))
  #kNet <- 1
  NFiles <- NDays[kNet];
  # Look for the number of indices 
  matINDEX = as.matrix(read.table(paste0(nameFolderINDEX,nameNetworks[kNet],"XXX.txt"))) 
  NVertices <- matINDEX[dim(matINDEX)[1],1]
  
  # For method (start) -------------
  total_Time <- NFiles
  Time_sets <- 1:total_Time
  total_node <- NVertices
  node_sets <- 1:total_node
  adj_matrix_complete_daily <- list() 
  # For method (end) -------------
  
  for(time_index in 1:NFiles){
    nameFileNet <- paste0(pathEdgeListFile,nameNetworks[kNet],as.character(time_index),'.txt')
    #print(nameFileNet)
    
    if(file.info(nameFileNet)$size==0){ # No edges were found in the network...
      matEL = c()
    }else{
      #... matEdgeList contains: from, to, weight... each column...
      matEL = as.matrix(read.table(nameFileNet)) 
    }
    
    # HERE!!! You must use the graph loaded in 'matEdgeList'...
    # For method (start) -------------
    # Build Symmetric matrix (old school)
    tmp <- matrix(0, total_node, total_node) 
    colnames(tmp) <- rownames(tmp) <- node_sets
    if(!is.null(dim(matEL)[1])){ # No empty file...
      for (i in 1:(dim(matEL)[1])){
        tmp[matEL[i,1], matEL[i,2]] <- matEL[i,3]
        tmp[matEL[i,2], matEL[i,1]] <- matEL[i,3]
      }
    }
    adj_matrix_complete_daily[[time_index]] <- tmp
    # For method (end) -------------
  }# End FOR daily network...
  plot(sapply(adj_matrix_complete_daily, sum), main=nameNetworks[kNet])
  
  ## Network CP Detection
  adj_matrix_complete <- adj_matrix_complete_daily
  Time <- length(adj_matrix_complete)
  n <- total_node
  h <- floor(sqrt(Time)) # sqrt(T)
  #h <- 14 # 2 weeks
  # h <- 7 # 1 week
  #h <- 4 # days
  #h <- 3 # days 
  #h <- floor(sqrt(Time/3)) # sqrt(T)
  
  ## ChenZhang(2015) via Binary Segmentation
  E <- adj_matrix_complete
  E <- t(sapply(E, as.vector))
  dist_frobenius <- dist(E)
  similarity_measure <- mstree(dist_frobenius)
  #r1 <- gseg1(Time, similarity_measure, n0=h, n1=Time-h+1, statistics="original")
  cz_result <- binaryS(E, h)
  if(!is.na(cz_result)){
    cz_result <- sort(cz_result[!is.na(cz_result)])
  }
  
  ## MNBS
  T0 <- h
  T1 <- Time-h
  E <- adj_matrix_complete
  d2_MNBS <- rep(0, Time) # store the SaRa statistics
  beta <- 0.5 # recommendated tuning parameters
  delta <- 0.1
  C0 <- 3
  D0 <- 0.25
  threshold <- D0*(log(n))^(1/2+delta)/sqrt(n)/h^beta # thresholding rule
  for (t in T0:T1){
    E1 <- E[(t-h+1):t]
    E2 <- E[(t+1):(t+h)]
    
    ## MNBS 
    p_MNBS_1 <- mnbs(beta,E1,C=C0)$P 
    p_MNBS_2 <- mnbs(beta,E2,C=C0)$P
    d2_MNBS[t] <- (d2infty_norm(p_MNBS_1-p_MNBS_2))^2 # d2infty norm
  }
  MNBS_result <- SaRa(d2_MNBS, h, threshold)
  MNBS_lm <- local_maximizer(d2_MNBS, h)
  
  ## Plot the result
  # png(paste0('MITnetwork.png'), width=10, height=10, units='in', res=300)
  par(mfrow=c(2,1))
  plot(Time_sets, d2_MNBS, type='l', main='Estimated change-points by MNBS',
       ylab=expression(paste(D(t,h))), xlab='', xaxt="n")
  #axis(1, Time_sets[seq(1,Time,length.out=12)], format(Time_sets[seq(1,Time,length.out=12)], "%m/%d/%Y"), cex.axis = .7)
  abline(h=threshold)
  points(Time_sets[MNBS_lm], d2_MNBS[MNBS_lm], col='red')
  abline(v=Time_sets[MNBS_result], col='blue', lty=2)
  legend_names <- c('local maximizer', 'change-point', 'threshold')
  legend('topright', legend_names, lty=c(NA, 2, 1), pch=c(1, NA, NA), col=c('red', 'blue', 'black'))
  
  plot(Time_sets, sapply(adj_matrix_complete_daily, sum)/2, ylab='Total links', xlab='', xaxt="n",
       main='Total daily links of the MIT dynamic network')
  #axis(1, Time_sets[seq(1,Time,length.out=12)], format(Time_sets[seq(1,Time,length.out=12)], "%m/%d/%Y"), cex.axis = .7)
  abline(v=Time_sets[MNBS_result], col='blue', lty=2)
  legend('topright', 'change-point', lty=2, col=c('blue'))
  # abline(v=Time_sets[cz_result], col='red', lty=2)
  # dev.off()
  
  print(cz_result)
  print(MNBS_result)
  
  write.table(cz_result, file = paste0('IDXTimeResul/',unlist(strsplit(pathEdgeListFile,'/'))[1],'_',nameNetworks[kNet],'_CZ.txt'), sep = " ", row.names = F, col.names = F, quote=F) 
  write.table(MNBS_result, file = paste0('IDXTimeResul/',unlist(strsplit(pathEdgeListFile,'/'))[1],'_',nameNetworks[kNet],'_MNBS.txt'), sep = " ", row.names = F, col.names = F, quote=F) 
  
}# End FOR TOKEN




