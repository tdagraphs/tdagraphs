#---------------------------------------------------------------------------------#
#TAD code: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks#
#---------------------------------------------------------------------------------#

# Given graph-files in edgelist format, this R code produces the corresponding    #
# Geodesic Distance Matrix (Geodesic Densification).                              #

####################################################################################################################################################################

rm(list = ls()) # Remove all the objects we created so far.

###################
#=== Libraries ===#
###################

library(MASS)
library(dodgr)
library(stringr)

setwd("/..../Ethereum/") # working directory #

##########################
#=== Input (Edgelist) ===#
##########################

pathEdgeListINDEX = 'FILES/ETH100_INDEX/'
#pathEdgeListFileAout = 'FILES/ETH100_AGD_EdgeList_A/' # Unweighted 
pathEdgeListFileWout = 'FILES/ETH100_AGD_EdgeList_W/'  # Weighted

####################################################################################################################################################################

lisFilesA <- list.files(pathEdgeListFileA)
lisFilesW <- list.files(pathEdgeListFileW)
NFiles <- length(lisFilesA) # We assume length(lisFilesA)=length(lisFilesW)
for(iii in 1:NFiles){
  #iii <- 1029 
  print(paste0("Files[",as.character(iii),"]: ", lisFilesA[iii], "  and  ",lisFilesW[iii]))
  
  if(file.info(paste0(pathEdgeListFileA,lisFilesA[iii]))$size==0 && file.info(paste0(pathEdgeListFileW,lisFilesW[iii]))$size==0){
    file.create(paste0(pathEdgeListFileAout,'AGD',lisFilesA[iii])) 
    file.create(paste0(pathEdgeListFileWout,'AGD',lisFilesW[iii])) 
  }else{
  
  #--- Read Original Edge Lists ---#
  matELA = as.matrix(read.table(paste0(pathEdgeListFileA,lisFilesA[iii]))) 
  matELW = as.matrix(read.table(paste0(pathEdgeListFileW,lisFilesW[iii]))) 

  #--- To find the number of vertices (Uunweighted) ---#
  auxStringA <- str_sub(lisFilesA[iii], 2, nchar(lisFilesA[iii])-4) # Letters and numbers
  #--- Delete last numbers ---#
  while(is.na(as.numeric( str_sub(auxStringA, nchar(auxStringA), nchar(auxStringA)) ))==FALSE){
    auxStringA <- str_sub(auxStringA, 1, nchar(auxStringA)-1)
  }
  #--- Open the INDEX file ---#
  matELINDEX = as.matrix(read.table(paste0(pathEdgeListINDEX,'INDEXnetwork',auxStringA,'XXX.txt'))) 
  NVerticesA <- matELINDEX[dim(matELINDEX)[1]]
  
  #--- To find the number of vertices (Weighted) ---#
  auxStringW <- str_sub(lisFilesW[iii], 2, nchar(lisFilesW[iii])-4) # Letters and numbers
  # Delete last numbers
  while(is.na(as.numeric( str_sub(auxStringW, nchar(auxStringW), nchar(auxStringW)) ))==FALSE){
    auxStringW <- str_sub(auxStringW, 1, nchar(auxStringW)-1)
  }
  #--- Open the INDEX file ---#
  matELINDEX = as.matrix(read.table(paste0(pathEdgeListINDEX,'INDEXnetwork',auxStringW,'XXX.txt'))) 
  NVerticesW <- matELINDEX[dim(matELINDEX)[1]]
  
  #--- Building the matrix ---#
  matAGD_A <- matrix(0, nrow=NVerticesA, ncol=NVerticesA) # For Weigths based on Distances  
  matAGD_W <- matrix(0, nrow=NVerticesW, ncol=NVerticesW) # For Weigths based on Distances  

  ###########################
  #=== GEODESIC DISTANCE ===#
  ###########################
  
  #--- Unweighted ---#  
  ddnetworkEdgeListA <- data.frame(matELA) 
  colnames(ddnetworkEdgeListA) <- c("from", "to", "d")
  matDisDD_A <- dodgr_dists(ddnetworkEdgeListA) # All distances
  # Set the distances (integers)
  vecColsA<-as.numeric(colnames(matDisDD_A))
  vecRowsA<-as.numeric(rownames(matDisDD_A))
  for(kc in 1:length(vecColsA)){
    for(kf in 1:length(vecRowsA)){
      if(is.na(matDisDD_A[kc, kf])==FALSE){ # To avoid NA
        matAGD_A[vecColsA[kc],vecRowsA[kf]] <- matDisDD_A[kc, kf]
      } # End if
    } # End for kf
  } # End for cf

  #--- Weighted ---#
  ddnetworkEdgeListW <- data.frame(matELW) 
  colnames(ddnetworkEdgeListW) <- c("from", "to", "d")
  matDisDD_W <- dodgr_dists(ddnetworkEdgeListW) # All distances
  # Set the distances (integers)
  vecColsW<-as.numeric(colnames(matDisDD_W))
  vecRowsW<-as.numeric(rownames(matDisDD_W))
  for(kc in 1:length(vecColsW)){
    for(kf in 1:length(vecRowsW)){
      if(is.na(matDisDD_W[kc, kf])==FALSE){ # To avoid NA
        matAGD_W[vecColsW[kc],vecRowsW[kf]] <- matDisDD_W[kc, kf]
      } # End if
    } # End for kf
  } # End for cf

  # Geodesic Distance (normalization - last step)
  matAGD_A <- matAGD_A/NVerticesA 
  matAGD_W <- matAGD_W/NVerticesW 
  
  # Save as EdgeList file  (Unweighted)
  AGD_A_EdgeList <- c()
  for(irow in 1:NVerticesA){
    for(icol in 1:NVerticesA){
      if(matAGD_A[irow,icol]>0){
        AGD_A_EdgeList <- rbind(AGD_A_EdgeList, c(irow, icol, matAGD_A[irow,icol]))
      }
    }
  }
  # Save as EdgeList file  (Weighted)
  AGD_W_EdgeList <- c()
  for(irow in 1:NVerticesW){
    for(icol in 1:NVerticesW){
      if(matAGD_W[irow,icol]>0){
        AGD_W_EdgeList <- rbind(AGD_W_EdgeList, c(irow, icol, matAGD_W[irow,icol]))
      }
    }
  }
  
  # Save EdgeLists
  write.table(AGD_A_EdgeList, file = paste0(pathEdgeListFileAout,'AGD',lisFilesA[iii]), sep = " ", row.names = F, col.names = F) 
  write.table(AGD_W_EdgeList, file = paste0(pathEdgeListFileWout,'AGD',lisFilesW[iii]), sep = " ", row.names = F, col.names = F) 
  }
}  

####################################################################################################################################################################
