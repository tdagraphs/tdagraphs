#-------------------------------------------------------------------------------#
#S-TAD: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks #
#-------------------------------------------------------------------------------#

# This R code compute the anomaly detection on unilayer dynamic networks        # 
# using the files generated in previous steps.                                  #
# See the submitted paper for complete description (S-TAD)                      #

##############################################################################################################################################

rm(list = ls()) # Remove all the objects we created so far.

###################
#=== Libraries ===#
###################

library(R.matlab)
library(TDA)
#=== To use approach from the Twitter company (2015) [S-ESD] ===#
#install.packages("devtools")
#devtools::install_github("twitter/AnomalyDetection")

library(AnomalyDetection)

##############################################################################################################################################

setwd("/...PATH_RIPPLE.../") #=== Working Directory ===#

#####################
#=== Input Files ===#
#####################

nameSummary = 'RIP5Y90summary.txt'
tblSummary = read.table(nameSummary, header = FALSE, sep = "", dec = ".")
daysStart <- as.Date(tblSummary[,3])
daysEnd <- as.Date(tblSummary[,4])
NDays <- as.numeric(tblSummary[,5])
nameNetworks <- as.character(tblSummary[,6])

# ****** Methods (NORMAL) --- To choose One of them ******** # K-Clique-Commmunity (C++) 

#=== Unweighted ===#
pathPDs = 'RIPPLE_PD_A_GD/GDA'; # GD on Adjacency (PD)
nameFileIDX <- 'IDXTimeResul/RIPPLE_AGD_EdgeList_A_'

#=== Weighted ===#
#pathPDs = 'RIPPLE_PD_W_GD/GDW' # GD on Weighted-Adjacency (PD)
#nameFileIDX <- 'IDXTimeResul/RIPPLE_AGD_EdgeList_W_' 

#=== Output ===#
splitPathPD <- unlist(strsplit(pathPDs, "/"))
fileOutputTable <- paste0('TableOutliers/',splitPathPD[1],'_OD2014.txt') # OD2014  or OD1993

##############################################################################################################################################

##############################
#=== Distance between PDs ===#
##############################

maxNetworks <- 5; 
tableOutliers2014 <- c() 
# WTDA, N_Outliers, dateOutlier, indexOutlier, NEvents_found, dateEveFound, indexEveFound
# N_Ground_Truth, Dates_GT, index_GT
allKCCPTDA <- matrix(list(), maxNetworks, 7) 
dateLow <- as.Date("2016-10-09") 
dateUp <- as.Date("2020-03-24") 
lisSelTokens <- c(1, 2, 3, 4, 5) 
for(kNet in lisSelTokens){
  #kNet <- 1 # Index of Network

  #==== For every file (PD) ===#
  wasser_dist_01 =  c()
  NFiles <- NDays[kNet];
  for(iii in 1:NFiles){
    #iii <- 2
    #print(as.character(iii))
    
    nameFileNet <- paste0(pathPDs,nameNetworks[kNet],as.character(iii),'_PD.txt') 
    P = as.matrix(read.table(nameFileNet)) 
    P[P[,3] == -1, 3] = 1 # Replace -1=Infinity by...  1 (which is the maximum possible value)
    
    #=== Computing distance (Pn-1 vs Pn) ===#
    if (iii == 1){ P_ant = P }
    wasser_dist_01 = c(wasser_dist_01, wasserstein(P_ant, P, dimension = c(0,1)))
    P_ant = P
    
  } # End For
  
  WTDA = wasser_dist_01
  NORM_WTDA = wasser_dist_01/max(wasser_dist_01) #=== Normalize distances ===# 

######################################################################################################################################################

###################################
#=== Outlier Detection (S-ESD) ===#
###################################

  # 1: AnomalyDetection(Twitter Co, 2015) 
  # A period of 2 months and only taking at most 5% of anomalies
  #resOutlier = AnomalyDetectionVec(data.frame(WTDA), max_anoms=0.05, period=60, direction='both', only_last=FALSE, plot=TRUE) 
  # A period of 2 month and only taking at most 10% of anomalies
  #resOutlier = AnomalyDetectionVec(data.frame(WTDA), max_anoms=0.1, period=60, direction='both', only_last=FALSE, plot=TRUE) 
  
  #=== The new ones ===#
  resOutlier = AnomalyDetectionVec(data.frame(WTDA), max_anoms=0.1, period=480, alpha=(0.05/5), direction='both', only_last=FALSE, plot=TRUE) 
  
   #=== Plots and plot summaries ===#
  # Days in DATE format
  xDays <- seq(daysStart[kNet], daysEnd[kNet], by="days")
  xDays <- xDays[1:length(xDays)-1]

  #=== Add Outliers ===#
  allKCCPTDA[[kNet,1]] <- WTDA # Time serie
  allKCCPTDA[[kNet,2]] <- dim(resOutlier$anoms)[1] # Number of anomalies detected
  if(dim(resOutlier$anoms)[1]>0){

    #=== Building outlier's table ===#
    tableOutliers2014 <- rbind(tableOutliers2014, cbind(rep(nameNetworks[kNet], length(resOutlier$anoms[,1])), as.character(xDays[resOutlier$anoms[,1]])))
   
    #=== Adding ===#
    allKCCPTDA[[kNet,3]] <- xDays[resOutlier$anoms[,1]] # Days detected
    allKCCPTDA[[kNet,4]] <- resOutlier$anoms[,1] # Indices detected
  }
} # End FOR


#####################################
#=== Save tables ( S-ESD: 2014 ) ===#
#####################################

colnames(tableOutliers2014)<-c("Token", "DateAnomaly")
write.table(tableOutliers2014, file=fileOutputTable, append = FALSE, quote=FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)


######################################################################################################################################################

pathGroundTruth <- 'gt.txt' # 

#######################
#=== Open Outliers ===#
####################### 

matOutliers <- tableOutliers2014  
NOutliers<- dim(matOutliers)[1]  
nameToken <- as.character(matOutliers[,1])  
dayOutlier <- as.character(matOutliers[,2])  

########################################
#=== Open Ground Truth (Old School) ===#
########################################

linesGT <- readLines(pathGroundTruth)
NEvents <- length(linesGT)
dayGT <- c() # To save dates
eventGT <- c() # To save events
for(i in 1:NEvents){
  #i<-1
  lineGTsplit <- unlist(strsplit(linesGT[i], split = " "))
  lineGTsplit <- lineGTsplit[lineGTsplit!=""] # Remove all empty elements
  NSplit <- length(lineGTsplit)
  # Build date
  dayGT <- rbind(dayGT, as.character(as.Date(lineGTsplit[1], format='%d-%m-%Y')))
  # Build comments
  eventGT <- rbind(eventGT, paste(lineGTsplit[2:NSplit], sep = '', collapse = ' '))
}

#############################################
#=== Search Outliers in the Ground Truth ===#
#############################################

foundEvent <- vector(mode="character", length=NOutliers)
NODetected <- 0
for (i in 1:NOutliers) {
  for (k in 1:NEvents) {
    subsDays <- as.numeric(as.Date(dayGT[k])-as.Date(dayOutlier[i]))
    if(abs(subsDays)<=2){ # Without direction
      #if(subsDays==0){ # The very same day....
      #if(subsDays>0 & subsDays<=2){ # Event AFTER the outlier...
      #if(subsDays>=-2 & subsDays<0){ # Event BEFORE the outlier...
      foundEvent[i] <- paste0(foundEvent[i], dayGT[k], ' ', eventGT[k],'* ')
      NODetected <- NODetected + 1
    }
  }# En FOR K
}# End FOR I

##################################
#=== Save Table with outliers ===#
##################################

tableDetectedOutliers<-cbind(nameToken, dayOutlier, foundEvent)
colnames(tableDetectedOutliers)<-c("Token", "Date_Anomaly", "Outliers_Detected") 
pathDetected <- paste0(unlist(strsplit(fileOutputTable,".txt")),'_Detected.txt')
write.table(tableDetectedOutliers, file=pathDetected, append = FALSE, quote=FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)


#########################
#=== Some statistics ===#
#########################

uniqueNetworks<-unique(nameToken)
NUniqueNet <- length(uniqueNetworks)
print("Network / Number Outliers found / Number Match Ground Truth")
for (i in 1:NUniqueNet) {
  tblDetected <- tableDetectedOutliers[tableDetectedOutliers[,1]==uniqueNetworks[i],]
  if(!is.null(dim(tblDetected)[1])){
    print(paste(uniqueNetworks[i], as.character(dim(tblDetected)[1]), as.character(sum(tblDetected[,3]!=""))))
  }
}
print(paste(as.character(sum(tableDetectedOutliers[,3]!="")), "out of", as.character(dim(tableDetectedOutliers)[1]))) 


######################################################################################################################################################

##############################
#=== Graphs and Summaries ===#
##############################

for(kNet in lisSelTokens){
  #kNet <- 1
  # Events inside the time window 
  flagEventWindow <- as.Date(dayGT)>=daysStart[kNet] & as.Date(dayGT)<daysEnd[kNet] 
  # Ground Truth 
  vecGroundTruth <- as.Date(dayGT[flagEventWindow])
  # xDays 
  xDays <- seq(daysStart[kNet], daysEnd[kNet], by="days")
  xDays <- xDays[1:length(xDays)-1]
  # Method 1... 
  WTDA <- allKCCPTDA[[kNet,1]]
  plot(xDays,WTDA, type = "l", col="blue", lwd=2, main = paste0(nameNetworks[kNet]," - Outlier Detection 2014"), xlab = "Week", ylab = expression(paste(W[2](D[P-1],D[P]))) )
  for(i in 1:length(vecGroundTruth)){
    abline(v=vecGroundTruth[i], col="green")
  }# End For
  matplot(xDays,WTDA, type = "l", col="blue", lwd=2, main = paste0(nameNetworks[kNet]," - Outlier Detection 2014"), xlab = "Week", ylab = expression(paste(W[2](D[P-1],D[P]))) , add=T)
  # Add Outliers 
  if(dim(resOutlier$anoms)[1]>0){
    for(i in 1:length(resOutlier$anoms[,1])){
      points(xDays[resOutlier$anoms[i,1]], WTDA[resOutlier$anoms[i,1]], pch = 18, col = "red", bg = "yellow", cex = 1.5)
    }
  }
  
   #=== Save information ===#
  allKCCPTDA[[kNet,5]] <- length(vecGroundTruth)
  allKCCPTDA[[kNet,6]] <- vecGroundTruth
  allKCCPTDA[[kNet,7]] <- match(vecGroundTruth, xDays)
  
  #=== To save outliers indices ===#
  indicesOutliers <- allKCCPTDA[[kNet,4]]
  write.table(indicesOutliers, file = paste0(nameFileIDX,nameNetworks[kNet],'_KCCPTDA.txt'), sep = " ", row.names = F, col.names = F, quote=F) 
}


######################################################################################################################################################


