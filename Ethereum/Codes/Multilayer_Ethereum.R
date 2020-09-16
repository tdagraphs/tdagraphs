#---------------------------------------------------------------------------------#
#TAD code: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks#
#---------------------------------------------------------------------------------#

# This R code compute the anomaly detection on multilayer dynamic networks        # 
# using the files generated in previous steps.                                    #
# See the submitted paper for thorough procedure description. (TAD)               #

####################################################################################################################################################################

rm(list = ls()) # Remove all the objects we created so far

###################
#=== LIBRARIES ===#
###################

library(R.matlab)
library(TDA)
library(caret)
library(e1071)

#--- To use an approach from the Twitter company (2015) ---#
#install.packages("devtools")
#devtools::install_github("twitter/AnomalyDetection")

library(AnomalyDetection)


#=== Working Directory ===#
setwd("/...PATH_ETHEREUM.../")

####################################################################################################################################################################

#####################
#=== Input Files ===#
#####################
#pathPDs = 'FILES/ETH100_AGD_A_PD/AGDA' # Graph based on the Average Geodesic Distance - Unweighted
pathPDs = 'FILES/ETH100_AGD_W_PD/AGDW' # Graph based on the Average Geodesic Distance - Weighted

#############################
#=== Necessary summaries ===#
############################# 
nameSummary = 'FILES/ETH100summary.txt' 
tblSummary = read.table(nameSummary, header = FALSE, sep = "", dec = ".")
daysStart <- as.Date(tblSummary[,3])
daysEnd <- as.Date(tblSummary[,4])
NDays <- as.numeric(tblSummary[,5])
nameNetworks <- as.character(tblSummary[,6])
maxNetworks <- 31; # Number of Tokens

#####################
#=== Output file ===#
#####################
splitPathPD <- unlist(strsplit(pathPDs, "/"))
fileOutputTable <- paste0('TableOutliers/',splitPathPD[2],'_OD2014.txt') # OD2014  or OD1993

############################################################################
#=== The strongest connections (6 tokens with the highest connectivity) ===#
############################################################################
lisSC <- c()
lisSC <- c(lisSC, 'decentraland')
lisSC <- c(lisSC, 'tierion')
lisSC <- c(lisSC, 'vechain') 
lisSC <- c(lisSC, 'zrx') 
lisSC <- c(lisSC, 'cybermiles') 
lisSC <- c(lisSC, 'bytom') 
NCToken <- length(lisSC); # Number of Tokens
#lisSC[sample(1:NCToken)] # Random order > it produces the same...

##############################################################################################################################################

###################################
#=== Search indices in summary ===#
###################################

indFindTK <- c()
for (kTok in 1:NCToken) {
  indFindTK[kTok] <- match(lisSC[kTok], nameNetworks)
}

############################
#=== Search time window ===#
############################

limLowDay <- max(daysStart[indFindTK]) 
limUppDay <- min(daysEnd[indFindTK]) 

##############################################################################################################################################

##############################
#=== Distance between PDs ===#
##############################

# 1. Load PD's in subsequent graphs and measure #
difDays <- as.numeric(limUppDay-limLowDay)-1
vecKNet <- 1+as.numeric(rep(limLowDay,NCToken)-daysStart[indFindTK])
tableOutliers2014 <- c()
wasser_dist_012n <- c()
for (sDay in 0:difDays) { # For each day in the time winddow
  # To build the Big Persistent diagram
  #sDay <- 0
  PDBig <- c()
  for (tk in 1:NCToken) {
    nameFileNet <- paste0(pathPDs,nameNetworks[indFindTK[tk]],as.character(vecKNet[tk]+sDay),'_PD.txt')
    P = as.matrix(read.table(nameFileNet)) 
    P[P[,3] == -1, 3] = 1 # Replace -1=Infinity by...  1 (which is the maximum possible value)
    P[,1] <- tk-1
    PDBig <- rbind(PDBig, P)
  }
  # 2. Computing distances (Pn-1 vs Pn) #
  if(sDay==0){ PDBig_ant = PDBig }
  wasser_dist_012n = c(wasser_dist_012n, wasserstein(PDBig_ant, PDBig, dimension=seq(0,(NCToken-1)) ))
  PDBig_ant = PDBig
  
} # End FOR

#=== Normalize Distances ===#
WTDA = wasser_dist_012n
NORM_WTDA = wasser_dist_012n/max(wasser_dist_012n)

##############################################################################################################################################

###################################
#=== Anomaly detection (S-ESD) ===#
###################################

# alpha=0.05
resOutlier = AnomalyDetectionVec(data.frame(WTDA), max_anoms=0.1, period=60, direction='both', only_last=FALSE, plot=TRUE) 

#=== Days in DATE format ===#
xDays <- seq(limLowDay, limUppDay, by="days")
xDays <- xDays[1:length(xDays)-1]

#=== Add Outliers ===# 
if(dim(resOutlier$anoms)[1]>0){
  # Building outlier's table... 
  tableOutliers2014 <- as.character(xDays[resOutlier$anoms[,1]])
}

#=== Save tables ( Outlier method: 2014 ) ===#
write.table(tableOutliers2014, file=fileOutputTable, append = FALSE, quote=FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)


##############################################################################################################################################

###############################################
#=== Match anomalies with the Ground Truth ===#
###############################################

pathGroundTruth <- 'gt.txt' # 
NOutliers<- length(tableOutliers2014)
dayOutlier <- tableOutliers2014 

#=== Open Ground Truth (Old School) ===#
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

#=== Events inside the time window ===#
flagEventWindow <- as.Date(dayGT)>=limLowDay & as.Date(dayGT)<limUppDay 

#=== Search Outliers in the Ground Truth ===#
foundEvent <- vector(mode="character", length=NOutliers)
NODetected <- 0
indexGTFound <- c() 
for (i in 1:NOutliers) {
  for (k in 1:NEvents) {
    subsDays <- as.numeric(as.Date(dayGT[k])-as.Date(dayOutlier[i]))
    if(abs(subsDays)<=2){ # Without direction
    #if(subsDays==0){ # The very same day....
    #if(subsDays>0 & subsDays<=2){ # Event AFTER the outlier...
    #if(subsDays>=-2 & subsDays<0){ # Event BEFORE the outlier...
      foundEvent[i] <- paste0(foundEvent[i], dayGT[k], ' ', eventGT[k],'* ')
      NODetected <- NODetected + 1
      indexGTFound <- c(indexGTFound, k)
    }
  }# En FOR K
}# End FOR I

##################################
#=== Save Table with outliers ===#
##################################

tableDetectedOutliers<-cbind(dayOutlier, foundEvent)
colnames(tableDetectedOutliers)<-c("Date_Anomaly", "Outliers_Detected") 
pathDetected <- paste0(unlist(strsplit(fileOutputTable,".txt")),'_Detected.txt')
write.table(tableDetectedOutliers, file=pathDetected, append = FALSE, quote=FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)

##############################################################################################################################################

##############################
#=== Graphs and Summaries ===#
##############################
#############
# 1. Graphs #
#############
vecGroundTruth <- as.Date(dayGT[flagEventWindow])

#=== Add Outliers (Method 2014) ===#
plot(xDays,WTDA, type = "l", col="blue", lwd=2, main = paste0("Multilayer Network"), xlab = "Days", ylab = expression(paste(W[](D[P-1],D[P]))) )
for(i in 1:length(vecGroundTruth)){
  abline(v=vecGroundTruth[i], col="green")
}# End For

matplot(xDays,WTDA, type = "l", col="blue", lwd=2, main = paste0("Multilayer Network"), xlab = "Days", ylab = expression(paste(W[](D[P-1],D[P]))) , add=T)

#=== Add Outliers  ===#
if(dim(resOutlier$anoms)[1]>0){
  # Adding points....
  for(i in 1:length(resOutlier$anoms[,1])){
    points(xDays[resOutlier$anoms[i,1]], WTDA[resOutlier$anoms[i,1]], pch = 18, col = "red", bg = "yellow", cex = 1.5)
  }
  # Building outlier's table # 
  #tableOutliers2014 <- as.character(xDays[resOutlier$anoms[,1]])
}

################
# 2. Summaries #
################
print(paste0('*Time Window: ',limLowDay,'  ',limUppDay))
print(paste0('*Total days: ',length(xDays)))
print(paste0('*True Events: ',length(vecGroundTruth)))
print(paste0('*Outliers obtained: ',NOutliers))
print(paste0('*Events Detected: ',NODetected))
print(paste0('*Read file: ',pathDetected))
print('---')
cmTP <- NODetected
cmFP <- NOutliers-NODetected
cmTN <- (length(xDays)-length(vecGroundTruth))-(NOutliers-NODetected)
cmFN <- length(vecGroundTruth)-NODetected
print(paste0('*True Positive: ',cmTP))
print(paste0('*False Positive: ',cmFP))
print(paste0('*True Negative: ',cmTN))
print(paste0('*False Negative: ',cmFN)) 


#=== To create binary vectors ===#
vecTruth <- rep(0, times = length(xDays))
vecTruth[1:length(vecGroundTruth)] <- 1 # Assign 1 on Ground truth years  
vecPred <- rep(0, times = length(xDays))
vecPred[(cmFN+1):(cmFN+NOutliers)] <- 1 # Assign 1 on predicted years  

#=== To compute statistcs ===#
tablePT <- table(factor(vecPred,levels=c(1,0)), factor(vecTruth,levels=c(1,0)))

#=== Confusion matrix (compute) ===#
resulCM <- confusionMatrix(tablePT)

#=== Matrix with statistics ===#
NMethods <- 1
NMetrics <- 24
matMetrics <- matrix(-1, nrow = NMetrics, ncol = NMethods)

#=== Fill out the statistics ===#
iM <- 1
matMetrics[1, iM] <- length(vecGroundTruth)
matMetrics[2, iM] <- NOutliers
matMetrics[3, iM] <- resulCM$table[1,1]
matMetrics[4, iM] <- resulCM$table[2,2]
matMetrics[5, iM] <- resulCM$table[1,2]
matMetrics[6, iM] <- resulCM$table[2,1]
matMetrics[7, iM] <- resulCM$overall[1]
matMetrics[8, iM] <- resulCM$overall[2]
matMetrics[9, iM] <- resulCM$overall[3]
matMetrics[10, iM] <- resulCM$overall[4]
matMetrics[11, iM] <- resulCM$overall[5]
matMetrics[12, iM] <- resulCM$overall[6]
matMetrics[13, iM] <- resulCM$overall[7]
matMetrics[14, iM] <- resulCM$byClass[1]
matMetrics[15, iM] <- resulCM$byClass[2]
matMetrics[16, iM] <- resulCM$byClass[3]
matMetrics[17, iM] <- resulCM$byClass[4]
matMetrics[18, iM] <- resulCM$byClass[5]
matMetrics[19, iM] <- resulCM$byClass[6]
matMetrics[20, iM] <- resulCM$byClass[7]
matMetrics[21, iM] <- resulCM$byClass[8]
matMetrics[22, iM] <- resulCM$byClass[9]
matMetrics[23, iM] <- resulCM$byClass[10]
matMetrics[24, iM] <- resulCM$byClass[11]

#=== Matrix to save performace metrics ===#
charMetrics <- vector(mode = "character", length = NMetrics)
charMetrics[1] <- "True Anomalies"
charMetrics[2] <- "# Anomaly Found"
charMetrics[3] <- "True Positive"
charMetrics[4] <- "True Negative"
charMetrics[5] <- "False Positive"
charMetrics[6] <- "False Negative"
charMetrics[7] <- "Accuracy"
charMetrics[8] <- "Kappa"
charMetrics[9] <- "Accuracy Lower"
charMetrics[10] <- "Accuracy Upper"
charMetrics[11] <- "Accuracy Null"
charMetrics[12] <- "Accuracy P-Value"
charMetrics[13] <- "Mcnemar P-Value"
charMetrics[14] <- "Sensitivity"
charMetrics[15] <- "Specificity"
charMetrics[16] <- "Pos Pred Value"
charMetrics[17] <- "Neg Pred Value"
charMetrics[18] <- "Precision"
charMetrics[19] <- "Recall"
charMetrics[20] <- "F1"
charMetrics[21] <- "Prevalence"
charMetrics[22] <- "Detection Rate"
charMetrics[23] <- "Detection Prevalence"
charMetrics[24] <- "Balanced Accuracy"

#=== Save metrics ===#
tblOut <- matMetrics
rownames(tblOut) <- charMetrics
colnames(tblOut) <- c('Multiplex')
# Save table Metrics
fileOutMetrics <- paste0('Results/',splitPathPD[2],'_ALLMetrics.csv') # Metrics 
write.csv(tblOut,fileOutMetrics)

#=== Extra (Save the found days) ===#
listFoundEventsTXT <- vector(mode="character", length=length(xDays))
indexXDaysFound <- match(as.Date(dayGT[indexGTFound]), xDays) # Match Days from ground truth with real days
listFoundEventsTXT[indexXDaysFound] <- eventGT[indexGTFound] 

##############################################################################################################################################

#####################
#=== Final table ===#
#####################

tblAllFoundEvents <- cbind(as.character(xDays),WTDA,listFoundEventsTXT) 
#tblAllFoundEvents <- cbind(WTDA,listFoundEventsTXT) 
#rownames(tblAllFoundEvents) <- as.character(xDays)
colnames(tblAllFoundEvents) <- c("Days","Wasserstein","EventsFound")
fileOutFoundEvents <- paste0('Results/',splitPathPD[2],'_FoundEvents.csv') # Found Events
write.csv(tblAllFoundEvents, fileOutFoundEvents)

##########################
#=== To open the file ===#
##########################

dataAD <- read.table(fileOutFoundEvents,sep=",",fill=TRUE,header = TRUE)
# Indices where an event was found
indexFoundEvent <- dataAD$X[dataAD$EventsFound!=""]
fileOutFoundEventsShort <- paste0('Results/',splitPathPD[2],'_FoundEvents_Short.csv') # Found Events
write.csv(dataAD[indexFoundEvent,c(2,4)], fileOutFoundEventsShort)
  
  
