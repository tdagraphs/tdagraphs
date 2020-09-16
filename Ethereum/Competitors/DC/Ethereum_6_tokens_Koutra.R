#----------------------------------------------------------------------------
#TAD: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks
#----------------------------------------------------------------------------

# This script loads all dynamic graphs....
# Database: 6 Ethereum Tokens 
# 

rm(list = ls()) # Remove all the objects we created so far.
library(igraph)
library(visNetwork)
library(networkD3)


library(astsa)
library(igraph)
library(network)
library(tseries)
library(sandwich)
library(forecast)
library(cointReg)
library(R.matlab)



# Path
setwd("/.../Ethereum")


source("DeltaCon.R")


script.name1 = "Ethereum_TDA_Koutra_ADJ_new6"
foo.ncol = 0
colnames1 = c("subj", "point")
write.table(t(colnames1), paste(script.name1, ".csv", sep = ""),
            append = F, sep = ",", row.names = F, col.names = F)

Index_days = c()

# ----- To load the summary... 
nameSummary = 'ETH100summary_new6.txt' 
tblSummary = read.table(nameSummary, header = FALSE, sep = "", dec = ".")
daysStart <- as.Date(tblSummary[,3])
daysEnd <- as.Date(tblSummary[,4])
NDays <- as.numeric(tblSummary[,5])
nameNetworks <- as.character(tblSummary[,6])

# **** To choose a type of graph ****
#pathEdgeListFile = 'ETH100_AGD_EdgeList_A_6tokens/AGDA' # Binary Graph based on the Average Geodesic Distance 
pathEdgeListFile = 'ETH100_AGD_EdgeList_W_6tokens/AGDW' # Weighted Graph based on the Average Geodesic Distance

maxNetworks <- 6; # Number of Tokens
nameFolderINDEX = 'ETH100_INDEX_new6/INDEXnetwork' 

pathEdgeSaveFile = "ETH100_AGD_EdgeList_A_6tokens/"

#Sim_matrix = matrix(0, ncol = maxNetworks, nrow = 450)



# ----- To Open every network 
for(kNet in 1:maxNetworks){
  #kNet <- 1
  NFiles <- NDays[kNet];
  # Look for the number of indices 
  matINDEX = as.matrix(read.table(paste0(nameFolderINDEX,nameNetworks[kNet],"XXX.txt"))) 
  NVertices <- matINDEX[dim(matINDEX)[1],1]
  
  
  
  trm = 0
  Sim_values = c()
  
  
  for(iii in 1:NFiles){
    #print(iii)
    nameFileNet <- paste0(pathEdgeListFile,nameNetworks[kNet],as.character(iii),'.txt')
    #print(nameFileNet)
    
    if(file.info(nameFileNet)$size==0){ # No edges were found in the network...
      next
    }else{
      #... matEdgeList contains: from, to, weight... each column...
      matEdgeList = as.matrix(read.table(nameFileNet)) 
      #print(matEdgeList)
      # To draw graph... 
      nodes <- data.frame("id"=1:NVertices, "label"=as.character(1:NVertices)) 
      edges <- data.frame("from"=matEdgeList[,1], "to"=matEdgeList[,2], "weight"=matEdgeList[,3]) 
      routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE) 
      # coords <- layout_in_circle(routes_igraph, V(routes_igraph)) 
      # V(routes_igraph)$shape <- "none" 
      
      net_A = routes_igraph
      edgesA = as.data.frame(get.edgelist(net_A, names = F))
      
      trm = trm + 1
      
      write.csv(edgesA, paste0(pathEdgeSaveFile,nameNetworks[kNet],"/","edges_",as.character(trm),".csv"))
    }
    
    
    
    
    
  }# End FOR daily network...
  
  Index_days[kNet] = trm
} 
write.csv(Index_days, paste0(pathEdgeSaveFile,"no_of_days",".csv"))


no_days = read.csv(paste0(pathEdgeSaveFile,"no_of_days",".csv"))


script.name1 = "Ethereum_TDA_Koutra_ADJ_new6"

alp = 1 - 0.95 # level of significance

# ----- To Open every network 
for(kNet in 1:maxNetworks){
  #kNet <- 1
  print(paste0("network is :",nameNetworks[kNet]))
  NFiles <- no_days[kNet,2];
  
  trm = 2
  Sim_values = c()
  
  for(iii in 1:(NFiles - 1)){
    print(iii)
    nameFileNet <- paste0(pathEdgeSaveFile,nameNetworks[kNet],"/","edges_",as.character(trm-1),".csv")
    edgesA = as.data.frame(read.csv(nameFileNet))[,-1]
    
    
    nameFileNet <- paste0(pathEdgeSaveFile,nameNetworks[kNet],"/","edges_",as.character(trm),".csv")
    edgesB = as.data.frame(read.csv(nameFileNet))[,-1]
    
    # Similiarity scores #
    Sim_values[iii] = delta_con(edgesA, edgesB, max(edgesA, edgesB), symmetrical = TRUE)
    
    
    trm  = trm + 1
    
    
    
  }# End FOR daily network...
  
  Sim_values = as.numeric(Sim_values)
  
  mSim.values = Sim_values[!is.na(Sim_values)]#-Remove all NA values
  len = length(mSim.values)
  row.size = len
  
  
  norm.test = shapiro.test(Sim_values) # normality test
  
  #==========================================================================#
  
  Mean.value = mean(Sim_values)  
  
  if(norm.test$p.value >0.05){Median.value = Mean.value}else{
    Median.value = median(Sim_values)}
  
  
  MRi = c()
  
  t = 0
  
  for(i in 2:row.size){
    
    t = t + 1
    MRi[t] = abs(Sim_values[i] - Sim_values[i-1])
    
  }
  
  MRconst = mean(MRi)
  
  conf.values_Lower = qnorm(alp/(2*maxNetworks),mean = Median.value, sd = MRconst, lower.tail = T)
  
  # Anomalous event detection 3
  ff_vals = which(Sim_values<=conf.values_Lower)
  
  
  foo = data.frame(subj = nameNetworks[kNet],
                   point = as.data.frame(t(ff_vals))
  )
  write.table(foo, paste(script.name1, "_",alp, ".csv", sep = ""), 
              append = T, sep = ",", row.names = F, col.names = F)
  foo.ncol = max(foo.ncol, ncol(foo))
  
  
  
}
