#---------------------------------------------------------------------------------#
#TAD code: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks#
#---------------------------------------------------------------------------------#

# Using the information computed by ReducedDataset.m                              #    
# this Python's code produces a concise summary of each layer.                    #

####################################################################################################################################################################


###################
#=== Libraries ===#
###################

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os 
import glob 

####################
#=== Parameters ===#
####################
 
nameFolderINDEX = 'FILES/ETH100_INDEX/' 
nameFolderEdgeList = 'FILES/ETH100_RI/'  
nameFileOutput = 'FILES/ETH100summary.txt'

####################################################################################################################################################################


###############################################
#=== To extract the data from the files... ===#
###############################################

lisAllFiles = glob.glob(nameFolderINDEX+'/*.txt') 
lisAllFiles.sort() 
NFiles = len(lisAllFiles) 
summaryFile = open(nameFileOutput, 'w') # Create new file
summaryFile.close()
summaryFile = open(nameFileOutput, 'a') # Summary to append
for nf in range(NFiles):
    print("File["+str(nf+1)+"]... name: "+lisAllFiles[nf]+"\n") 
    splitFName = lisAllFiles[nf].split('/')
    networkName = splitFName[-1][12:-7]
    #print(nameFolderEdgeList+"RInetwork"+networkName+"XXX.txt\n")
    data = np.loadtxt(nameFolderEdgeList+"RInetwork"+networkName+"XXX.txt") # Load data
    # First and last date (timestamp)
    tstampStart = data[0][2] 
    tstampEnd = data[-1][2] 
    # First and last date (date)
    dtobjStart = datetime.fromtimestamp(tstampStart)
    dtobjEnd = datetime.fromtimestamp(tstampEnd)
    # First and last date (day:format)
    dayIni = datetime(dtobjStart.year,dtobjStart.month,dtobjStart.day,0,0,0)
    dayLast = datetime(dtobjEnd.year,dtobjEnd.month,dtobjEnd.day,0,0,0)
    dayLast = dayLast + timedelta(days=1) # To count the extra day
    # Save date  #
    summaryFile.write(str(int(tstampStart))+" "+str(int(tstampEnd))+" "+str(dayIni)[:-9]+" "+str(dayLast)[:-9]+" "+str(dayLast-dayIni)[:-14]+" "+networkName+"\n") 

summaryFile.close()

####################################################################################################################################################################

