#---------------------------------------------------------------------------------#
#TAD code: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks#
#---------------------------------------------------------------------------------#

# Given a dynamic network dataset, this Python's code creates daily graphs in     #
# edgelist format.                                                                #

####################################################################################################################################################################

###################
#=== Libraries ===#
###################

import numpy as np
import matplotlib.pyplot as plt
import glob 
from datetime import datetime, timedelta

#####################
#=== Input Files ===#
#####################

namePath = 'ReducedDataset100/' 
lisAllFiles = glob.glob(namePath+'RI/*.txt') 
NFiles = len(lisAllFiles) 

####################################################################################################################################################################

###############################
#=== For each file (Token) ===#
############################### 
for nf in range(NFiles):
    print("File: "+lisAllFiles[nf]+" (working)... # "+str(nf+1)+"\n") 
    #-- To open the file and to load the data ---#
    #nf = 0 
    splitFName = lisAllFiles[nf].split('/') # splitFName[-1][:-4]  # Filename (without extension)
    data = np.loadtxt(lisAllFiles[nf]) 
    NRows = len(data) 
    NCols = len(data[0]) 

    #--- To find the start and end dates ---#
    # In timestamp
    tstampStart = data[0][2] 
    tstampEnd = data[-1][2] 
    # In date
    dtobjStart = datetime.fromtimestamp(tstampStart)
    dtobjEnd = datetime.fromtimestamp(tstampEnd)
    # Initial and Final dates (timestamp)
    fechaIni = datetime(dtobjStart.year,dtobjStart.month,dtobjStart.day,0,0,0)
    fechaFin = fechaIni + timedelta(days=1)
    fechaLast = datetime(dtobjEnd.year,dtobjEnd.month,dtobjEnd.day,0,0,0)
    fechaLast = fechaLast + timedelta(days=1) 
    timestampIni = datetime.timestamp(fechaIni) 
    timestampFin = datetime.timestamp(fechaFin) 
    timestampLast = datetime.timestamp(fechaLast) 

    #--- Find the number of nodes ---#
    matINDEX = np.loadtxt(splitFName[0]+'/'+'INDEX/INDEX'+splitFName[-1][2:]) 
    NNodes = int(matINDEX[-1][0]) # There are *** individuals 

    #--- Compute the maximum number of accumulated edges for W ---#
    maxNumContact = 0 
    while timestampFin<=timestampLast:
        # Data in this week: fechaIni to fechaFin
        periodData = data[(data[:,2]>=timestampIni) & (data[:,2]<timestampFin)]
        nRperiod = np.shape(periodData)[0]
        # Fill out Matrices 
        matAD = np.zeros((NNodes, NNodes))
        matWG = np.zeros((NNodes, NNodes))
        for i in range(nRperiod):
            matAD[int(periodData[i,0]-1), int(periodData[i,1]-1)] = 1
            matWG[int(periodData[i,0]-1), int(periodData[i,1]-1)] = matWG[int(periodData[i,0]-1), int(periodData[i,1]-1)] + 1
            # To count edges
            if(maxNumContact<matWG[int(periodData[i,0]-1), int(periodData[i,1]-1)]):
                maxNumContact = matWG[int(periodData[i,0]-1), int(periodData[i,1]-1)] 
        # To continue...
        fechaIni = fechaFin
        fechaFin = fechaFin + timedelta(days=1)
        timestampIni = datetime.timestamp(fechaIni)
        timestampFin = datetime.timestamp(fechaFin)

    #############################
    #=== Creating the graphs ===#
    #############################
    #--- To find the start and end dates (Again...) ---#
    # Initial and Final dates (timestamp)
    fechaIni = datetime(dtobjStart.year,dtobjStart.month,dtobjStart.day,0,0,0)
    fechaFin = fechaIni + timedelta(days=1)
    fechaLast = datetime(dtobjEnd.year,dtobjEnd.month,dtobjEnd.day,0,0,0)
    fechaLast = fechaLast + timedelta(days=1)
    timestampIni = datetime.timestamp(fechaIni)
    timestampFin = datetime.timestamp(fechaFin)
    timestampLast = datetime.timestamp(fechaLast)
     
    #########################
    #=== Create Networks ===# 
    #########################
    idNetwork = 1
    while timestampFin<=timestampLast:
        # Data in this week: fechaIni to fechaFin
        periodData = data[(data[:,2]>=timestampIni) & (data[:,2]<timestampFin)]
        nRperiod = np.shape(periodData)[0]
        # Fill out Matrices #
        matAD = np.zeros((NNodes, NNodes))
        matWG = np.zeros((NNodes, NNodes))
        for i in range(nRperiod):
            matAD[int(periodData[i,0]-1), int(periodData[i,1]-1)] = 1
            matWG[int(periodData[i,0]-1), int(periodData[i,1]-1)] = matWG[int(periodData[i,0]-1), int(periodData[i,1]-1)] + 1
        
        #=== Create Edge List ===#
        fnameAD = splitFName[0]+'/'+'EdgeList_A/A'+splitFName[-1][9:-7]+str(idNetwork)+'.txt'
        fnameWG = splitFName[0]+'/'+'EdgeList_W/W'+splitFName[-1][9:-7]+str(idNetwork)+'.txt'
        fptAD = open(fnameAD, "w") # create new file
        fptWG = open(fnameWG, "w") # create new file
        for c in range(NNodes): 
            for f in range(NNodes): 
                if(matAD[c,f]>0):
                    sAD = str(c+1)+" "+str(f+1)+" "+str(matAD[c,f])+"\n"
                    fptAD.write(sAD) 
                if(matWG[c,f]>0):
                    sWG = str(c+1)+" "+str(f+1)+" "+str(matWG[c,f]/maxNumContact)+"\n"
                    fptWG.write(sWG)     
        fptAD.close()
        fptWG.close()
        #print("Added: "+str(idNetwork)+"\n") 

        #=== To continue... ===#
        fechaIni = fechaFin
        fechaFin = fechaFin + timedelta(days=1)
        timestampIni = datetime.timestamp(fechaIni)
        timestampFin = datetime.timestamp(fechaFin)
        idNetwork = idNetwork + 1
    print("Added: "+str(idNetwork)+" Networks\n") 

####################################################################################################################################################################
