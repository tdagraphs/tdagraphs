#--------------------------------------------------------------------------------------------#
#TAD code: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks           #
#--------------------------------------------------------------------------------------------#

# Clique Community Persistence: A Topological Visual Analysis Approach for Complex Networks: #
# This Python script generates persistence diagrams using the Clique Community Persistence   #
# Implementation of Aleph library (see 'Requirements.txt).                                   #
 
######################################################################################################################################################

 
######################################################################################################################################################

###################
#=== Libraries ===#
###################

import numpy as np
import networkx as nx 
import dionysus as d 
import matplotlib.pyplot as plt
import os 
import glob 


####################
#=== Parameters ===#
####################

# 1. Unweighted #
#nameFolderIn = 'EdgeList_A_GD' 
#nameFolderOut = 'RIPPLE_PD_A_GD' 

# 2. Weighted #
nameFolderIn = 'EdgeList_W_GD' 
nameFolderOut = 'RIPPLE_PD_W_GD' 

######################################################################################################################################################

KMaxClique = '4' # Usually 
routeAllFiles = '/.../RIPPLE_5Y90/'

# Where is Aleph installed?
routeCMD = '/.../Aleph-master/build/tools/./clique_persistence_diagram'


###################
#=== Full code ===#
###################

#lisAllFiles = glob.glob('Aleph-master/build/tools/'+nameFolderIn+'/*.txt') 
lisAllFiles = glob.glob(routeAllFiles+nameFolderIn+'/*.txt') 
lisAllFiles.sort()
NFiles = len(lisAllFiles) 
for nf in range(NFiles):
    print("File["+str(nf+1)+"]... name: "+lisAllFiles[nf]+"\n") 
    data = np.loadtxt(lisAllFiles[nf]) # Load data
    splitFName = lisAllFiles[nf].split('/')
    if not data.any(): # Empty file...
        maxVal = 0 
        #fnamePD = 'Aleph-master/build/tools/'+nameFolderOut+'/'+splitFName[-1][:-4]+'_PD.txt' # Filename
        fnamePD = routeAllFiles+nameFolderOut+'/'+splitFName[-1][:-4]+'_PD.txt' # Filename
        matP = np.reshape([0, 0, 0], (-1,3))
        np.savetxt(fnamePD, matP) # Save PD 
    else: # There is elements in the edge list
        if(data.ndim==1): # It is an array
            maxVal = data[2] # Maximum weight 
        else: # It is a matrix
            maxVal = max(data[:,2]) # Maximum weight 
        
        #strCmd = 'Aleph-master/build/tools/./clique_persistence_diagram --ignore-empty --invert-weights '+lisAllFiles[nf]+' 20' 
        strCmd = routeCMD+' --ignore-empty --invert-weights '+lisAllFiles[nf]+' '+KMaxClique 
        os.system(strCmd)
        #fnamePD = 'Aleph-master/build/tools/'+nameFolderOut+'/'+splitFName[-1][:-4]+'_PD.txt' # Filename
        fnamePD = routeAllFiles+nameFolderOut+'/'+splitFName[-1][:-4]+'_PD.txt' # Filename
        lisFiles = glob.glob('/tmp/'+splitFName[-1][:-4]+'_k*.txt') 
        lisFiles.sort() 
        matr = [] 
        valD = 0 # Always 0... the 'intervals' shows when a k-community arise and die
        for fnameK in lisFiles: # Filling PDs
            fpt = open(fnameK, 'r') 
            flagPD = -1
            for line in fpt: 
                lisBreaks = line.split()
                if(lisBreaks[0]!=str("#")):
                    flagPD = 1
                if(flagPD==1):
                    valBirth = float(lisBreaks[0])
                    valDeath = float(lisBreaks[1])
                    if(valDeath>(1.5*maxVal)):
                        valDeath = -1
                    matr = np.concatenate((matr, [valD, valBirth, valDeath]))   
        matP = np.reshape(matr, (-1,3))
        np.savetxt(fnamePD, matP) # Save PD
        # Deleting all files...
        for fnameK in lisFiles: # Filling PDs    
            os.remove(fnameK) # To delete file #
        os.remove('/tmp/'+splitFName[-1])

######################################################################################################################################################
