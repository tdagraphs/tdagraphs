----------------------------------------------------------------------------
TAD: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks
----------------------------------------------------------------------------

Our source code was developed using the following list of software.

1) * Matlab *
A basic, standard, installation of Matlab is required to reduce the network size via maximum weight subgraph approximation.  

2) * Python *
We uses python 3.7.5 in several source codes, also we have tested such scripts in the latest python version 3.8.3.
We only use common libraries such as: numpy, matplotlib, os, glob, datetime.

4) * Aleph - A library for exploring persistent homology *
Compute the clique community persistence diagrams using the Alephs software. See python's examples in the website https://pseudomanifold.github.io/Aleph/Rieck17d.html.

3) * R *
Our computations of Wasserstein distance, anomaly detection and geodesic densification were done via R. We installed the anomaly detection library via install.packages("devtools") and devtools::install_github("twitter/AnomalyDetection"); and following libraries via the R Studio's packaged handler: MASS, dodgr, stringr, R.matlab, TDA, caret, e1071. 


README.txt includes the description of each source code and folder attached in this submission.
