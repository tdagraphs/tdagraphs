mnbs <- function(beta, A, C=1){
  # A function to execute MNBS algorithm (computational cost O(n^3))
  Time <- length(A)
  n <- nrow(A[[1]])
  A_bar = Reduce("+",A)/Time
  D_2 <- matrix(0, n, n) # distance matrix 'D_2'
  
  ## compute distance matrix: MNBS 'D_2'
  A_sq <- (A_bar %*% A_bar)/n
  for (i in 1:(n-2)){
    va1 <- sweep(A_sq[(i+1):n,], 2, A_sq[i,])
    va1[,i] <- 0
    diag(va1[,(i+1):n]) <- 0
    D_2[(i+1):n,i] <- apply(abs(va1),1,max)
  }
  va1 <- A_sq[n,]-A_sq[(n-1),]
  va1[c((n-1),n)] <- 0
  D_2[n,(n-1)] <- max(abs(va1))
  D_2 <- forceSymmetric(t(D_2))
  
  a <- n^(1/2)*(Time^beta)   # neighborhood
  q <- C*(log(n))^0.5/a
  
  ## estimate 'p_tuta'
  qv <- function(x){
    quantileVector <- x<=quantile(x,q)
    return(quantileVector)
  }
  quantileMatrix <- t(apply(D_2, 1, qv))
  neighborNumber <- apply(quantileMatrix,1,sum)
  p_tuta <- t(t(quantileMatrix %*% A_bar)/neighborNumber)

  ## outputs
  res <- vector("list")
  res$P <- p_tuta
  return(res)
}

gmodel.P0 <- function (P, rep = 1, noloop = TRUE){
  # generate adjacency matrix A based on link matrix P
  cond1 = ((all(P >= 0)) && (all(P <= 1)))
  cond2 = (nrow(P) == ncol(P))
  if (!(cond1 && cond2)) {
    stop("* gmodel.P : P is not a valid probability matrix.")
  }
  n = nrow(P)
  if (rep == 1) {
    tmpmat = matrix(runif(n^2), nrow = n)
    G = (tmpmat < P) * 1
    G[lower.tri(G)] <- 0
    G <- t(G) + G - diag(diag(G))
    if (noloop) {
      diag(G) = 0
    }
  }
  else {
    G = list()
    for (i in 1:rep) {
      tmpmat = matrix(runif(n^2), nrow = n)
      tmpG = 1 * (tmpmat < P)
      tmpG[lower.tri(tmpG)] <- 0
      tmpG <- t(tmpG) + tmpG - diag(diag(tmpG))
      if (noloop) {
        diag(tmpG) = 0
      }
      G[[i]] = tmpG
    }
  }
  return(G)
}

d2infty_norm <- function(A){
  # calculate d_{2,infty} norm
  sqrt(max(apply(A^2,1,sum))/dim(A)[1])
}

local_maximizer <- function(D, h){
  # A function to find the local-maximizer of D-measure
  Time <- length(D)
  T0 <- h; T1 <- Time-h
  LM <- rep(0, Time)
  for (t in T0:T1){
    if (D[t]==max(D[max(1,(t-h+1)):min(Time,(t+h-1))])){
      LM[t] <- 1
    }
  }
  LM <- LM==1
  return(LM)
}

SaRa <- function(D, h, threshold){
  # A function to find perform change-point detection of D-measure via SaRa
  # SaRa algorithm: Screening and Ranking algorithm
  Time <- length(D)
  T0 <- h; T1 <- Time-h
  cp <- rep(0, Time)
  for (t in T0:T1){
    if ((D[t]==max(D[max(1,(t-h+1)):min(Time,(t+h-1))]))&(D[t]>threshold) ){
      cp[t] <- 1
    }
  }
  cp <- which(cp==1)
  if(length(cp)==0){
    cp <- NULL
  }
  return(cp)
}

binaryS <- function(E, h, reference_point=0){
  # A binary segmentation for ChenZhang(2015) with Type I error 0.05
  # BS with ChenZhang(2015) for multiple change-point detection
  dist_frobenius <- dist(E)
  similarity_measure <- mstree(dist_frobenius)
  Time <- dim(E)[1]
  if((Time-h+1)<=h){
    return(NA)
  }
  # print(c(h,Time-h+1))
  r1 <- gseg1(Time, similarity_measure, n0=h, n1=Time-h+1, statistics="original")
  #if(r1$pval.appr$ori<0.05){
  if(r1$pval.appr$ori<(0.05/6)){
    tauhat <- r1$scanZ$ori$tauhat
    return(c(reference_point+tauhat, binaryS(E[1:tauhat,], h, reference_point=reference_point), 
             binaryS(E[(tauhat+1):Time,], h, reference_point=reference_point+tauhat)))
  }else{
    return(NA)
  }
}

boysen_dist <- function(est_cp, true_cp){
  # Boysen distance for evaluating estimated change-point accuracy
  if(length(true_cp)==0){
    if(length(est_cp)==0||is.na(est_cp)){
      return(c(dist1=0, dist2=0))
    }
    return(c(dist1=max(est_cp), dist2=NA))
  }
  if(length(est_cp)==0||is.na(est_cp)){
    return(c(dist1=NA, dist2=max(true_cp)))
  }
  dist1 <- dist2 <- 0
  for(cp in est_cp){ # measure error by over-estimation of change-point
    current_dist1 <- min(abs(true_cp-cp))
    dist1 <- max(dist1, current_dist1)
  }
  for(cp in true_cp){ # measure error by under-estimation of change-point
    current_dist2 <- min(abs(est_cp-cp))
    dist2 <- max(dist2, current_dist2)
  }
  return(c(dist1,dist2))
}

change_scenario <- function(n, Time, scenario){
  # Generate different types of Dynamic SBM change-point models
  # A function to simulate different scenarios of change-points in networks
  p1 <- matrix(0, ncol=n, nrow=n)
  p2 <- matrix(0, ncol=n, nrow=n)
  if(scenario=='NoCP1'){ # ZhangLevinaZhu(2017) Graphon-I
    K <- floor(log(n))
    block_size <- floor(n/K)
    for(no_block in 1:(K-1)){
      p1[((no_block-1)*block_size+1):(no_block*block_size), ((no_block-1)*block_size+1):(no_block*block_size)] <- 
        no_block/(K+1)
    }
    p1[((no_block*block_size)+1):n, ((no_block*block_size)+1):n] <- K/(K+1)
    p1[p1==0] <- 0.3/(K+1)
    p2 <- p1
    return(list(p1=p1, p2=p2))
  }
  
  if(scenario=='NoCP2'){ # ZhangLevinaZhu(2017) Graphon-II
    for(i in 1:n){
      for(j in 1:n){
        u <- i/(n+1)
        v <- j/(n+1)
        p1[i,j] <- sin(5*pi*(u+v-1)+1)/2+0.5
      }
    }
    p2 <- p1
    return(list(p1=p1, p2=p2))
  }
  
  if(scenario=='NoCP3'){ # ZhangLevinaZhu(2017) Graphon-III
    for(i in 1:n){
      for(j in 1:n){
        u <- i/(n+1)
        v <- j/(n+1)
        p1[i,j] <- (u^2+v^2)/3*cos(1/(u^2+v^2))+0.15
      }
    }
    p2 <- p1
    return(list(p1=p1, p2=p2))
  }

  if(scenario==3){ # DSBM switch community (one CP) with two blocks
    # DSBM-IV
    # ~ 0.3^2 for d2infty^2
    # ~ 0.3^2*1/T^(1/4)/n^(1/3) for Frobenius^2
    deltaT <- 1/4 # control the probability change
    deltan <- 1/3 # control the size of block that changes
    
    switch_portion <- 2*1/Time^(deltaT)/n^(deltan)
    block_1 <- floor(2*n/3) # block 1 size
    p1_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    
    # p1
    p1[1:block_1, 1:block_1] <-  p1_block[1,1]
    p1[1:block_1, (block_1+1):n] <- p1[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p1[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    
    # p2
    block_1 <- floor(2*n/3*(1-switch_portion)) # block 1 changes
    p2[1:block_1, 1:block_1] <-  p1_block[1,1]
    p2[1:block_1, (block_1+1):n] <- p2[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p2[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario==4){ # DSBM change connectivity (one CP) with two blocks
    # DSBM-V
    # ~ 1/T^(1/4)/n^(1/4) for d2infty^2
    # ~ 1/T^(1/4)/n^(1/2) for Frobenius^2
    deltaT <- 1/8 # control the probability change
    deltan <- 3/4 # control the size of block that changes
    
    block_1 <- floor(n^{deltan}) # block 1 size
    p1_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    
    # p1
    p1[1:block_1, 1:block_1] <-  p1_block[1,1]
    p1[1:block_1, (block_1+1):n] <- p1[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p1[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    
    # p2
    p2[1:block_1, 1:block_1] <-  p1_block[1,1] - Time^(-deltaT) # the only change
    p2[1:block_1, (block_1+1):n] <- p2[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p2[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario==5){ # DSBM community merging from three to two blocks (or splitting from two to three blocks)
    # DSBM-I
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/3)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p1_block <- cbind(c(0.6, 0.6-change, 0.3), c(0.6-change, 0.6, 0.3), c(0.3, 0.3, 0.6)) 
    p2_block <- cbind(c(0.6, 0.6, 0.3), c(0.6, 0.6, 0.3), c(0.3, 0.3, 0.6)) 
    
    # p1
    p1[1:block_1, 1:block_1] <- p1_block[1,1]
    p1[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p1_block[2,2]
    p1[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p1_block[3,3]
    p1[1:block_1, (block_1+1):(block_1+block_2)] <- p1[(block_1+1):(block_1+block_2), 1:block_1] <- p1_block[1,2]
    p1[1:block_1, (block_1+block_2+1):n] <- p1[(block_1+block_2+1):n, 1:block_1] <- p1_block[1,3]
    p1[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p1[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p1_block[2,3]
    
    # p2
    p2[1:block_1, 1:block_1] <- p2_block[1,1]
    p2[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p2_block[2,2]
    p2[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p2_block[3,3]
    p2[1:block_1, (block_1+1):(block_1+block_2)] <- p2[(block_1+1):(block_1+block_2), 1:block_1] <- p2_block[1,2]
    p2[1:block_1, (block_1+block_2+1):n] <- p2[(block_1+block_2+1):n, 1:block_1] <- p2_block[1,3]
    p2[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p2[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p2_block[2,3]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario==6){ # DSBM change connectivity
    # DSBM-II
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/3)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p1_block <- cbind(c(0.6, 0.6-change, 0.3), c(0.6-change, 0.6, 0.3), c(0.3, 0.3, 0.6)) 
    p2_block <- cbind(c(0.6+change, 0.6-change, 0.3), c(0.6-change, 0.6+change, 0.3), c(0.3, 0.3, 0.6)) 
    
    # p1
    p1[1:block_1, 1:block_1] <- p1_block[1,1]
    p1[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p1_block[2,2]
    p1[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p1_block[3,3]
    p1[1:block_1, (block_1+1):(block_1+block_2)] <- p1[(block_1+1):(block_1+block_2), 1:block_1] <- p1_block[1,2]
    p1[1:block_1, (block_1+block_2+1):n] <- p1[(block_1+block_2+1):n, 1:block_1] <- p1_block[1,3]
    p1[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p1[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p1_block[2,3]
    
    # p2
    p2[1:block_1, 1:block_1] <- p2_block[1,1]
    p2[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p2_block[2,2]
    p2[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p2_block[3,3]
    p2[1:block_1, (block_1+1):(block_1+block_2)] <- p2[(block_1+1):(block_1+block_2), 1:block_1] <- p2_block[1,2]
    p2[1:block_1, (block_1+block_2+1):n] <- p2[(block_1+block_2+1):n, 1:block_1] <- p2_block[1,3]
    p2[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p2[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p2_block[2,3]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario==7){ # DSBM change connectivity
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/3)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p1_block <- cbind(c(0.6, 0.6, 0.3), c(0.6, 0.6, 0.3), c(0.3, 0.3, 0.6)) 
    p2_block <- cbind(c(0.6+change, 0.6, 0.3), c(0.6, 0.6+change, 0.3), c(0.3, 0.3, 0.6)) 
    
    # p1
    p1[1:block_1, 1:block_1] <- p1_block[1,1]
    p1[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p1_block[2,2]
    p1[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p1_block[3,3]
    p1[1:block_1, (block_1+1):(block_1+block_2)] <- p1[(block_1+1):(block_1+block_2), 1:block_1] <- p1_block[1,2]
    p1[1:block_1, (block_1+block_2+1):n] <- p1[(block_1+block_2+1):n, 1:block_1] <- p1_block[1,3]
    p1[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p1[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p1_block[2,3]
    
    # p2
    p2[1:block_1, 1:block_1] <- p2_block[1,1]
    p2[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p2_block[2,2]
    p2[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p2_block[3,3]
    p2[1:block_1, (block_1+1):(block_1+block_2)] <- p2[(block_1+1):(block_1+block_2), 1:block_1] <- p2_block[1,2]
    p2[1:block_1, (block_1+block_2+1):n] <- p2[(block_1+block_2+1):n, 1:block_1] <- p2_block[1,3]
    p2[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p2[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p2_block[2,3]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario==8){ # Switch community even and odd
    # DSBM-VI
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/2)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p1_block <- cbind(c(0.6, 0.6-change), c(0.6-change, 0.6)) 
    
    # p1
    p1[1:block_1, 1:block_1] <- p1_block[1,1]
    p1[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p1_block[2,2]
    p1[1:block_1, (block_1+1):(block_1+block_2)] <- p1[(block_1+1):(block_1+block_2), 1:block_1] <- p1_block[1,2]
    
    # p2
    p2[which((1:n%%2==1)%*%t((1:n%%2==1))==1)] <- p1_block[1,1]
    p2[which((1:n%%2==0)%*%t((1:n%%2==0))==1)] <- p1_block[2,2]
    p2[which((1:n%%2==1)%*%t((1:n%%2==0))==1)] <- p2[which((1:n%%2==0)%*%t((1:n%%2==1))==1)] <-
      p1_block[1,2]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario==9){ # DSBM switch community (one CP) with two blocks
    # DSBM-III
    # ~ 1/n^(1/3)/T^(1/4) for d2infty^2
    # ~ 1/n^(4/3)/T^(1/4) for Frobenius^2
    
    switch_no_node <- 1
    block_1 <- floor(2*n/3) # block 1 size
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p1_block <- cbind(c(0.6, 0.6-change), c(0.6-change, 0.6)) # change even more
    
    # p1
    p1[1:block_1, 1:block_1] <-  p1_block[1,1]
    p1[1:block_1, (block_1+1):n] <- p1[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p1[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    
    # p2
    block_1 <- floor(2*n/3) - switch_no_node # block 1 changes
    p2[1:block_1, 1:block_1] <-  p1_block[1,1]
    p2[1:block_1, (block_1+1):n] <- p2[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p2[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    return(list(p1=p1,p2=p2))
  }
  
  if(scenario=='MCP3'){ # Multiple change-point with three change-points
    p1 <- matrix(0, ncol=n, nrow=n)
    p2 <- matrix(0, ncol=n, nrow=n)
    p3 <- matrix(0, ncol=n, nrow=n)
    p4 <- matrix(0, ncol=n, nrow=n)
    
    # 1st change-point DSBM switch community with two blocks
    # ~ 0.3^2 for d2infty^2
    # ~ 2*0.3^2*1/T^(1/4)/n^(1/3) for Frobenius^2
    deltaT <- 1/4 # control the probability change
    deltan <- 1/3 # control the size of block that changes
    
    # p1
    switch_portion <- 2*1/Time^(deltaT)/n^(deltan)
    block_1 <- 2*floor(n/3*(1-switch_portion)) # block 1 changes
    p1_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    p1[1:block_1, 1:block_1] <-  p1_block[1,1]
    p1[1:block_1, (block_1+1):n] <- p1[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p1[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    
    # p2
    block_1 <- 2*floor(n/3) # block 1 size
    p2_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    p2[1:block_1, 1:block_1] <-  p2_block[1,1]
    p2[1:block_1, (block_1+1):n] <- p2[(block_1+1):n, 1:block_1] <- p2_block[1,2]
    p2[(block_1+1):n, (block_1+1):n] <- p2_block[2,2]
    
    # 2nd change-point DSBM splitting community from two blocks to three blocks
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/3)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p3_block <- cbind(c(0.6, 0.6-change, 0.3), c(0.6-change, 0.6, 0.3), c(0.3, 0.3, 0.6)) 
    p3[1:block_1, 1:block_1] <- p3_block[1,1]
    p3[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p3_block[2,2]
    p3[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p3_block[3,3]
    p3[1:block_1, (block_1+1):(block_1+block_2)] <- p3[(block_1+1):(block_1+block_2), 1:block_1] <- p3_block[1,2]
    p3[1:block_1, (block_1+block_2+1):n] <- p3[(block_1+block_2+1):n, 1:block_1] <- p3_block[1,3]
    p3[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p3[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p3_block[2,3]
    
    # 3rd change-point DSBM change connectivity
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/3)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p4_block <- cbind(c(0.6+change, 0.6-change, 0.3), c(0.6-change, 0.6+change, 0.3), c(0.3, 0.3, 0.6)) 
    p4[1:block_1, 1:block_1] <- p4_block[1,1]
    p4[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p4_block[2,2]
    p4[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p4_block[3,3]
    p4[1:block_1, (block_1+1):(block_1+block_2)] <- p4[(block_1+1):(block_1+block_2), 1:block_1] <- p4_block[1,2]
    p4[1:block_1, (block_1+block_2+1):n] <- p4[(block_1+block_2+1):n, 1:block_1] <- p4_block[1,3]
    p4[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p4[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p4_block[2,3]
    return(list(p1=p1,p2=p2,p3=p3,p4=p4))
  }
  
  if(scenario=='MCP4'){ # Multiple change-point with three change-points
    p1 <- matrix(0, ncol=n, nrow=n)
    p2 <- matrix(0, ncol=n, nrow=n)
    p3 <- matrix(0, ncol=n, nrow=n)
    p4 <- matrix(0, ncol=n, nrow=n)
    p5 <- matrix(0, ncol=n, nrow=n)
    
    # 1st change-point DSBM switch community with two blocks
    # ~ 0.3^2 for d2infty^2
    # ~ 2*0.3^2*1/T^(1/4)/n^(1/3) for Frobenius^2
    deltaT <- 1/4 # control the probability change
    deltan <- 1/3 # control the size of block that changes
    
    # p1
    switch_portion <- 2*1/Time^(deltaT)/n^(deltan)
    block_1 <- 2*floor(n/3*(1-switch_portion)) # block 1 changes
    p1_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    p1[1:block_1, 1:block_1] <-  p1_block[1,1]
    p1[1:block_1, (block_1+1):n] <- p1[(block_1+1):n, 1:block_1] <- p1_block[1,2]
    p1[(block_1+1):n, (block_1+1):n] <- p1_block[2,2]
    
    # p2
    block_1 <- 2*floor(n/3) # block 1 size
    p2_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    p2[1:block_1, 1:block_1] <-  p2_block[1,1]
    p2[1:block_1, (block_1+1):n] <- p2[(block_1+1):n, 1:block_1] <- p2_block[1,2]
    p2[(block_1+1):n, (block_1+1):n] <- p2_block[2,2]
    
    # 2nd change-point DSBM splitting community from two blocks to three blocks
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- block_2 <- floor(n/3)
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    p3_block <- cbind(c(0.6, 0.6-change, 0.3), c(0.6-change, 0.6, 0.3), c(0.3, 0.3, 0.6)) 
    p3[1:block_1, 1:block_1] <- p3_block[1,1]
    p3[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p3_block[2,2]
    p3[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p3_block[3,3]
    p3[1:block_1, (block_1+1):(block_1+block_2)] <- p3[(block_1+1):(block_1+block_2), 1:block_1] <- p3_block[1,2]
    p3[1:block_1, (block_1+block_2+1):n] <- p3[(block_1+block_2+1):n, 1:block_1] <- p3_block[1,3]
    p3[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p3[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p3_block[2,3]
    
    # 3rd change-point DSBM merging community from three blocks to two blocks
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    block_1 <- 2*floor(n/3) # block 1 size
    p4_block <- cbind(c(0.6, 0.3), c(0.3, 0.6))
    p4[1:block_1, 1:block_1] <-  p4_block[1,1]
    p4[1:block_1, (block_1+1):n] <- p4[(block_1+1):n, 1:block_1] <- p4_block[1,2]
    p4[(block_1+1):n, (block_1+1):n] <- p4_block[2,2]
    
    # 4th change-point DSBM splitting community
    # ~ 1/3/T^(1/4)/n^(1/3) for d2infty^2
    # ~ 2/9/T^(1/4)/n^(1/3) for Frobenius^2
    deltaT <- 1/8 # control the probability change
    deltan <- 1/6 # control the size of block that changes
    change <- 1/n^deltan/Time^deltaT
    # block_1 <- 2*floor(n/3) # block 1 size
    # p5_block <- cbind(c(0.6, 0.3), c(0.3, 0.6-change))
    # p5[1:block_1, 1:block_1] <-  p5_block[1,1]
    # p5[1:block_1, (block_1+1):n] <- p5[(block_1+1):n, 1:block_1] <- p5_block[1,2]
    # p5[(block_1+1):n, (block_1+1):n] <- p5_block[2,2]
    
    block_1 <- block_2 <- floor(n/3)
    p5_block <- cbind(c(0.6+change, 0.6, 0.3), c(0.6, 0.6+change, 0.3), c(0.3, 0.3, 0.6)) 
    p5[1:block_1, 1:block_1] <- p5_block[1,1]
    p5[(block_1+1):(block_1+block_2), (block_1+1):(block_1+block_2)] <- p5_block[2,2]
    p5[(block_1+block_2+1):n, (block_1+block_2+1):n] <- p5_block[3,3]
    p5[1:block_1, (block_1+1):(block_1+block_2)] <- p5[(block_1+1):(block_1+block_2), 1:block_1] <- p5_block[1,2]
    p5[1:block_1, (block_1+block_2+1):n] <- p5[(block_1+block_2+1):n, 1:block_1] <- p5_block[1,3]
    p5[(block_1+1):(block_1+block_2), (block_1+block_2+1):n] <- p5[(block_1+block_2+1):n, (block_1+1):(block_1+block_2)] <- p5_block[2,3]
    
    return(list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5))
  }
}
