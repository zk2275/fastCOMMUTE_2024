woodbury = function (A, U, VT){
  first<-ginv(A)
  second<-ginv(A)%*%U%*%ginv(diag(nrow(VT))+VT%*%ginv(A)%*%U)%*%VT%*%ginv(A)
  result<-first-second
  return(result)
}
# U= a constant is not suggested


svd_transfer = function(X,y,Ratio=0.1){
  
  
    svd_result <- svd(X)
    X.temp <- X
    
    U <- svd_result$u
    d <- svd_result$d
    V <- svd_result$v
    m <- nrow(X)
    n <- ncol(X)
    kept <- Ratio*n
    
    ### Compression
    U<-U[,1:kept]
    
    d<-d[1:kept]
    
    Sigma <- diag(d)
    
    V<-V[,1:kept]
    ###
    

    X_new <-  Sigma%*%t(V)
    y_new <-  ginv(Sigma)%*%t(V) %*% t(X.temp) %*% y
    results <- list(X_new,y_new)

  return(results)
}

beta_fastCOMMUTE_generation = function(K, method = "fastcommute" ,nov = "small",X.tar,X.src,y.tar,y.src,delta,wthreshold,beta.true) {
  
  if(nov == "large")
  { 
    Ratio <- 1
    for (i in 1:K) {
      # Compression
      distillation_result <- svd_transfer(X.src[[i]],y.src[[i]],Ratio = 0.1)
      X.src[[i]] <- distillation_result[[1]]
      print("X.src distillation ready")
      y.src[[i]] <- distillation_result[[2]]
      print("y.src distillation ready")
    }
  }
  
 
  # After Compression
  
  nZ_score <- Z_score <- GKoff_Z_score <- GKoff <- list()
  XTX_tar <- t(X.tar)%*%X.tar
  XTY_tar <- t(X.tar)%*%y.tar
  scale_src_tar <- XTX_src <- XTY_src <- XTX_src_delta <- 
  scale_XTY_tar<- XTX_tar_delta <- scale_XTX_tar_delta<- list()
  
  # beta.tar
  beta.tar = ginv(t(X.tar)%*%X.tar)%*%t(X.tar)%*%y.tar
  
  # prepare for matrices
    for (i in 1:K) {
      XTX_src[[i]] <- t(X.src[[i]])%*%X.src[[i]]
      XTY_src[[i]] <- t(X.src[[i]])%*%y.src[[i]]
      scale_src_tar[[i]] <- (nrow(X.src[[i]])-1)/(nrow(X.tar)-1)
      XTX_src_delta[[i]] <- t(X.src[[i]])%*%X.src[[i]]%*%delta.TL[[i]]
      XTX_tar_delta[[i]] <- t(X.tar)%*%X.tar%*%delta.TL[[i]] # may be unnecessary
      
      scale_XTX_tar_delta[[i]] <- scale_src_tar[[i]]*t(X.tar)%*%X.tar%*%delta.TL[[i]]
      scale_XTY_tar[[i]]  <- scale_src_tar[[i]]*t(X.tar)%*%y.tar
    }
# Common situation with a small number of p    
if(nov == "small")    { 
  #--------fastCOMMUTE--------#
  if(method == "fastcommute"){
    
    for (i in 1:K) {
      nZ_score[[i]] <-  scale_src_tar[[i]]*t(X.tar)%*%X.tar%*%wthreshold[[i]]
      Z_score[[i]] <- 1/sqrt(nrow(X.tar))*nZ_score[[i]]
      
    }
    
  }
  else{
    print("Other methods need to be expored")
  }
  #--------fastCOMMUTE--------#
  
  
  ##beta_commute0
  sum0 <- list()
  for (i in 1:K) {
    sum0[[i]] <- nZ_score[[i]] + scale_src_tar[[i]]*XTX_tar_delta[[i]] #1129 scale_src_tar[[i]]*XTX_tar_delta
  }
  
  #sum_of_list <- Reduce(`+`, my_list)
  wb0<-woodbury(t(X.tar)%*%X.tar, Reduce(`+`, scale_src_tar)*t(X.tar), X.tar)
  beta_commute0 <- wb0%*% # ginv((1+Reduce(`+`, scale_src_tar)  )*XTX_tar  )
    ( XTY_tar +  Reduce(`+`, sum0) )
  
  ##beta_commute1-K
  beta_commute <- list()
  for (k in 1:K) {
    
    # calibrated
    sum_calibrated <- Reduce(`+`, XTY_src[1:k])+
      Reduce(`+`, XTX_src_delta[1:k])
    
    if(k==K){
      #combine
      
      beta_commute[[k]] <- ginv(XTX_tar + Reduce(`+`, XTX_src[1:k]) )%*%
        ( XTY_tar + sum_calibrated)
    }
    else{ 
      
      # syn
      sum_syn <- Reduce(`+`, nZ_score[(k+1):K])+ # 1129 nz_score
        Reduce(`+`, scale_XTX_tar_delta[(k+1):K]) # 1129 1129 scale_src_tar[[i]]*XTX_tar_delta
      
      #combine
      wbk<-woodbury(t(X.tar)%*%X.tar+Reduce(`+`, XTX_src[1:k]), Reduce(`+`, scale_src_tar[(k+1):K])*t(X.tar), X.tar)
      beta_commute[[k]] <- wbk%*% # ginv((1+Reduce(`+`, scale_src_tar[(k+1):K]) )*XTX_tar + Reduce(`+`, XTX_src[1:k]) )
        (XTY_tar + sum_calibrated + sum_syn)
    }
    #+Reduce(`+`, scale[(k+1):K]) )*XTX_tar
    
  }
  ### !!!put k=0 to the last one
  beta_commute[[K+1]]<-beta_commute0
  return (beta_commute)
# The end of a common situation  
}
    # Complex situation with a large number of p    
else if(nov == "large")  {
  # creating lambda*I(Ridge Regression)
  I = diag(ncol(X.tar))
  lambda = 0.0001
  adjmatrix <- lambda*I
  #--------fastCOMMUTE--------#
  if(method == "fastcommute"){
    
    for (i in 1:K) {
      nZ_score[[i]] <-  scale_src_tar[[i]]*t(X.tar)%*%X.tar%*%wthreshold[[i]]
      Z_score[[i]] <- 1/sqrt(nrow(X.tar))*nZ_score[[i]]
      
    }
    
  }
  else{
    print("Other methods need to be expored")
  }
  #--------fastCOMMUTE--------#
  
  
  # beta_commute0
  sum0 <- list()
  for (i in 1:K) {
    sum0[[i]] <- nZ_score[[i]] + scale_src_tar[[i]]*XTX_tar_delta[[i]] #1129 scale_src_tar[[i]]*XTX_tar_delta
  }
  
  ## sum_of_list <- Reduce(`+`, my_list)
  A = (1+Reduce(`+`, scale_src_tar))*t(X.tar)%*%X.tar
  
  ## Get lambda for beta_commute0
  x<-X.tar;y<-y.tar
  lambda = cv.glmnet(x,y)$lambda.min;adjmatrix <- lambda*I
  
  
  wb0<-woodbury(A, adjmatrix , I )
  
  beta_commute0 <- wb0%*% ( XTY_tar +  Reduce(`+`, sum0) )
    
  
  ##beta_commute1-K
  beta_commute <- list()
  
  ## Get lambda for beta_commute0-K
  x<-rbind(X.tar,X.src[[1]],X.src[[2]],X.src[[3]]);y<-rbind(y.tar,y.src[[1]],y.src[[2]],y.src[[3]])
  lambda = cv.glmnet(x,y)$lambda.min;adjmatrix <- lambda*I
  
  for (k in 1:K) {
    
    # calibrated
    sum_calibrated <- Reduce(`+`, XTY_src[1:k])+
      Reduce(`+`, XTX_src_delta[1:k])
    
    if(k==K){
      # combine

      
      beta_commute[[k]] <- ginv(XTX_tar + Reduce(`+`, XTX_src[1:k])+adjmatrix )%*%
        ( XTY_tar + sum_calibrated)
    }
    else{ 
      
      # syn
      sum_syn <- Reduce(`+`, nZ_score[(k+1):K])+ # 1129 nz_score
        Reduce(`+`, scale_XTX_tar_delta[(k+1):K]) # 1129 1129 scale_src_tar[[i]]*XTX_tar_delta
      
      #combine
      A = t(X.tar)%*%X.tar+Reduce(`+`, XTX_src[1:k])+ Reduce(`+`, scale_src_tar[(k+1):K])*t(X.tar)%*%X.tar
      
      
      
      wbk<-woodbury(A,adjmatrix,I)
      beta_commute[[k]] <- wbk%*% # ginv((1+Reduce(`+`, scale_src_tar[(k+1):K]) )*XTX_tar + Reduce(`+`, XTX_src[1:k]) )
        (XTY_tar + sum_calibrated + sum_syn)
    }
    #+Reduce(`+`, scale[(k+1):K]) )*XTX_tar
    
  }
 
  
  
  #####################
  ### aggregation
  #####################
  X.pool2 = X.pool[,1:p]
  X.til = Data.gen.one(n.tar=100, X.pool2)$X.tar
  # include intercepts
  X.til = cbind(rep(1,nrow(X.til)),X.til)
  y.til <- X.til %*% as.matrix(beta.true) + rnorm(nrow(X.til))
  
  
  ### fastCOMMUTE M=0 (federated)
  B.syn.M0 = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute0)
  wt.syn.M0.agg.new <- Agg.fun.new(B.syn.M0, X.til, y.til, const=1)
  beta_commute0 <- B.syn.M0%*%wt.syn.M0.agg.new
  
  ### COMMUTE M=1
  B.syn.M1 = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute[[1]])
  wt.syn.M1.agg.new <- Agg.fun.new(B.syn.M1, X.til, y.til, const=1)
  beta_commute[[1]] <- B.syn.M1%*%wt.syn.M1.agg.new
  
  ### COMMUTE M=2
  B.syn.M2 = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute[[2]])
  wt.syn.M2.agg.new <- Agg.fun.new(B.syn.M2, X.til, y.til, const=1)
  beta_commute[[2]] <- B.syn.M2%*%wt.syn.M2.agg.new
  
  ### COMMUTE M=3 (pooledTL)
  B.TL = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute[[3]])
  wt.TL.agg.new <- Agg.fun.new(B.TL, X.til, y.til, const=1)
  beta_commute[[3]] <- B.TL%*%wt.TL.agg.new
  
  
  ### !!!put k=0 to the last one
  beta_commute[[K+1]]<-beta_commute0
  return (beta_commute)
}
   
else{
    print("Please correctly specify the nov, large or small.")
}
    
    
}


# beta_fastDD_generation = function(K,X.tar,X.src,y.tar,y.src,delta,wthreshold,beta.true) {
#   # Compression
#  
#     Ratio <- 1
#     for (i in 1:K) {
#       ###
#       svd_result <- svd(X.src[[i]])
#       
#       U <- svd_result$u
#       d <- svd_result$d
#       V <- svd_result$v
#       m <- nrow(X.src[[i]])
#       n <- ncol(X.src[[i]])
#       kept <- Ratio*n
#       
#       # Compression
#       # U[,80:min(n,m)] <- 0
#       U<-U[,1:kept]
#       # d[80:min(n,m)] <- 0
#       d<-d[1:kept]
#       
#       Sigma <- diag(d)
#       
#       # V[,80:min(n,m)] <- 0
#       V<-V[,1:kept]
#       ###
#       X.temp <- X.src[[i]]
#       X.src[[i]] <- Sigma%*%t(V)
#       y.src[[i]] <-  ginv(Sigma)%*%t(V) %*% t(X.temp) %*% y.src[[i]]
#     }
#   
#   
#   # After Compression
#   
#   nZ_score <- Z_score <- GKoff_Z_score <- GKoff <- list()
#   XTX_tar <- t(X.tar)%*%X.tar
#   XTY_tar <- t(X.tar)%*%y.tar
#   scale_src_tar <- XTX_src <- XTY_src <- XTX_src_delta <- 
#     scale_XTY_tar<- XTX_tar_delta <- scale_XTX_tar_delta<- list()
#   
#   # beta.tar
#   beta.tar = ginv(t(X.tar)%*%X.tar)%*%t(X.tar)%*%y.tar
#   
#   # prepare for matrices
#   for (i in 1:K) {
#     XTX_src[[i]] <- t(X.src[[i]])%*%X.src[[i]]
#     XTY_src[[i]] <- t(X.src[[i]])%*%y.src[[i]]
#     XTX_src_delta[[i]] <- t(X.src[[i]])%*%X.src[[i]]%*%delta.TL[[i]]
#     XTX_tar_delta[[i]] <- t(X.tar)%*%X.tar%*%delta.TL[[i]]
#   }
# 
# 
#     # creating lambda*I(Ridge Regression)
#     I = diag(ncol(X.tar))
#     lambda = 0.0001
#     adjmatrix <- lambda*I
#     
#       for (i in 1:K) {
#         nZ_score[[i]] <-  t(X.src[[i]])%*%X.src[[i]]%*%wthreshold[[i]]
#        
#       }
#     
#     
#     # beta_commute0
#     sum0 <- list()
#     for (i in 1:K) {
#       sum0[[i]] <- nZ_score[[i]] + scale_src_tar[[i]]*XTX_tar_delta[[i]] #1129 scale_src_tar[[i]]*XTX_tar_delta
#     }
#     
#     ## sum_of_list <- Reduce(`+`, my_list)
#     A = (1+Reduce(`+`, scale_src_tar))*t(X.tar)%*%X.tar
#     
#     ## Get lambda for beta_commute0
#     x<-X.tar;y<-y.tar
#     lambda = cv.glmnet(x,y)$lambda.min;adjmatrix <- lambda*I
#     
#     
#     wb0<-woodbury(A, adjmatrix , I )
#     
#     beta_commute0 <- wb0%*% ( XTY_tar +  Reduce(`+`, sum0) )
#     
#     
#     ##beta_commute1-K
#     beta_commute <- list()
#     
#     ## Get lambda for beta_commute0-K
#     x<-rbind(X.tar,X.src[[1]],X.src[[2]],X.src[[3]]);y<-rbind(y.tar,y.src[[1]],y.src[[2]],y.src[[3]])
#     lambda = cv.glmnet(x,y)$lambda.min;adjmatrix <- lambda*I
#     
#     for (k in 1:K) {
#       
#       # calibrated
#       sum_calibrated <- Reduce(`+`, XTY_src[1:k])+
#         Reduce(`+`, XTX_src_delta[1:k])
#       
#       if(k==K){
#         # combine
#         
#         
#         beta_commute[[k]] <- ginv(XTX_tar + Reduce(`+`, XTX_src[1:k])+adjmatrix )%*%
#           ( XTY_tar + sum_calibrated)
#       }
#       else{ 
#         
#         # syn
#         sum_syn <- Reduce(`+`, nZ_score[(k+1):K])+ # 1129 nz_score
#           Reduce(`+`, scale_XTX_tar_delta[(k+1):K]) # 1129 1129 scale_src_tar[[i]]*XTX_tar_delta
#         
#         #combine
#         A = t(X.tar)%*%X.tar+Reduce(`+`, XTX_src[1:k])+ Reduce(`+`, scale_src_tar[(k+1):K])*t(X.tar)%*%X.tar
#         
#         
#         
#         wbk<-woodbury(A,adjmatrix,I)
#         beta_commute[[k]] <- wbk%*% # ginv((1+Reduce(`+`, scale_src_tar[(k+1):K]) )*XTX_tar + Reduce(`+`, XTX_src[1:k]) )
#           (XTY_tar + sum_calibrated + sum_syn)
#       }
#       #+Reduce(`+`, scale[(k+1):K]) )*XTX_tar
#       
#     }
#     
#     
#     
#     #####################
#     ### aggregation
#     #####################
#     X.pool2 = X.pool[,1:p]
#     X.til = Data.gen.one(n.tar=100, X.pool2)$X.tar
#     # include intercepts
#     X.til = cbind(rep(1,nrow(X.til)),X.til)
#     y.til <- X.til %*% as.matrix(beta.true) + rnorm(nrow(X.til))
#     
#     
#     ### fastCOMMUTE M=0 (federated)
#     B.syn.M0 = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute0)
#     wt.syn.M0.agg.new <- Agg.fun.new(B.syn.M0, X.til, y.til, const=1)
#     beta_commute0 <- B.syn.M0%*%wt.syn.M0.agg.new
#     
#     ### COMMUTE M=1
#     B.syn.M1 = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute[[1]])
#     wt.syn.M1.agg.new <- Agg.fun.new(B.syn.M1, X.til, y.til, const=1)
#     beta_commute[[1]] <- B.syn.M1%*%wt.syn.M1.agg.new
#     
#     ### COMMUTE M=2
#     B.syn.M2 = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute[[2]])
#     wt.syn.M2.agg.new <- Agg.fun.new(B.syn.M2, X.til, y.til, const=1)
#     beta_commute[[2]] <- B.syn.M2%*%wt.syn.M2.agg.new
#     
#     ### COMMUTE M=3 (pooledTL)
#     B.TL = cbind(beta.tar, wthreshold[[1]], wthreshold[[2]], wthreshold[[3]], beta_commute[[3]])
#     wt.TL.agg.new <- Agg.fun.new(B.TL, X.til, y.til, const=1)
#     beta_commute[[3]] <- B.TL%*%wt.TL.agg.new
#     
#     
#     ### !!!put k=0 to the last one
#     beta_commute[[K+1]]<-beta_commute0
#     return (beta_commute)
#   }



data_generation <- function(sim, p, s, K, n0, nt, nov = "large",h, sig.delta, intercept, exact, r, n.test){
  
  set.seed(sim)
  n.tar <- n0 
  n.src <- nt 
  coef.all <- Coef.gen(s=s, h=h, K=K, sig.delta, p=p, exact, intercept) #exact = False?
  beta.true <- as.numeric(coef.all$beta)
  B=cbind(beta.true, coef.all$w) #for target then for source(w for source)
  

  # Choose p variables for simulation temporarily
  X.pool2 = X.pool[,1:p]
  # There are 20000 rows, just sample n.tar(100) rows for target
  sample.tar = sample(1:nrow(X.pool2), n.tar, replace = FALSE)
  # Sample out
  X.tar = X.pool2[sample.tar,] 
  # include intercepts
  X.tar = cbind(rep(1,nrow(X.tar)),X.tar)
  # Generate y.tar using MLR, without intercept, don't forget the epsilon. 
  y.tar<- X.tar %*% beta.true + rnorm(nrow(X.tar))
  
  # temp generate values for the targeted response variable
  # y.tar.fast.M3 <- X.tar%*%beta.true
  
  # renew X.pool
  X.pool2 = X.pool2[-sample.tar,] 
  X.src <- y.src <- list()
  for(k in 1:K){
    # from the rest of X.pool, get all of the source data
    sample.k = sample(1:nrow(X.pool2), n.src[k], replace = FALSE) 
    X.src[[k]] = X.pool2[sample(1:nrow(X.pool2), n.src[k], replace = FALSE),] 
    # include intercepts
    X.src[[k]] = cbind(rep(1,nrow(X.src[[k]])),X.src[[k]])
    y.src[[k]]<- X.src[[k]]%*%B[, k+1]+ rnorm(nrow(X.src[[k]]))
    X.pool = X.pool[-sample.k,]
  }
  

  
  # Test data
  X.test = X.pool2[sample(1:nrow(X.pool2), n.test, replace = FALSE),]
  # include intercepts
  X.test = cbind(rep(1,nrow(X.test)),X.test)
  y.test <- X.test%*%B[,1]+rnorm(n.test)
  
  
  ###directly use cv.glmnet to fit lasso
  beta.tar <- as.numeric(ST.init(X.tar, y.tar)$beta0)
  
  w = list()
  TL = list()
  for(k in 1:K){
    w[[k]] <- as.numeric(ST.init(X.src[[k]], y.src[[k]])$beta0)
    print(paste0('w',k))
    TL[[k]] <- TL.init(X.tar, y.tar, X.src=X.src[[k]], y.src=y.src[[k]], w=w[[k]])
    print(paste0('beta.TL',k))
  }
  
  ###calculate delta for each source
  delta.TL = list()
  for(k in 1:K){
    delta.TL[[k]] = thres(TL[[k]]$delta0, n.tar, p) ###add threshold
  }
  
  
  wthreshold = list()
  #source-only LASSO + threshold
  for (i in 1:K) {
    wthreshold[[i]] = thres(w[[i]], sqrt(n.src[i]), p)
  }

  
  return(list(X.tar=X.tar,X.src=X.src,
              y.tar=y.tar,y.src=y.src,
              delta.TL=delta.TL,beta.true=beta.true,
              beta.tar.lasso = beta.tar,
              wthreshold = wthreshold))
}