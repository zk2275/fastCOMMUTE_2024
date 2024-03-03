# time test for fastCOMMUTE

# M=0



nZ_score <- Z_score <- GKoff_Z_score <- GKoff <- list()
XTX_tar <- t(X.tar)%*%X.tar
XTY_tar <- t(X.tar)%*%y.tar




scale_src_tar <- XTX_src <- XTY_src <- XTX_src_delta <- 
  scale_XTY_tar<- XTX_tar_delta <- scale_XTX_tar_delta<- list()
for (i in 1:K) {
  XTX_src[[i]] <- t(X.src[[i]])%*%X.src[[i]]
  XTY_src[[i]] <- t(X.src[[i]])%*%y.src[[i]]
  scale_src_tar[[i]] <- (nrow(X.src[[i]])-1)/(nrow(X.tar)-1)
  XTX_src_delta[[i]] <- t(X.src[[i]])%*%X.src[[i]]%*%delta.TL[[i]]
  XTX_tar_delta[[i]] <- t(X.tar)%*%X.tar%*%delta.TL[[i]] # may be unnecessary
  
  scale_XTX_tar_delta[[i]] <- scale_src_tar[[i]]*t(X.tar)%*%X.tar%*%delta.TL[[i]]
  scale_XTY_tar[[i]]  <- scale_src_tar[[i]]*t(X.tar)%*%y.tar
}


  
  for (i in 1:K) {
    nZ_score[[i]] <-  scale_src_tar[[i]]*t(X.tar)%*%X.tar%*%wthreshold[[i]]
    Z_score[[i]] <- 1/sqrt(nrow(X.tar))*nZ_score[[i]]
    
  }



##beta_commute M=0
sum0 <- list()
for (i in 1:K) {
  sum0[[i]] <- nZ_score[[i]] + scale_src_tar[[i]]*XTX_tar_delta[[i]] #1129 scale_src_tar[[i]]*XTX_tar_delta
}



#sum_of_list <- Reduce(`+`, my_list)
 
wb0<-woodbury(t(X.tar)%*%X.tar, Reduce(`+`, scale_src_tar)*t(X.tar), X.tar)
beta_commute0 <- wb0%*% # ginv((1+Reduce(`+`, scale_src_tar)  )*XTX_tar  )
  ( XTY_tar +  Reduce(`+`, sum0) )

start.time <- Sys.time() 
k=3
##beta_commute1-K
beta_commute <- list()

  
  # calibrated
  sum_calibrated <- Reduce(`+`, XTY_src[1:k])+
    Reduce(`+`, XTX_src_delta[1:k])


    
  beta_commute[[k]] <- ginv(XTX_tar + Reduce(`+`, XTX_src[1:k]) )%*%
    ( XTY_tar + sum_calibrated)

  #+Reduce(`+`, scale[(k+1):K]) )*XTX_tar
  



end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

### !!!put M=0 to the last one
beta_commute[[K+1]]<-beta_commute0

