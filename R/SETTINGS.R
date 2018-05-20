#' calculation of best intial settings

SETTINGS <- function(VARIANCE,N,BETA_LAMDBA,alpha,beta){
  
    if (length(alpha!=N)) alpha = rep(1,N)
    if (length(beta!=N))  beta  = rep(1,N)
    
    R     <- rep(0,2)
    ROH   <- rep(0,2)
    
    CONTAINER  <- VARIANCE[,2]
    
    for (i in 3:length(VARIANCE[1,]-1)){
      CONTAINER= c(CONTAINER,VARIANCE[,i])
    }
    

    PHI   <- MASS::fitdistr(1/CONTAINER, "gamma")
    
    
    ALPHA <- rep(1,N)*PHI[[1]][1]*alpha
    BETA  <- rep(1,N)*(PHI[[1]][2])*BETA_LAMDBA*beta
    R[2]      <- 1000
    ROH[2]    <- 10
    R[1]      <- 1000
    ROH[1]    <- 10
    
    TAU     <- rep(1,N)
    LAMBDA1 <- rep(1,N)
    LAMBDA2 <- 1 
    
    LIST <- list(R,ROH,ALPHA,BETA,TAU,LAMBDA1,LAMBDA2)
    names(LIST) <- c("R","ROH","ALPHA","BETA","TAU","LAMBDA1","LAMBDA2")
    return(LIST)

}
