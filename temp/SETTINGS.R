SETTINGS <- function(VARIANCE,N,BETA_LAMDBA,NAME='JAKSTAT'){
  if (NAME == 'JAKSTAT') {
    R     <- rep(0,2)
    ROH   <- rep(0,2)
    PHI   <- MASS::fitdistr(1/(c(VARIANCE[["sd_STAT5p_cyt"]],VARIANCE[["sd_STAT5ptot_cyt"]])), "gamma")
    ALPHA <- rep(1,N)*PHI[[1]][1]
    BETA  <- c(1,1,1,.1)*(PHI[[1]][2])*BETA_LAMDBA
    R[2]      <- 1000
    ROH[2]    <- 10
    R[1]      <- 1000
    ROH[1]    <- 10
    
    TAU     <- c(1,1,1,1)
    LAMBDA1 <- c(1,1,1,1)
    LAMBDA2 <- 1 
    
    LIST <- list(R,ROH,ALPHA,BETA,TAU,LAMBDA1,LAMBDA2)
    names(LIST) <- c("R","ROH","ALPHA","BETA","TAU","LAMBDA1","LAMBDA2")
    return(LIST)
  }
}