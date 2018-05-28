#' calculates the Gibbs Update

GIBBS_update  <- function(D,EPS_inner,R,ROH,SGIMA_0,n,SIGMA,LAMBDA2,LAMBDA1,TAU){
  SIGMA_A       <- (n/2)
  SIGMA_B       <- 1/((SGIMA_0*0.5)+0.5*((EPS_inner[2,]-EPS_inner[1,])*(D+diag(1,n))^-1*(EPS_inner[2,]-EPS_inner[1,])))
  SIGMA   <- 1/rgamma(1, shape = SIGMA_A, scale = SIGMA_B)
  LAMBDA2_A    <- n/2+R[2]
  LAMBDA2_B    <- 1/((2*SIGMA)^-1*sum((EPS_inner[2,]-EPS_inner[1,])^2)+ROH[2])
  LAMBDA2     <- rgamma(1, shape = LAMBDA2_A, scale = LAMBDA2_B)
  for (j in 1:n){
    if ((EPS_inner[2,j]-EPS_inner[1,j]) < 0.001){
      EPS_diff = 0.001
    } else{
      EPS_diff <- EPS_inner[2,j]-EPS_inner[1,j]
    }
    TAU_A  <- sqrt((LAMBDA1[j]*SIGMA)/EPS_diff^2)
    TAU_B  <- LAMBDA1[j] 
    TAU[j] <- statmod::rinvgauss(1,TAU_A,TAU_B);                            
    if ((TAU[j] == 0)||is.na(TAU[j])){
      TAU[j] <- 1;  
    }
  }
  for (j in 1:n){
    LAMBDA1_A  <- 1+R[1]
    LAMBDA1_B  <- 1/(0.5*(TAU[j]^-1)+ROH[1])  
    LAMBDA1[j] <- rgamma(1, shape = LAMBDA1_A, scale = LAMBDA1_B)
  }       
  LIST <- list(SIGMA,LAMBDA1,LAMBDA2,TAU)
  names(LIST) <- c("SIGMA","LAMBDA1","LAMBDA2","TAU")
  return(LIST)
}
