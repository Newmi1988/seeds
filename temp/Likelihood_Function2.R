#log[L(G|x)P(G)]
LOGLIKELIHOOD_func  <- function(pars,Step,OBSERVATIONS,x_0,parameters,EPS_inner,INPUT,D,GIBBS_PAR,k,MU_JUMP,SIGMA_JUMP,eps_new){
  

  EPS_inner=rbind(EPS_inner,eps_new)
  EPS_PARS <- EPS_inner
  
  EPS_PARS[2,k] <- pars

  
  RATIO_partial_new <- PARTIALLIKELIHOOD_func(Step,OBSERVATIONS,x_0,parameters,INPUT,EPS_PARS,GIBBS_PAR[["BETA"]],GIBBS_PAR[["ALPHA"]]) 
  
  if (!is.na(RATIO_partial_new)){
  
  r                 <- RATIO_partial_new+mvtnorm::dmvnorm(EPS_PARS[2,],EPS_inner[1,],D,log=TRUE)
  return(r)
  } 
  
  return(NA)
  
}

#partial
PARTIALLIKELIHOOD_func <- function(STEP,OBSERVATIONS,x_0,parameters,input,W,BETA,ALPHA){

  TIME <-  c(OBSERVATIONS[STEP-1,1],OBSERVATIONS[STEP,1])
  
  
  X <- ode_solv(TIME,x_0,parameters,input,W,OBSERVATIONS[(STEP-1):STEP,1])

  if(is.na(X)) return(NA)
  
  SUM <-  sum(-log((1+(1/(2*BETA))*(((OBSERVATIONS[2,-1]-sapply(1:4,objective,y=tail(X,1),parameter=parameters[5:6],USE.NAMES = TRUE))^2)))^(ALPHA+0.5)))
        

  
  return(SUM)
  
  
}
