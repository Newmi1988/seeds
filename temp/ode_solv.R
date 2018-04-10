ode_solv <- function(TIME,x_0,parameter,input,w_estimate,STEP){




  times = TIME
  
  if(!is.null(input)){
    inputApprox <- apply(X = input[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approx(x = input[,1], y = x, xout = times, rule = 2))

    inputApprox = list(cbind(times,inputApprox$EpoRp$y))
  } else {
    inputApprox <- list(cbind(times,rep(0,length(times))))
  }
  
  newW <- w_estimate
  
  wSplit <- split(newW, rep(1:ncol(newW), each = nrow(newW)))
  wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
  forcings <- c(inputApprox, wList)
  

  
  sol = deSolve::ode(y = x_0, time=times, func = "derivsc",
                      parms = parameter, dllname = "model", initforc="forcc",
                      forcings = forcings, initfunc = "parmsc")
  
  
  
   
   
    if (!is.null(sol)) sol[sol> -0.00001&sol<0] <-0
    
    
   if (is.null(sol)|((sum(sol< 0)!=0))){
    return(NA)}
    
    return(as.data.frame(sol[,1:length(x_0)+1]))
 

  
}

