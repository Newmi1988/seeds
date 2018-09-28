#' Wrapper to the ODE Solver


ode_solv <- function(TIME,x_0,parameter,input,w_estimate,LogTransform){
  
  times = TIME
  
  if(LogTransform){
    x_0[x_0!=0] <- log(x_0[x_0!=0])
    x_0[x_0==0] <- log(0.00001)
    }
  print(x_0)
  
  if(!is.null(input)){
    inputApprox <- apply(X = input[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approx(x = input[,1], y = x, xout = times, rule = 2))
    
    inputApprox = list(cbind(times,inputApprox$u$y))
  } else {
    inputApprox <- list(cbind(times,rep(0,length(times))))
  }
  

  newW <- as.matrix(w_estimate[2,],1)

  # erweitere den Vektor um diese Werte, damit diese bei t0 zugewiesen werden
  parametersW = c("t0"= TIME[1],"w1t0" = w_estimate[1,1], 
                  "w2t0" = w_estimate[1,2], 
                  "w3t0" = w_estimate[1,3], 
                  "w4t0" = w_estimate[1,4])
  w <- matrix(rep(0,length(x_0)*length(times)), ncol = length(x_0))


  # konstante hidden inputs über die Zeit
  w[,1] = w_estimate[2,1]
  w[,2] = w_estimate[2,2]
  w[,3] = w_estimate[2,3]
  w[,4] = w_estimate[2,4]

  # aufteilung der hidden input Matrix in einzelene vectoren; Spalten sind die Werte von w_i zu bestimmten Zeitpunkten
  wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
  
  wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
  
  
  forcings <- c(inputApprox, wList)
  

  
  parameters = c(parameter, parametersW)
  
  runSilent <- function() {
    options(warn = -1)
    on.exit(options(warn = 0))
    
    
    capture.output(sol <- deSolve::ode(y = x_0, time=times, func = "derivsc",
                                       parms = parameters, dllname = "model", initforc="forcc",
                                       forcings = forcings, initfunc = "parmsc",atol = c(1e-6,1e-16,1e-6,1e-16)))
    sol
  }
  

  sol <- runSilent()
  
  if(LogTransform){
    
    sol[sol!=0]        <- exp(sol[sol!=0])
    sol[sol<=0.00001] <- 0
    
    if (any(is.na(sol))|is.null(sol)|((sum(sol< 0)!=0))){
      
      print('Log transformation failed :: No success :: skip sample')
      
      return(NA)}
    
    }
  else{

  if (!is.null(sol)){
    sol[sol> -0.1&sol<0] <-0
  } 
  
  if (any(is.na(sol))|is.null(sol)|((sum(sol< 0)!=0))){
    
    runSilent <- function() {
      options(warn = -1)
      on.exit(options(warn = 0))
      
      
      capture.output(sol <- deSolve::ode(y = x_0, time=times, func = "derivsc",
                                         parms = parameters, dllname = "model", initforc="forcc",
                                         forcings = forcings, initfunc = "parmsc",method = 'euler',atol = c(1e-6,1e-16,1e-6,1e-16)))
      sol
    }
    
    
    sol <- runSilent()
    #print('Bad parametrization -> using Euler method')

    
    }
  


  if (!is.null(sol)){
    sol[sol> -0.1&sol<0] <-0
  } 
  

  
  if (any(is.na(sol))|is.null(sol)|((sum(sol< 0)!=0))){
    
    runSilent <- function() {
      options(warn = -1)
      on.exit(options(warn = 0))
      
      
      capture.output(sol <- deSolve::ode(y = x_0, time=times, func = "derivsc",
                                         parms = parameters, dllname = "model", initforc="forcc",
                                         forcings = forcings, initfunc = "parmsc",method = 'radau'))
      sol
    }
    
    
    sol <- runSilent()
    #print('Bad parametrization -> using Euler method')

    
  }

  if (!is.null(sol)){
    sol[sol<0] <-0
  } 
  
  

  if (any(is.na(sol))|is.null(sol)|((sum(sol< 0)!=0))){
    
    print('Euler method failed :: No success :: skip sample')
    
return(NA)}
  
  }
  
  
  return(as.data.frame(sol[,1:length(x_0)+1]))
  
  
  
}

