#' estimating the optimal control using the dynamic elastic net
#'
#' @param alphaStep starting value of the stepsize for the gradient descent, will be calculate to minimize the costfunction by backtracking algorithm
#' @param armijoBeta  scaling of the alphaStep to find a approximatly optimal value for the stepsize
#' @param x0 initial state of the ode system
#' @param parameters parameters of the ODE-system
#' @param times measurement times at which the estimations will be evaluated
#' @param alpha1 L1 cost term skalar
#' @param alpha2 L2 cost term skalar
#' @param measData measured values of the experiment
#' @param STD standard deviation of the experiment; leave empty if unknown
#' @param modelFunc function that describes the ODE-system of the model
#' @param measFunc function that maps the states to the outputs
#' @param optW vector that indicated at which knots of the network the alogrithm should estimate the hidden inputs
#' @param origAUC AUCs of the first optimisation; only used by the algorithm
#' @param plotEsti boolean that controlls of the current estimates should be plotted
#' @param modelInput an dataset that discribes the external input of the system
#' @param conjGrad boolean that indicates the usage of conjugate gradient method over the normal steepest descent
#' @param constStr  a string that represents constrains, can be used to calculate a hidden input for a komponent that gradient is zero
#' @param maxIteration a upper bound for the maximal number of iterations
#' @param eps citeria for stopping the algorithm
#'
#' @return A list containing the estimated hidden inputs, the AUCs, the estimated states and resulting measurements and the costfunction
#' @export
dynElasticNet <- function(alphaStep,armijoBeta,x0,parameters,times,alpha1,alpha2,measData, constStr,
                          STD,modelFunc,measFunc,modelInput,optW,origAUC,maxIteration,plotEsti, conjGrad, eps) {

  source('stateHiddenInput.R')
  source('costate.R')
  
  costate <- get('costate', envir = environment())
  hiddenInputState <- get('hiddenInputState', envir = environment())
  

  startOptW <- optW

  # Check for the optW vector, only states with an index of 1 will be calclated to enable the 'greedy' approach
  if (missing(optW)) {
    optW <- rep(1,length(x0))
  }

  if (missing(plotEsti)) {
    plotEsti <- FALSE
  }





  # Initial Options
  optInit <- sum(optW)  # Anzahl der zu optimierenden Inputs

  N <- 100
  t0 <- times[1]
  tf <- utils::tail(times, n=1)
  times <- seq(from = t0, to = tf, length.out = N)
  tInt <- c(t0,tf)

  measureTimes <- measData[,1] # select the time points of the given data
  measureData <- measData[,-1] # select the data

  # setting alpha1 = 0 for this approach
  alphaDynNet <- list(a1 = alpha1, a2 = alpha2) # list of the alpha_1 and alpha_2 values

  # interpolation of the standard deviations
  if(is.null(STD)) {
    # if the standard deviations are missing the weights are set to 1 for every point of measurement, no weights will be applied
    measureData <- as.matrix(measureData)
    Q <- matrix(data = rep(1,length(measureData)), ncol = ncol(measureData))
    ## normalization of the weights based on the meassured data as suggested by Dominik Kahl
    Q <- Q / abs(measureData)
    Q[is.infinite(Q)] = 0
    Q[is.na(Q)] = 0
    interpQ <- apply(X = Q, MARGIN = 2, FUN = function(t) stats::approx(x=measureTimes, y=t, xout = times))
    interpQ = do.call(cbind, lapply(interpQ, FUN = function(t) cbind(t$y)) )
    Q <- interpQ
  }
  else {
    # interpolating the given standard deviations to be used as weights in the costete equation
    interpSTD <- apply(X = STD, MARGIN = 2, FUN = function(t) stats::approx(x = measData[,1], y = t, xout = times))
    interpSTD = do.call(cbind, lapply(interpSTD, FUN = function(t) cbind(t$y)))
    Q <- apply(X = interpSTD, MARGIN = 2, FUN = function(t)  (1/t^2)/length(t) )
  }

  # get the trajectories of the nominal mondel without inputs
  # using deSolve 'ode' function
  # the first evaluation of the costfunction is based on the nominal model with hidden inputs constant at zero for all times
  if(all(is.null(names(x0)))) {
    names(x0) <- paste0(rep("x",length(x0)),1:length(x0))
  }
  
  if(!is.null(modelInput)){
    inputInterp <- list()
    inputInterp <- apply(X = modelInput[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approxfun(x = modelInput[,1], y = x, rule = 2, method = 'linear'))
    solNominal <- as.data.frame(deSolve::ode(y = x0, times = times, func = modelFunc, parms = parameters, input = inputInterp))  
  } else {
    solNominal <- as.data.frame(deSolve::ode(y = x0, times = times, func = modelFunc, parms = parameters))
  }
  
  
  ### with c code #########################################################
  if(!is.null(modelInput)){
    inputApprox <- apply(X = modelInput[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approx(x = modelInput[,1], y = x, xout = times, rule = 2))
    inputApprox = list(cbind(times,inputApprox$u$y))
  } else {
    inputApprox <- list(cbind(times,rep(0,length(times))))
  }
  
  if(grepl("Rtools",Sys.getenv('PATH'))){
    w <- matrix(rep(0,length(x0)*length(times)), ncol = length(x0))
    wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
    wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
    forcings <- c(inputApprox, wList)
    
    solNominal = deSolve::ode(y = x0, times, func = "derivsc",
                              parms = parameters, dllname = "model", initforc="forcc",
                              forcings = forcings, initfunc = "parmsc")
  }

  
  

  # Tx the timepoints of the solution of the nominal model
  Tx <- solNominal[,1]
  # X trajectories of the states
  x <- solNominal[,-1, drop=FALSE]

  w = matrix(rep(0,nrow(x)*ncol(x)), nrow = nrow(x)) # is initialized as constant zero
  colnames(w) <- paste0(rep("w",ncol(w)),1:ncol(w))

  getMeassures <- function(x,measFunc) {
    if(missing(measFunc)) {
      # if there is no function given for measurements it is assumed that all states are measurable
      cat('No meassurement function defined. Assuming all states are observable.\n')
      y <- x[,-1, drop=FALSE]
    } else {
      # combine the listed measurements to a matrix using do.call and cbind
      y <- do.call(cbind,measFunc(x[,-1,drop=FALSE]))
    }
    y = as.data.frame(cbind(x[,1],y))
    names(y)[1] <- 't'
    names(y)[-1] <- paste0(rep("y",ncol(y)-1),1:(ncol(y)-1))
    return(y)
  }

  # function to calculate the stepsize for the gradient descent
  getAlphaBacktracking <- function(oldW,W,Q,y,gradStep,J,currIter,alphaDynNet,alphaS,stepBeta,optW,para,tInt,Tp,measFunc,input,measureTimes) {


    iter = 200
    alpha = alphaS
    arrayJ = rep(0,iter)
    
    time <- seq(from = tInt[1], to = tInt[2], length.out = 100)
    solX <- matrix( rep(0,length(time)*(length(optW)+1)),ncol = length(optW)+1)
    Tx <- rep(0,length(time))
    x <- matrix(0,length(optW)*length(time), ncol = length(optW))
    yHat <- matrix(rep(0,length(measData)), ncol = ncol(measData))
    
   
    for (i in 1:iter) {
      newW = oldW + alpha*gradStep
      
      if(grepl("Rtools",Sys.getenv('PATH'))){
        wSplit <- split(newW, rep(1:ncol(newW), each = nrow(newW)))
        wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
        forcings <- c(inputApprox, wList)
        solX = deSolve::ode(y = x0, time, func = "derivsc",
                            parms = parameters, dllname = "model", initforc="forcc",
                            forcings = forcings, initfunc = "parmsc")
      } else {
        input$optW = optW
        input$w = apply(X = newW, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tp, y = x, method = 'linear', rule=2))
        solX = deSolve::ode(y = x0, times = time,func = hiddenInputState, parms = parameters, input=input)
      }


      
      Tx = solX[,1]
      x = solX[,-1, drop=FALSE]

      yHat = getMeassures(solX,measFunc)
      input$interpX = apply(X = x, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tx, y = x, rule=2, method = 'linear'))
      input$interpyHat = apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) stats::approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
      
      arrayJ[i] = costFunction(measureTimes,input,alphaDynNet)
      
      # cat(paste0('i=',i,' J[w]=',arrayJ[i], ' alpha=', alpha,' Beta=',stepBeta ,'\n'))

      if ( i>1 && (arrayJ[i]>arrayJ[i-1]) && (arrayJ[i] < J[currIter])) {
        alpha = alphaS*stepBeta^(i-2)
        break
      }
      # beta = stepBeta^(i)
      alpha = alpha*stepBeta
    }
    # quadratic interpolation to find the minimum
    
    intAlpha1 <- alpha
    intAlpha2 <- alpha*stepBeta
    costAlpha1 <- arrayJ[i-1]
    costAlpha2 <- arrayJ[i]
  
    
    cubicInterpolMin <- function(alphaA,alphaB,jA,jB){
      alpha3 <- 0.5*(alphaA+alphaB)
      newW = oldW + alpha3*gradStep
      
      if(grepl("Rtools",Sys.getenv('PATH'))){
        wSplit <- split(newW, rep(1:ncol(newW), each = nrow(newW)))
        wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
        forcings <- c(inputApprox, wList)
        solX = deSolve::ode(y = x0, time, func = "derivsc",
                            parms = parameters, dllname = "model", initforc="forcc",
                            forcings = forcings, initfunc = "parmsc")
      } else {
        input$optW <- optW
        input$w <- apply(X = newW, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tp, y = x, method = 'linear', rule=2))
        time <- seq(from = tInt[1], to = tInt[2], length.out = 300)
        solX <- deSolve::ode(y = x0, times = time,func = hiddenInputState, parms = parameters, input=input)
      }

      
      
      
      Tx <- solX[,1]
      x <- solX[,-1, drop=FALSE]
      
      yHat <- getMeassures(solX,measFunc)
      
      input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tx, y = x, rule=2, method = 'linear'))
      input$interpyHat <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) stats::approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
      
      j3 = costFunction(measureTimes,input,alphaDynNet)

      alphaT = alpha3 - (alphaB-alphaA)/4 * (jB - jA)/(jB-2*j3+jA)
      
      if(is.nan(alphaT)){
        return(alphaA)
      } else {
        if(alphaT > 0 ){
          alpha = alphaT
          return(alpha)
        } else {
          alpha <- cubicInterpolMin(alphaA = alphaA, alphaB = alpha3, jA = jA, jB = j3)
          return(alpha)
        }
      }
    }
    
    # return(alpha)

    #check if the cubicInterpolation gives a lower value as the last iteration
    # 
    alphaTemp <- cubicInterpolMin(alphaA = intAlpha1, alphaB = intAlpha2, jA = costAlpha1, jB = costAlpha2)
    newW = oldW + alphaTemp*gradStep
    
    if(grepl("Rtools",Sys.getenv('PATH'))){
      wSplit <- split(newW, rep(1:ncol(newW), each = nrow(newW)))
      wList <- lapply(wSplit, FUN = function(x) cbind(time,x))
      forcings <- c(inputApprox, wList)
      solX = deSolve::ode(y = x0, time, func = "derivsc",
                          parms = parameters, dllname = "model", initforc="forcc",
                          forcings = forcings, initfunc = "parmsc")
    } else {
      input$optW <- optW
      input$w <- apply(X = newW, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tp, y = x, method = 'linear', rule=2))
      time <- seq(from = tInt[1], to = tInt[2], length.out = 300)
      solX <- deSolve::ode(y = x0, times = time,func = hiddenInputState, parms = parameters, input=input)
    }

    Tx <- solX[,1]
    x <- solX[,-1, drop=FALSE]

    yHat <- getMeassures(solX,measFunc)

    input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tx, y = x, rule=2, method = 'linear'))
    input$interpyHat <- apply(X = yHat[,-1], MARGIN = 2, FUN = function(x) stats::approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))

    alphaCubicCOst = costFunction(measureTimes,input,alphaDynNet)

    if(alphaCubicCOst > arrayJ[i-1]){
      return(alpha)
    } else {
      return(alphaTemp)
    }
    
    
  }

  showEstimates <- function(measureTimes,AUCs,input, alpha2, J, nomSol, STD){
    tPlot <- seq(from=measureTimes[1], to = measureTimes[length(measureTimes)], length.out = 50)

    y <- sapply(input$interpY, mapply, measureTimes)
    yhat <- sapply(input$interpyHat, mapply, tPlot)
    w <- sapply(input$w, mapply, tPlot)
    yNom <- sapply(nomSol, mapply, tPlot)

    J <- unlist(J)
    J = J[J!=0]
    
    width = 2
    numMeas <- ncol(y)
    if((numMeas+3)%% 3 == 0){
      n <- (numMeas + 3) %/% 3
    } else {
      n <- (numMeas + 3) %/% 3 + 1
    }
    m <- 3
    par(mfrow=c(n,m), ask=F)
    barplot(unlist(AUCs[1,]), col = 'red', xlab = 'hidden inputs', main = 'AUC (a.u)')
    for( i in 1:numMeas) {
      yLab <- paste0('y',as.character(i))
      yMax <- max(max(y[,i]),max(yhat[,i]),max(yNom[,i]))
      yMin <- min(min(y[,i]),min(yhat[,i]),min(yNom[,i]))
      if(is.null(STD)){
        plot(x = measureTimes, y = y[,i], type = 'p', pch = 20, col = 'black', xlab = 't', ylab = yLab, ylim = c(yMin, yMax), lwd = width)
      } else {
        Hmisc::errbar(x = measureTimes, y = y[,i], yplus = y[,i]+STD[,i], yminus = y[,i]-STD[,i], ylab = yLab, ylim = c(yMin, yMax), add = FALSE)
      }
      par(new=T)
      plot(x = tPlot, y = yhat[,i], type='l', col = 'red', xlab = 't', ylab = yLab, ylim = c(yMin, yMax), lwd = width)
      par(new=T)
      plot(x = tPlot, y = yNom[,i], type='l', col = 'blue', xlab = 't', ylab = yLab, ylim = c(yMin, yMax), lwd = width)
    }
    plot(J, type = 'l', xlab = 'iteration', ylab = 'J[w]', lwd = width)
    matplot(x = tPlot, y = w, type='l', col = 'red', lwd = width)
  }

  createConst <- function(constString,needGrad) {
    trim <- function(x) gsub(pattern = '\\s', replacement = "", x = x)
    cont = strsplit(x = trim(constString), split = "==")[[1]][2]
    
    eqs <- character(length = length(needGrad))
    for(i in 1:length(needGrad)) {
      str <- paste0('Solve({',trim(constString),'},{x',needGrad[i],'})\n')
      eqs[i] <- as.character(suppressWarnings(Ryacas::yacas(str))) #warnings are generated because of an error in the orphaned xml1 package / issu is known
    }
    eq <- trim(gsub(pattern = 'list\\(||\\)\\)', replacement = "", x = eqs))
    eq = gsub(pattern = '==', replacement = '=', x = eq)
    eq = gsub(pattern = "(x)([0-9])", replacement = 'P[,\\2]', x = eq)
    eq = gsub(pattern = cont, replacement = '', x = eq)
    return(eq)
  }
  
  evalGrad <- function(constStr,gradM, optW) {
    gradzero <- which(colSums(gradM) == 0)
    optCur <- which(optW > 0)
    
    nG <- optCur[optCur %in% gradzero]
    
    if(length(nG)>0) {
      cP <- createConst(constString = constStr, needGrad = nG)
    }
    return(cP)
  }

  # cost function that is to be optimized
  costFunction <- function(measureTimes,input,alphaDynNet) {
    y <- sapply(input$interpY, mapply, measureTimes)
    yhat <- sapply(input$interpyHat, mapply, measureTimes)
    q <- sapply(input$q, mapply, measureTimes)
    w <- sapply(input$w, mapply, measureTimes)

    yCost <- list()
    # cost of the deviation of the calculated measurements to the given data
    for (i in 1:ncol(yhat)) {
      yCost$Start[[i]] = sum((yhat[1,i]- y[1,i]) * q[1,i] * (yhat[1,i]- y[1,i]))
      yCost$Middle[[i]] = sum((yhat[,i]- y[,i]) * q[,i] * (yhat[,i]- y[,i]))
      yCost$End[[i]] = sum((yhat[nrow(yhat),i]- y[nrow(y),i]) * q[nrow(q),i] * (yhat[nrow(yhat),i]- y[nrow(y),i]))
    }

    # cost for the inputs
    wCost <- list(L1 = 0, L2 = 0)
    for (i in 1:ncol(w)) {
      wCost$L1 = wCost$L1 + sum(abs(w[,i]))
      wCost$L2 = wCost$L2 + sum(abs(w[,i]^2))
    }

    #combining the costs
    cost = sum(yCost$Start)  + sum(yCost$Middle) + sum(yCost$End) + alphaDynNet$a1*wCost$L1 + alphaDynNet$a2*wCost$L2
    return(cost)
  }

  yHat <- getMeassures(solNominal,measFunc)
  yNominal <- apply(X = yHat[,-1, drop=FALSE], MARGIN = 2, FUN = function(x) stats::approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
  
  
  # interpolation
  # linear approximation of the calculated values of x,y and yhat
  xInterp <- apply(X = x, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tx, y = x, rule=2, method = 'linear'))
  yInterp <- apply(X = measData[,-1, drop=FALSE], MARGIN = 2, FUN = function(x) stats::approxfun(x = measData[,1], y = x, rule=2, method = 'linear'))
  yHatInterp <- apply(X = yHat[,-1, drop=FALSE], MARGIN = 2, FUN = function(x) stats::approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
  qInterp <- apply(X = Q, MARGIN = 2, FUN = function(x) stats::approxfun(x =times, y = x, rule = 2, method = 'linear'))
  wInterp <- apply(X = w, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tx, y = x, rule = 2, method = 'linear'))


  # list of functions that approximate the data
  if(is.null(modelInput)) {
    input <- list(optW= optW,interpX =xInterp,interpY = yInterp, interpyHat= yHatInterp, q=qInterp, w = wInterp)
  } else {
    uInterp <- inputInterp
    input <- list(optW= optW,interpX =xInterp,interpY = yInterp, interpyHat= yHatInterp, q=qInterp, w = wInterp, u = uInterp)
  }

  offset <- 5
  usedAlphas <- rep(0,offset)
  cP <- NULL

  if(missing(maxIteration)) {
    maxIter <- 100
  } else {
    maxIter <- maxIteration
  }
  
  J <- rep(0,maxIter)
  J[1] = costFunction(measureTimes,input,alphaDynNet)
  cat('\n')
  cat(paste0('nominal cost J[w]= ',J[1],'\n'))
  
  #initialize objects (testing for speed)
  lT <- rep(0,ncol(x))
  timesCostate <- seq(from = tf, to= t0, length.out =N)
  solCostate <- matrix(rep(0,(length(optW)+1)*length(timesCostate)),ncol = length(optW)+1)
  Tp <- rep(0,length(timesCostate))
  P <- matrix(rep(0,length(timesCostate)*length(optW)), ncol=length(optW))
  oldW <- matrix(rep(0,length(optW)*length(timesCostate)),ncol = length(optW))
  inputState <- list()
  inputState$optW <- optW
  solX <- matrix(rep(0,length(optW)*length(times)))
  alphaS = alphaStep
  
  wApprox <- matrix(rep(0,length(times)*2), ncol=2)
  wApprox[,1] = times
  
  for (i in 1:maxIter) {
    solCostate = deSolve::ode(y = lT, times = timesCostate, func = costate, parms = parameters, input=input)
    solCostate = solCostate[nrow(solCostate):1,]

    Tp = solCostate[,1]
    P = solCostate[,-1, drop=FALSE]
    
    if(!is.null(constStr) && i == 1) {
      cP = evalGrad(gradM = P, optW = optW, constStr = constStr)
    }
    if(!is.null(cP)){
      eval(parse(text = cP))
    }

    oldW = w

    if(conjGrad){
      if(i==1){
        gNeg = P - alpha2*w
        oldGrad = -gNeg
        step = gNeg
      }
      else {
        gNeg = P - alpha2*w
        
        # newInt <- apply(X = -gNeg, MARGIN = 2, FUN = function(x) pracma::trapz(Tp, x^2))
        newGrad <- gNeg * (gNeg + oldGrad)
        newInt <- apply(X = newGrad, MARGIN = 2, FUN = function(x) pracma::trapz(Tp, x))

        oldInt <- apply(X = oldGrad, MARGIN = 2, FUN = function(x) pracma::trapz(Tp, x^2))

        newInt[is.nan(newInt)] <- 0
        oldInt[is.nan(oldInt)] <- 0
        betaTest <- sum(newInt)/sum(oldInt)
        step = gNeg + betaTest*step

        oldGrad = -gNeg

      }
    } else {
      step = P - alpha2*w
    }



    alpha = getAlphaBacktracking(oldW = oldW,W = w,Q = Q,y = measData,
                                 gradStep = step,J = J,currIter = i,alphaDynNet = alphaDynNet,
                                 alphaS = alphaS,stepBeta = armijoBeta,optW = optW,para = parameters,
                                 tInt = tInt,Tp = Tp,measFunc = measFunc,input = input,measureTimes = measureTimes)
    usedAlphas[pracma::mod(i,offset)+1] = alpha

    # calculate the new hidden inputs
    w = oldW + alpha*step

    # CALCULATION OF THE TRAJEKTORIES THAT RESULTS FROM THE NEW HIDDEN INPUT
    inputState$wInterp <- apply(X = w, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tp, y = x, method = 'linear', rule=2))
    
    if(grepl("Rtools",Sys.getenv('PATH'))){
    ### c solver
    wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
    wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
    
    forcings <- c(inputApprox, wList)
    solX = deSolve::ode(y = x0, times, func = "derivsc",
                              parms = parameters, dllname = "model", initforc="forcc",
                              forcings = forcings, initfunc = "parmsc")
    } else {
      ### normal solver
      if(!is.null(modelInput)){
        inputState$u <- inputInterp
      }
      solX = deSolve::ode(y = x0, times = times,func = hiddenInputState, parms = parameters, input=inputState)
    }
    ####
    Tx = solX[,1]
    x = solX[,-1, drop=FALSE]


    yHat <- getMeassures(solX,measFunc)
    # interp
    input$interpX <- apply(X = x, MARGIN = 2, FUN = function(x) stats::approxfun(x = Tx, y = x, rule=2, method = 'linear'))
    input$interpyHat <- apply(X = yHat[,-1, drop=FALSE], MARGIN = 2, FUN = function(x) stats::approxfun(x = yHat[,1], y = x, rule=2, method = 'linear'))
    input$w <- inputState$wInterp

    # calculate the new cosT
    J[i+1] = costFunction(measureTimes,input,alphaDynNet)

    tAUC <- measureTimes
    absW <- abs(sapply(input$w, mapply, tAUC))
    interpAbsW <- apply(X = absW, MARGIN = 2, FUN = function(x) stats::approxfun(x = tAUC, y = x, rule=2, method = 'linear'))

    AUCs <- sapply(X = interpAbsW, FUN = function(x) pracma::trapzfun(f = x, a = t0, b = tf))
    cat(paste0('Iteration ',i,' J[w]=',round(J[i+1],2),'     change J[w]: ',round((1-abs(J[i+1]/J[i]))*100,4),' % \t\talpha=',alpha,'\n'))


    if(plotEsti == TRUE) {
      showEstimates(measureTimes,AUCs,input,alpha2,J, yNominal,STD)
    }
    # if the change in the cost function is smaller that epsilon the algorithmus stops
    if (( abs(J[i+1]/J[i]) > 1-(eps/100)) && i>1) {
      break
    }
  }


  if (missing(origAUC)) {
    origAUC <- AUCs
  }

  greedySelection <- function(AUC, optW,origAUC) {
        orderAUCs <- order(-do.call(cbind,as.list(origAUC[1,])))
        tempOptW = rep(0,ncol(AUC))

        if(sum(unlist(origAUC[1,]))==sum(unlist(AUC[1,]))) {
          # if no aucs are given (first optimisation is running) -> select bigges AUCs
          tempOptW[orderAUCs[1]] = 1
        }
        else {
          tempOptW[orderAUCs[1:(sum(optW)+1)]] = 1
        }
        return(tempOptW)
  }

  rmse <- function(measureTimes,input){
    y <- sapply(input$interpY, mapply, measureTimes)
    yhat <- sapply(input$interpyHat, mapply, measureTimes)

    # feature scaling of the data
    y = apply(y, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    yhat = apply(yhat, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

    y[is.nan(y)] <- 0
    yhat[is.nan(yhat)] <- 0

    return((colSums((yhat - y)^2))/nrow(yhat))
  }

  colnames(yHat) <- append('t', paste0('y',1:(ncol(yHat)-1)))

  results <- list()
  results$w <- cbind(Tp,w)
  results$AUC <- do.call(cbind,AUCs[1,])
  results$optW <- greedySelection(AUCs, optW, origAUC)
  results$x <- solX
  results$y <- yHat
  results$rmse <- rmse(measureTimes,input)
  # lastJ <- unlist(J)
  # lastJ = lastJ[lastJ>0]
  # lastJ = lastJ[length(lastJ)]
  lastJ = J[J>0]
  lastJ = lastJ[length(lastJ)]
  results$J <- lastJ
  results$totalJ <- J
  
  cat(paste0('RMSE:',mean(results$rmse),'\n'))


  return(results)

}
