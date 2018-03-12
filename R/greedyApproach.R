#' Greedy Approach Algorithm
#'
#' calculates controls based on a first optimisation with gradient descent; should result in a sparse vector
#' of hidden inputs
#'
#' @param alphaStep      the starting stepsize for the gradient descent
#'                a fitting stepsize will be calculated based on a backtracking line search
#'                if the algorithm converges to slow use a bigger stepsize
#'
#' @param alpha1         L1-norm parameter of the dynamic elastic net approach, is set to zero for this algorithm
#'
#' @param alpha2         L2-norm parameter of the dynamic elastic net approach
#'                used for regulation purposes
#'                set to NULL for a approximation of alpha2 - will results in a longer runtime
#'
#' @param x0             inital state of the ODE system
#'
#' @param optW           a vector that indicates for which knots of the network a input should be calculated
#'
#' @param times          time sequence for which output is wanted; the first value of times must be the initial time
#'
#' @param measFunc       a R-Function that is used for measurement of the states if the system is not completly
#'                measurable; an empty argument will result in the assumption that the complete system is
#'                measurable
#'
#' @param measData       a table that containts the measurements of the experiment; used to calculate the needed inputs
#'
#' @param parameter      vector or named vector that contains the parameters of the ODE equation
#'
#' @param modelFunc      a R-Function that states the ODE system for which the hidden inputs should be calculated
#'
#' @param greedy         a boolean that states if the greedy approach should be used; if set to FALSE the algorithm
#'                will only use perform a calculation of the inputs for all knots without a sparse solution
#'
#' @param plotEstimates  boolean that indicated if the current estimate should be shown
#'
#' @return returns a results-object with default plot function. The plot shows the estimated best sparse fit
#' 
#' @examples
#' N = 10^0.31
#' x0 = c(N, 0, 0, 0)
#' y <- c(X = x0)
#' times <- c( 0.0208,  0.1098,   0.2696,    0.4999,    0.8002,    1.1697,    1.6077,    2.1129,    2.6843,    3.3205,    4.0200,    4.7811,    5.6020,    6.4808,    7.4154,    8.4035,    9.4429,   10.5310, 11.6653,   12.8431,   14.0616,   15.3179,   16.6090,   17.9319,   19.2834,   20.6603,   22.0594,   23.4773,   24.9107,   26.3561,   27.8102,   29.2695,   30.7305,   32.1898,   33.6439,   35.0893, 36.5227,   37.9406,   39.3397,   40.7166,   42.0681,   43.3910,   44.6821,   45.9384,   47.1569,   48.3347,   49.4690,   50.5571,   51.5965,   52.5846,   53.5192,   54.3980,   55.2189,   55.9800, 56.6795,   57.3157,   57.8871,   58.3923,   58.8303,   59.1998,   59.5001,   59.7304,   59.8902,   59.9792)
#' parameters = 10^c(0.31, -1, -0.49, 0.42, -0.21, -0.34)
#' 
#' inputData <- jakstatInput
#' measure <- jakstatMeasurement
#' 
#' 
#' modelJakStat  <- function(t, x, parameters, input) {
#'   with (as.list(parameters),{
#'     
#'     k1 = parameters[1]
#'     k2 = parameters[2]
#'     k3 = parameters[3]
#'     k4 = parameters[4]
#'     s1 = parameters[5]
#'     s2 = parameters[6]
#'     
#'     u <- input$u(t)
#'     
#'     dx1 = -k1 * x[1]  * u
#'     dx2 = k1 *  x[1]  * u - k2 * x[2]^2
#'     dx3 = -k3*x[3] + 0.5*k2*x[2]^2
#'     dx4 = k3 * x[3]
#'     
#'     list(c(dx1 ,dx2 ,dx3 ,dx4 ))
#'   })
#' }
#' 
#' measJakStat <- function(x) {
#'   
#'   s1 <- 10^(-0.21)
#'   s2 <- 10^(-0.34)
#'   
#'   y1 = s1*(x[,2]+ 2*x[,3])
#'   y2 = s2*(x[,1] + x[,2] + 2*x[,3])
#'   
#'   return(list(y1,y2))
#' }
#' 
#' y <- data.frame(measure$t, measure$y1, measure$y2)
#' sd <- data.frame(measure$y1sd, measure$y2sd)
#' 
#' JakStatConst <- '2*x4+ 2*x3 + x1 + x2 == N'
#' 
#' 
#' results <- greedyApproach(alphaStep = 0.01, alpha2 = 0.1,
#'                           x0 = x0, optW = c(1,1,1,1) , times=times,
#'                           measFunc= measJakStat,  measData = y, std = sd,
#'                           parameters = parameters, systemInput = inputData,
#'                           modelFunc = modelJakStat, plotEstimates = TRUE, conjGrad = FALSE, cString = JakStatConst)
#' 
greedyApproach <- function(alphaStep,Beta,alpha1, alpha2, x0, optW, times, measFunc, measData, std,
                           parameters, systemInput, modelFunc, greedyLogical, plotEstimates, conjGrad, cString) {

  if(missing(systemInput)) {
    systemInput <- NULL
  }
  
  if(missing(std)){
    std <- NULL
  }
  
  if(missing(cString)) {
    cString <- NULL
  }

  if(missing(greedyLogical)) {
    greedyLogical <- TRUE
  }

  if(missing(plotEstimates)) {
    plotEstimates <- FALSE
  }

  if(missing(alpha1)) {
    alpha1 <- 0
  }

  if(missing(Beta)) {
    Beta <- 0.8
  }

  checkSkalar <- function(argSkalar) {
    if(length(argSkalar)>1) {
      argName <- toString(deparse(substitute(argSkalar)))
      errortext <- ' has to be a skalar not a vector'
      stop(paste0(argName,errortext))
    }
  }

  if(missing(conjGrad)){
    conjGrad <- T
  }

  checkSkalar(alphaStep)
  checkSkalar(Beta)
  checkSkalar(alpha1)
  checkSkalar(alpha2)

  checkFunctions <- function(argFunc) {
    if(class(argFunc)!= "function") {
      argName <- toString(deparse(substitute(argFunc)))
      errorText <- ' has to be a function that can be used with deSolve. Type ??deSolve for examples and documentation.'
      stop(paste0(argName,errorText) )
    }
  }

  checkFunctions(modelFunc)
  checkFunctions(measFunc)

  checkLogical <- function(argLog) {
    if(!is.logical(argLog)) {
      argName <- toString(deparse(substitute(argLog)))
      errorText <- ' has to be a logical.'
      stop(paste0(argName,errorText) )
    }
  }

  checkLogical(plotEstimates)
  checkLogical(greedyLogical)


  checkDimensions <- function() {
    if(length(x0)!= length(optW)) {
      stop('The vectors x0 and optW must have the same dimensions')
    }
  }

  checkDimensions()

  #create the needed files
  odeEq <- odeEq()
  odeEq <- createModelEqClass(odeEq,modelFunc)
  odeEq <- setMeassureFunc(odeEq,measFunc)
  odeEq <- isDynElaNet(odeEq)
  odeEq <- calculateCostate(odeEq)
  createFunctions(odeEq)
  
  detach('package:Deriv')

  source('costate.R')
  source('stateHiddenInput.R')


  iter <- (sum(optW))

  estiAlpha2 <- list()

  if(is.null(alpha2) && requireNamespace('parallel', quietly = TRUE) && requireNamespace('doParallel', quietly = TRUE) && requireNamespace('foreach', quietly = TRUE)) {
    alpha2Start <- 1  # starting value for estimating alpha2

    steps <- 6        # number of values that are valuated for best fit
                      # alpha2Start * 10^(1-i)      i = 1:steps

    error <- matrix(rep(0,2),ncol=2)
    colnames(error) <- c('alpha','MSE')

    noCores <- detectCores() -1
    if(noCores > 1){
      # library('parallel')
      # library('doParallel')
      # library('foreach')
      print(costate)
      print('More than 1 core detected, using parallel computing.')
      exportVars <- c('dynElasticNet','measFunc', 'y', 'costate')

      cl <- parallel::makeCluster(noCores)
      registerDoParallel(cl)

      estiAlpha2 <- foreach::foreach(i = 1:steps, .export = exportVars) %dopar% {
        dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta,x0 = x0, optW = optW,
                      times=times, measFunc= measFunc, measData = y, STD = std, constStr = cString,
                      alpha1 = 0, alpha2 = alpha2Start*10^(1-i), modelInput = systemInput,
                      parameters = parameters, modelFunc = modelFunc,maxIteration=100, plotEsti = FALSE, conjGrad = conjGrad)
      }
      doParallel::stopCluster(cl)
      error[1,] = c( 10^(1-1),mean(estiAlpha2[[1]]$rmse))
      for( i in 2:steps) {
        error = rbind(error,c( alpha2Start*10^(1-i),mean(estiAlpha2[[i]]$rmse)))
      }
      print(error)


    } else {

      alpha1 = 0
      for (i in 1:steps) {

        alpha2 = alpha2Start*10^(1-i)
        estiAlpha2[[i]] <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta,x0 = x0, optW = optW,
                                         times=times, measFunc= measFunc, measData = y, STD = std,
                                         alpha1 = alpha1, alpha2 = alpha2, constStr = cString,
                                         parameters = parameters, modelFunc = modelFunc, modelInput = systemInput,
                                         maxIteration=100, plotEsti = plotEstimates, conjGrad = conjGrad)
        if (i==1){
          error[i,] = c(alpha2,mean(estiAlpha2[[i]]$rmse))
        } else {
          error = rbind(error,c(alpha2,mean(estiAlpha2[[i]]$rmse)))
        }
        print(error)
      }
    }

    slopeErr <- diff(error[,1]) / diff(error[,2])
    slopeErr = slopeErr[which(slopeErr >0 )]
    changeTresh <- min(which(slopeErr <0.5)) + 1

    alpha2 = alpha2Start*10^(-changeTresh)  # alpha2 is selected based on the squared error at the given measurement times
    results <- estiAlpha2[[changeTresh-1]]     # use the estimated results of the estimation for saving time


  } else {
    results <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta, x0 = x0, optW = optW,
                             times=times, measFunc= measFunc, measData = y, STD = std,
                             alpha1 = alpha1, alpha2 = alpha2, constStr = cString,
                             parameters = parameters, modelFunc = modelFunc, plotEsti = plotEstimates,
                             modelInput = systemInput, conjGrad = conjGrad)
  }

  if (!greedyLogical) {
    return(results)
  }
  else {
    orgOptW <- optW <- results$optW
    orgAUC <- results$AUC
    barplot(orgAUC)
    print(results$rmse)
    optWs <- list()
    resAlg <- list()
    costError <- cbind(rep(0,length(optW)),rep(0,length(optW)))
    colnames(costError) <- c('sum(MSE)','cost')

    # alphaStep = alphaStep*10

    for(i in 1:(iter-1)) {
      cat('_________________________________________\n')
      cat('selection done: starting new optimization\n')
      cat('selected inputs:\n')
      cat(which(optW > 0))
      optWs[[i]] <- optW
      resAlg[[i]] <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta, alpha1 = alpha1, alpha2 = alpha2,x0 = x0, optW = optW,
                                   times=times, measFunc= measFunc, measData = y, STD = std, modelInput = systemInput, constStr = cString,
                                   parameters = parameters, modelFunc = modelFunc, origAUC = orgAUC, plotEsti = plotEstimates, conjGrad = conjGrad)

      # costError[i,] = c(sum(resAlg[[i]]$rmse),resAlg[[i]]$J)
      costError[i,] = c(mean(resAlg[[i]]$rmse),resAlg[[i]]$J)


      ## use best fit inteads last iteration
      if(i > 1 && ( costError[i,1] > costError[i-1,1])  ) {
        cat('hidden inputs on knots:\n')
        cat(which(optWs[[i-1]] %in% 1))
        break
      }
      optW <- resAlg[[i]]$optW
    }

    if(length(resAlg)==(iter-1)) {
      cat('The algorithm did not find a minimal sparse solution. Returning last solution as best fit\n')
      resAlg$optimalSol <- i
      resAlg$measurements <- measData
      i = i+1
    } else {
      #return(resAlg[[i-1]])
      resAlg$optimalSol <- i-1
      resAlg$measurements <- measData
    }

    res <- results(modelFunction = odeEq@modelStr,
                   measureFunction = odeEq@measureStr,
                   hiddenInputs = resAlg[[i-1]]$w,
                   auc = resAlg[[i-1]]$AUC,
                   estimatedStates = resAlg[[i-1]]$x[],
                   estimatedMeasurements = resAlg[[i-1]]$y,
                   optimalSolution = resAlg$optimalSol,
                   allData = resAlg
                   )


    return(res)

  }

}
