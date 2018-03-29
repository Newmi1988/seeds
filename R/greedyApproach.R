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
#' @param parameters      vector or named vector that contains the parameters of the ODE equation
#'
#' @param modelFunc      a R-Function that states the ODE system for which the hidden inputs should be calculated
#'
#' @param greedyLogical         a boolean that states if the greedy approach should be used; if set to FALSE the algorithm
#'                will only use perform a calculation of the inputs for all knots without a sparse solution
#'
#' @param plotEstimates  boolean that indicated if the current estimate should be shown
#'
#' @param Beta          skaling parameter for the backtracking to approximate the stepsize of the gradient descent. Is set to  0.8
#'                 if no value is given to the function
#'
#' @param std     standard deviation of the measurement; used to weight the errors of the estimates in the cost function
#' 
#' @param conjGrad boolean that indicates the usage of conjugate gradient method over the normal steepest descent
#' 
#' @param cString  a string that represents constrains, can be used to calculate a hidden input for a komponent that gradient is zero
#'
#' @param systemInput an dataset that discribes the external input of the system
#' 
#' @param epsilon parameter that defines the stopping criteria for the algorithm, in this case percent change in cost function J[w]
#'
#' @return returns a results-object with default plot function. The plot shows the estimated best sparse fit
#'
#' @example examples/uvb.R
#' 
#' @export

greedyApproach <- function(alphaStep,Beta,alpha1, alpha2, x0, optW, times, measFunc, measData, std, epsilon,
                           parameters, systemInput, modelFunc, greedyLogical, plotEstimates, conjGrad, cString) {

  if(missing(systemInput)) {
    systemInput <- NULL
  }
  
  if(missing(epsilon)) {
    epsilon <- 0.25
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
  
  if(is.loaded('derivsc')){
    dyn.unload("model.dll")
  }

  # create the needed files
  odeEq <- new("odeEquations")
  odeEq <- createModelEqClass(odeEq,modelFunc)
  odeEq <- setMeassureFunc(odeEq,measFunc)

  
  numInputs = length(x0)+1
  createCFile(parameters = parameters,inputs = numInputs, odeEq)
  
  odeEq <- isDynElaNet(odeEq)
  odeEq <- calculateCostate(odeEq)
  createFunctions(odeEq)
  if(grepl("Rtools",Sys.getenv('PATH'))){
    cat('Rtools found. Using compiled code for more performance.\n')
    system("R CMD SHLIB model.c")
    # check system format for dynamic library
    ext <- .Platform$dynlib.ext
    compiledModel <- paste0('model',ext)
    dyn.load(compiledModel)
  } else {
    cat('No installation of Rtools detected using the normal solver.\n')
  }
  iter <- (sum(optW))
  estiAlpha2 <- list()

  if(is.null(alpha2) && requireNamespace('parallel', quietly = TRUE) && requireNamespace('doParallel', quietly = TRUE) && requireNamespace('foreach', quietly = TRUE)) {
    alpha2Start <- 1  # starting value for estimating alpha2

    steps <- 6        # number of values that are valuated for best fit
                      # alpha2Start * 10^(1-i)      i = 1:steps

    error <- matrix(rep(0,2),ncol=2)
    colnames(error) <- c('alpha','MSE')

    noCores <- parallel::detectCores() -1
    if(noCores > 1){

      cat('More than 1 core detected, using parallel computing.\n')
      exportVars <- c('dynElasticNet','measFunc', 'y', 'costate')

      cl <- parallel::makeCluster(noCores)
      doParallel::registerDoParallel(cl)

      estiAlpha2 <- foreach::foreach(i = 1:steps, .export = exportVars) %dopar% {
        dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta,x0 = x0, optW = optW, eps = epsilon,
                      times=times, measFunc= measFunc, measData = measData, STD = std, constStr = cString,
                      alpha1 = 0, alpha2 = alpha2Start*10^(1-i), modelInput = systemInput,
                      parameters = parameters, modelFunc = modelFunc,maxIteration=100, plotEsti = FALSE, conjGrad = conjGrad)
      }
      parallel::stopCluster(cl)
      error[1,] = c( 10^(1-1),mean(estiAlpha2[[1]]$rmse))
      for( i in 2:steps) {
        error = rbind(error,c( alpha2Start*10^(1-i),mean(estiAlpha2[[i]]$rmse)))
      }
      print(error)


    } else {

      alpha1 = 0
      for (i in 1:steps) {

        alpha2 = alpha2Start*10^(1-i)
        estiAlpha2[[i]] <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta,x0 = x0, optW = optW, eps = epsilon,
                                         times=times, measFunc= measFunc, measData = measData, STD = std,
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
    results <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta, x0 = x0, optW = optW, eps = epsilon,
                             times=times, measFunc= measFunc, measData = measData, STD = std,
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
    optWs <- list()
    resAlg <- list()
    costError <- cbind(rep(0,length(optW)),rep(0,length(optW)))
    colnames(costError) <- c('sum(MSE)','cost')

    # alphaStep = alphaStep*4

    for(i in 1:(iter-1)) {
      cat('_________________________________________\n')
      cat('selection done: starting new optimization\n')
      cat('optimizing states:\n')
      cat(which(optW > 0))
      optWs[[i]] <- optW
      resAlg[[i]] <- dynElasticNet(alphaStep = alphaStep,armijoBeta = Beta, alpha1 = alpha1, alpha2 = alpha2,x0 = x0, optW = optW, eps=epsilon,
                                   times=times, measFunc= measFunc, measData = measData, STD = std, modelInput = systemInput, constStr = cString,
                                   parameters = parameters, modelFunc = modelFunc, origAUC = orgAUC, plotEsti = plotEstimates, conjGrad = conjGrad)

      # costError[i,] = c(sum(resAlg[[i]]$rmse),resAlg[[i]]$J)
      costError[i,] = c(mean(resAlg[[i]]$rmse),resAlg[[i]]$J)


      ## use best fit inteads last iteration
      if(i > 1 && ( costError[i,1] > costError[i-1,1])  ) {
        cat('hidden inputs on knots:\n')
        cat(which(optWs[[i-1]] %in% 1))
        cat('\n')
        break
      }
      optW <- resAlg[[i]]$optW
    }
    
    dyn.unload("model.dll")

    if(length(resAlg)==(iter-1)) {
      cat('The algorithm did stop at the last combination of hidden inputs. Returning last solution as best fit\n')
      resAlg$optimalSol <- i
      resAlg$measurements <- measData
      i = i+1
    } else {
      resAlg$optimalSol <- i-1
      resAlg$measurements <- measData
    }
    res <- new('results',modelFunction = odeEq@modelStr,
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
