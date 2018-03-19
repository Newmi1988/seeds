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
#' @examples{
#' \dontrun{# usb network example
#' uvbParameter = c(  ks1=0.23,
#'                    ks2=4.0526,
#'                    kdr1=0.1,
#'                    kdr2=0.2118,
#'                    k1=0.0043,
#'                    k2=161.62,
#'                    ka1=0.0372,
#'                    ka2=0.0611,
#'                    ka3=4.7207,
#'                    kd1=94.3524,
#'                    kd2=50.6973,
#'                    kd3=0.5508,
#'                    ks3=0.4397,
#'                    kdr3=1.246,
#'                    uv=1,
#'                    ka4=10.1285,
#'                    kd4=1.1999,
#'                    n1=3,
#'                    n2=2,
#'                    n3=3.5,
#'                    kdr3a=0.9735,
#'                    kdr3b=0.406,
#'                    ksr=0.7537,
#'                    FHY3_s=5)
#' 
#' x0 = c(0.2,10,2,0,0,20,0,0,0,4.2,0.25,20,0)
#' 
#' uvbModel <- function(t,x,parameters) {
#'   ks1 = parameters[1]
#'   ks2 = parameters[2]
#'   kdr1 = parameters[3]
#'   kdr2 = parameters[4]
#'   k1 = parameters[5]
#'   k2 = parameters[6]
#'   ka1 = parameters[7]
#'   ka2 = parameters[8]
#'   ka3 = parameters[9]
#'   kd1 = parameters[10]
#'   kd2 = parameters[11]
#'   kd3 = parameters[12]
#'   ks3 = parameters[13]
#'   kdr3 = parameters[14]
#'   uv = parameters[15]
#'   ka4 = parameters[16]
#'   kd4 = parameters[17]
#'   n1 = parameters[18]
#'   n2 = parameters[19]
#'   n3 = parameters[20]
#'   kdr3a = parameters[21]
#'   kdr3b = parameters[22]
#'   ksr = parameters[23]
#'   fhy3_s = parameters[24]
#'   
#'   dx1 = ((-2) * ((ka1 * (x[1]^2) * (x[4]^2)) - (kd1 * x[5])) +
#'      (-2) * ((ka2 * (x[1]^2) * x[2]) - (kd2 * x[3])) + 
#'      ((ks1 *((1) + (uv * n3 * (x[11] + fhy3_s))))  - (kdr1 * ((1) 
#'      + (n1 * uv)) * x[1])))
#'   dx2 = ((-1) * ((ka2*(x[1]^2) * x[2]) - (kd2 * x[3])) +
#'      (-1) * ((ka4 * x[2] * x[12]) - (kd4 * x[13])))
#'   dx3 = (((ka2 * (x[1]^2) * x[2]) - (kd2*  x[3]))) 
#'   dx4 = ((-2) * (k1*(x[4]^2)) + (2) * (k2 * x[6]) + 
#'      (-2) * ((ka1 * (x[1]^2)* (x[4]^2)) - (kd1 * x[5])) 
#'      + (-1)* (ka3 * x[4] *x[7]))
#'   dx5 =  (((ka1 * (x[1]^2) * (x[4]^2)) -(kd1 * x[5])))
#'   dx6 = ((-1) * (k2 * x[6]) +  (k1 * (x[4]^2)) +(kd3 * (x[8]^2)))
#'   dx7 = ((-1) * (ka3 * x[4] * x[7]) + ((ks2 * ((1) + 
#'      (uv * x[5]))) -(kdr2 * x[7])) + (2) * (kd3 * (x[8]^2)))
#'   dx8 = ((-2) * (kd3 * x[8]^2) + (ka3 * x[4] * x[7])) 
#'   dx9  = 0 
#'   dx10 = 0
#'   dx11 =  (((ks3 * ((1) + (n2 * uv))) -
#'      (kdr3 * (((x[3] / (kdr3a + x[3])) +
#'       (x[13] / (kdr3b + x[13]))) -(x[5] / (ksr + x[5]))) *  x[11])))
#'   dx12 = ((-1) * (ka4 * x[2] * x[12]) + (kd4 * x[13]))
#'   dx13 =((ka4 * x[2] * x[12]) - (kd4 * x[13]))
#'   
#'   return(list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10,dx11,dx12,dx13)))
#' }
#' 
#' 
#' uvbMeasure <- function(x) {
#'   
#'   y1 = 2*x[,5] + x[,4] + x[,8]
#'   y2 = 2*x[,5] + 2* x[,3] + x[,1]
#'   y3 = x[,6]
#'   y4 = x[,11]
#'   y5 = x[,4]
#'   
#'   return(list(y1,y2,y3,y4,y5))
#' }
#' 
#' 
#' y <- uvbData[,1:6]
#' t <- uvbData$t
#' sd <- uvbData[,7:11]
#' 
#' res <- greedyApproach(alphaStep = 100, alpha2 = 0.0002, optW = rep(1,13), x0 = x0,
#'                       measFunc = uvbMeasure,times = t, measData = y, 
#'                      parameters = uvbParameter, modelFunc = uvbModel, 
#'                      plotEstimates = TRUE, conjGrad = TRUE)
#' }
#' }
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

  # create the needed files
  odeEq <- new("odeEquations")
  odeEq <- createModelEqClass(odeEq,modelFunc)
  odeEq <- setMeassureFunc(odeEq,measFunc)
  odeEq <- isDynElaNet(odeEq)
  odeEq <- calculateCostate(odeEq)
  createFunctions(odeEq)
  

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
        break
      }
      optW <- resAlg[[i]]$optW
    }

    if(length(resAlg)==(iter-1)) {
      cat('The algorithm did stop at the last combination of hidden inputs. Returning last solution as best fit\n')
      resAlg$optimalSol <- i
      resAlg$measurements <- measData
      i = i+1
    } else {
      #return(resAlg[[i-1]])
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
