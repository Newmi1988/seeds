#' Greedy method for estimating a sparse solution
#'
#' calculates controls based on a first optimisation with gradient descent; should result in a sparse vector
#' of hidden inputs.
#' 
#' @param odeModel a object of class \code{\link{odeModel}} from the package seeds. The class saves the details of an experiment for easier manipulation and analysis. 
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
#' @param x0             inital state of the ODE system. Can be supplied with  the odeModel class.
#'
#' @param optW           a vector that indicates for which knots of the network a input should be calculated. The default is all nodes.
#'
#' @param measFunc       a R-Function that is used for measurement of the states if the system is not completly
#'                measurable; an empty argument will result in the assumption that all states of the system are
#'                measurable. Can be supplied by the odeModel parameter.
#'
#' @param measData       a table that containts the measurements of the experiment. Used to calculate the needed inputs. Can be supplied with  the odeModel class.
#'
#' @param parameters      vector or named vector that contains the parameters of the ODE equation. Can be supplied with  the odeModel class.
#'
#' @param modelFunc      a R-Function that states the ODE system for which the hidden inputs should be calculated. Can be supplied with  the odeModel class.
#'
#' @param greedyLogical         a boolean that states if the greedy approach should be used;if set to FALSE the algorithm
#'                will only use perform a calculation of the inputs for all knots without a sparse solution
#'
#' @param plotEstimates  boolean that indicated if the current estimate should be plotted.
#'
#' @param Beta          skaling parameter for the backtracking to approximate the stepsize of the gradient descent. Is set to  0.8
#'                 if no value is given to the function
#'
#' @param sd      Standard deviation of the measurement. Is used to weight the errors of the estimates in the cost function. Optional parameter. Can be supplied with  the odeModel class.
#' 
#' @param conjGrad Boolean that indicates the usage of conjugate gradient method over the normal steepest descent. Defaults to true if not specified.
#' 
#' @param cString Optional parameter: A string that represents constants, can be used to calculate a hidden input for a component that gradient is zero.
#'
#' @param systemInput A dataset that discribes the external input of the system. The time steps should be given in the first column for the interpolation.
#' 
#' @param epsilon parameter that defines the stopping criteria for the algorithm, in this case percent change in cost function J[w]
#'
#' @param nnStates A bit vector indicating the states that should be non negative. Default behaviour will calculete positive and negative states. Can be supplied with  the odeModel class.
#'
#' @return returns a list of results objects. The default plot function can be used to plot the results.
#'
#' @export

sgdn <- function(odeModel, alphaStep, Beta, alpha1, alpha2, x0, optW, measFunc, measData, sd, epsilon,
                           parameters, systemInput, modelFunc, greedyLogical, plotEstimates, conjGrad, cString, nnStates) {

  #### new object implementation ####
  if (!missing(odeModel)) {
    modelFunc <- odeModel@func
    parameters <- odeModel@parms
    systemInput <- odeModel@input
    if (sum(colSums(systemInput)) == 0) {
      systemInput <- NULL
    }
    measFunc <- odeModel@measFunc
    x0 <- odeModel@y
    measData <- odeModel@meas
    sd <- odeModel@sd
  }

  if ("u" %in% names(parameters)) {
    stop('Variable name u is reserved for the inputs of a system')
  }

  if (missing(systemInput)) {
    systemInput <- NULL
  }

  if (missing(epsilon)) {
    epsilon <- 0.25
  }

  if (missing(sd) || nrow(sd) == 0) {
    sd <- NULL
  }

  if (missing(cString)) {
    cString <- NULL
  }

  if (missing(greedyLogical)) {
    greedyLogical <- TRUE
  }

  if (missing(plotEstimates)) {
    plotEstimates <- FALSE
  }

  if (missing(alpha1)) {
    alpha1 <- 0
  }

  if (missing(alpha2)) {
    alpha2 <- 0.01
  }

  if (missing(alphaStep)) {
    alphaStep <- 100
  }

  if (missing(Beta)) {
    Beta <- 0.8
  }

  if (missing(parameters)) {
    parameters <- c()
  }

  if (missing(nnStates)) {
    nnStates <- rep(0, length(x0))
  }


  checkSkalar <- function(argSkalar) {
    if (length(argSkalar) > 1) {
      argName <- toString(deparse(substitute(argSkalar)))
      errortext <- ' has to be a skalar not a vector'
      stop(paste0(argName, errortext))
    }
  }

  if (missing(conjGrad)) {
    conjGrad <- T
  }

  checkSkalar(alphaStep)
  checkSkalar(Beta)
  checkSkalar(alpha1)
  checkSkalar(alpha2)

  checkFunctions <- function(argFunc) {
    if (class(argFunc) != "function") {
      argName <- toString(deparse(substitute(argFunc)))
      errorText <- ' has to be a function that can be used with deSolve. Type ??deSolve for examples and documentation.'
      stop(paste0(argName, errorText))
    }
  }

  checkFunctions(modelFunc)
  checkFunctions(measFunc)

  checkLogical <- function(argLog) {
    if (!is.logical(argLog)) {
      argName <- toString(deparse(substitute(argLog)))
      errorText <- ' has to be a logical.'
      stop(paste0(argName, errorText))
    }
  }

  checkLogical(plotEstimates)
  checkLogical(greedyLogical)

  if (missing(optW)) {
    optW <- rep(1, length(x0))
  }

  checkDimensions <- function() {
    if (length(x0) != length(optW)) {
      stop('The vectors x0 and optW must have the same length')
    }
  }

  checkDimensions()


  #### Creation of C-files for use with deSolve ####
  # extract the equations of the model and save them in an odeEq object
  odeEq <- new("odeEquations")
  odeEq <- createModelEqClass(odeEq, modelFunc)
  odeEq <- setMeassureFunc(odeEq, measFunc)

  numInputs = length(x0) + 1

  # the model equations will be written wo a C file
  createCFile(parameters = parameters, inputs = numInputs, Eq = odeEq, nnStates = nnStates)

  odeEq <- isDynElaNet(odeEq)
  odeEq <- calculateCostate(odeEq)

  createFunctions(odeEq)

  # check the operating system
  #   Windows uses Rtools for compilation
  #   Unix systems should be distributed with a C compiler
  if (grepl("Rtools", Sys.getenv('PATH')) || (.Platform$OS.type != "windows")) {
    if (.Platform$OS.type != "windows") {
      cat('Using compiled code for more speed.')
    } else {
      cat('Rtools found. Using compiled code for more performance.\n')
    }
    # check the compiled file extensions
    #   Platform    filename extension
    #   Windoes     .dll
    #   Unix        .so 
    ext <- .Platform$dynlib.ext
    compiledModel <- paste0('model', ext)

    if (.Platform$OS.type != "windows"){
      temp_compiled_model <- paste0(tempdir(),'/',compiledModel)
    } else {
      temp_compiled_model <- paste0(tempdir(),'\\',compiledModel)
      temp_compiled_model = gsub('\\\\','/', temp_compiled_model)
    }
    # check if the library is loaded, so changes can be applied
    if (is.loaded('derivsc')) {
      dyn.unload(temp_compiled_model)
    }
    
    if (.Platform$OS.type != "windows"){
      temp_file_path <- paste0(tempdir(),'/','model.c')
    } else {
      temp_file_path <- paste0(tempdir(),'\\','model.c')
      temp_file_path = gsub('\\\\', '/', temp_file_path)
    }

    # compile the C function of the system
    system(paste0("R CMD SHLIB ",temp_file_path))
    # load the dynamic link library
    
    
    dyn.load(temp_compiled_model)
  } else {
    cat('No installation of Rtools detected using the normal solver.\n')
  }
  iter <- (sum(optW))

  #### initialize start of alpha2 estimation ####
  estiAlpha2 <- list()
  alpha2Start <- 1 # starting value for estimating alpha2
  steps <- 6 # number of values that are valuated for best fit
  numCores <- parallel::detectCores() - 1
  error <- matrix(rep(0, 2), ncol = 2)
  colnames(error) <- c('alpha', 'MSE')

  #### parallel estimation of a fitting alpha2 value ####
  if (is.null(alpha2) && requireNamespace('parallel', quietly = TRUE) && requireNamespace('doParallel', quietly = TRUE) && requireNamespace('foreach', quietly = TRUE) && numCores > 1) {


    # load the dynamic linked shared object library
    worker.init <- function() {
      dyn.load(temp_compiled_model)
    }

    cat('More than 1 core detected, using parallel computing.\n')
    exportVars <- c('dynElasticNet', 'y', 'costate')

    cl <- parallel::makeCluster(numCores)
    if (grepl("Rtools", Sys.getenv('PATH'))) {
      parallel::clusterCall(cl, worker.init)
    }
    doParallel::registerDoParallel(cl)

    estiAlpha2 <- foreach::foreach(i = 1:steps, .export = exportVars)
    estiAlpha2 = foreach::'%dopar%'(estiAlpha2,
        dynElasticNet(alphaStep = alphaStep, armijoBeta = Beta, x0 = x0, optW = optW, eps = epsilon,
                      measFunc = measFunc, measData = measData, SD = sd, constStr = cString,
                      alpha1 = 0, alpha2 = alpha2Start * 10 ^ (1 - i), modelInput = systemInput,
                      parameters = parameters, modelFunc = modelFunc, maxIteration = 100, plotEsti = FALSE, conjGrad = conjGrad, nnStates = nnStates)
      )
    parallel::stopCluster(cl)
    closeAllConnections()
    error[1,] = c(10 ^ (1 - 1), mean(estiAlpha2[[1]]$rmse))
    for (i in 2:steps) {
      error = rbind(error, c(alpha2Start * 10 ^ (1 - i), mean(estiAlpha2[[i]]$rmse)))
    }

  } else if (is.null(alpha2)) {

    cat('\nNo installation of package doParallel found:\n')
    cat('Using sequencial optimisation to find a fitting value of alpha2.\n')

    alpha1 = 0
    for (i in 1:steps) {

      alpha2 = alpha2Start * 10 ^ (1 - i)
      cat('\nOptimization with alpha2=', alpha2, '\n')
      estiAlpha2[[i]] <- dynElasticNet(alphaStep = alphaStep, armijoBeta = Beta, x0 = x0, optW = optW, eps = epsilon,
                                         measFunc = measFunc, measData = measData, SD = sd,
                                         alpha1 = alpha1, alpha2 = alpha2, constStr = cString,
                                         parameters = parameters, modelFunc = modelFunc, modelInput = systemInput,
                                         maxIteration = 100, plotEsti = plotEstimates, conjGrad = conjGrad, nnStates = nnStates)
      if (i == 1) {
        error[i,] = c(alpha2, mean(estiAlpha2[[i]]$rmse))
      } else {
        error = rbind(error, c(alpha2, mean(estiAlpha2[[i]]$rmse)))
      }

    }

    slopeErr <- abs(diff(error[, 1]) / diff(error[, 2]))
    #slopeErr = slopeErr[which(slopeErr >0 )]
    changeTresh <- min(which(slopeErr < 0.5))

    alpha2 = error[changeTresh + 1, 1] # alpha2 is selected based on the squared error at the given measurement times
    results <- estiAlpha2[[changeTresh + 1]] # use the estimated results of the estimation for saving time

    cat('Conservativ estimated alpha2=', alpha2)


  } else {
    results <- dynElasticNet(alphaStep = alphaStep, armijoBeta = Beta, x0 = x0, optW = optW, eps = epsilon,
                             measFunc = measFunc, measData = measData, SD = sd,
                             alpha1 = alpha1, alpha2 = alpha2, constStr = cString,
                             parameters = parameters, modelFunc = modelFunc, plotEsti = plotEstimates,
                             modelInput = systemInput, conjGrad = conjGrad, nnStates = nnStates)
  }

  resAlg <- list()
  if (!greedyLogical || (sum(optW) == 1)) {
    resAlg[[1]] <- results
    i = 2
  } else {

    orgOptW <- optW <- results$optW
    orgAUC <- results$AUC
    optWs <- list()
    costError <- cbind(rep(0, length(optW)), rep(0, length(optW)))
    colnames(costError) <- c('sum(MSE)', 'cost')


    for (i in 1:(iter - 1)) {
      cat('_________________________________________\n')
      cat('selection done: starting new optimization\n')
      cat('optimizing states:\n')
      cat(which(optW > 0))
      optWs[[i]] <- optW
      resAlg[[i]] <- dynElasticNet(alphaStep = alphaStep, armijoBeta = Beta, alpha1 = alpha1, alpha2 = alpha2, x0 = x0, optW = optW, eps = epsilon,
                                   measFunc = measFunc, measData = measData, SD = sd, modelInput = systemInput, constStr = cString,
                                   parameters = parameters, modelFunc = modelFunc, origAUC = orgAUC, plotEsti = plotEstimates, conjGrad = conjGrad, nnStates = nnStates)



      costError[i,] = c(mean(resAlg[[i]]$rmse), resAlg[[i]]$J)


      # use best fit inteads last iteration
      if (i > 1 && (costError[i, 1] > costError[i - 1, 1])) {
        cat('hidden inputs on knots:\n')
        cat(which(optWs[[i - 1]] %in% 1))
        cat('\n')
        break
      }

      if (sum(colSums(resAlg[[i]]$w[, -1])) == 0) {
        orgAUC[which(optW > 0)] = 0
        optW <- resAlg[[i]]$optW - optW
      } else {
        optW <- resAlg[[i]]$optW
      }

    }

    # unload the dynamic linked shared object library
    # has to be unleaded to makes changes
    if (grepl("Rtools", Sys.getenv('PATH')) || (.Platform$OS.type != "windows")) {
      dyn.unload(temp_compiled_model)
    }


    if ((length(resAlg) == (iter - 1))) {
      cat('Best solution for the given Problem.\n Returning solution with best fit\n')
      costError <- costError[, 2]
      costError <- costError[costError > 0]
      i <- which(costError == min(costError))
      cat('Best solution in interation: ', i)
      resAlg$optimalSol <- i
      resAlg$measurements <- measData
      i = i + 1
    } else {
      resAlg$optimalSol <- i - 1
      resAlg$measurements <- measData
    }

  }

  res <- list()
  resLength <- i - 1

  #### reformating and return####
  for (i in 1:resLength) {

    states <- as.data.frame(resAlg[[i]]$x[])
    colnames(states)[1] <- "t"
    stateUnsc <- states
    stateUnsc[, 2:ncol(stateUnsc)] = NaN



    hiddenInp <- as.data.frame(resAlg[[i]]$w)
    colnames(hiddenInp)[1] <- "t"
    hiddenInpUnsc <- hiddenInp
    hiddenInpUnsc[, 2:ncol(hiddenInpUnsc)] = NaN

    outputMeas <- as.data.frame(resAlg[[i]]$y)
    outPutUnsc <- outputMeas
    outPutUnsc[, 2:ncol(outPutUnsc)] = NaN

    if (is.null(sd)) {
      emptyStd <- matrix(rep(0, length(measData[, -1, drop = FALSE])), ncol = ncol(measData[, -1, drop = FALSE]))
      dataError <- data.frame(t = measData[, 1], emptyStd)
      colnames(dataError) <- c("t", paste0('y', 1:(ncol(emptyStd))))
    } else {
      dataError <- cbind(t = measData[, 1], sd)
      colnames(dataError) <- c("t", paste0('y', 1:(ncol(sd))))
    }

    colnames(measData) <- c("t", paste0('y', 1:(ncol(measData[, -1, drop = FALSE]))))

    nomStates <- as.data.frame(resAlg[[i]]$nomX)
    colnames(nomStates)[1] = "t"

    res[[i]] <- resultsSeeds(stateNominal = nomStates,
                      stateEstimates = states,
                      stateUnscertainLower = stateUnsc,
                      stateUnscertainUpper = stateUnsc,
                      hiddenInputEstimates = hiddenInp,
                      hiddenInputUncertainLower = hiddenInpUnsc,
                      hiddenInputUncertainUpper = hiddenInpUnsc,
                      outputEstimates = outputMeas,
                      outputEstimatesUncLower = outPutUnsc,
                      outputEstimatesUncUpper = outPutUnsc,
                      Data = measData,
                      DataError = dataError
  )

  }


  return(res)

}
