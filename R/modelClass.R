### a S4 class to save the important information of a experiment

#' A class to store the important information of an experiment. 
#' 
#' The slots are used to store the important information of an experiment. The class is used to create object for the
#' two algorithms implemented in seeds. Methods are implemented to easily calculate the nominal solution of the model and
#' change the details of the saved model.
#' 
#' The numerical solutions are calculated using the \pkg{deSolve} - package. 
#' 
#' @slot func A funtion containing the ode-equations of the model. For syntax look at the given examples of the \pkg{deSolve} package.
#' @slot times timesteps on which the model should be evaluated
#' @slot parms the parameters of the model
#' @slot input matrix containing the inputs with the time points
#' @slot measFunc function that converts the output of the ode solution
#' @slot y initial (state) values of the ODE system, has to be a vector
#' @slot meas matrix with the (experimental) measurements of the system
#' @slot sd optional standard deviations of the measurements, is used by the algorithms as weights in the costfunction
#' @slot custom customized link function
#' @slot nnStates bit vector that indicates if states should be observed by the root function
#' @slot nnTollerance tollerance at which a function is seen as zero
#' @slot resetValue value a state should be set to by an event
#' 
#' @export odeModel
#' @exportClass odeModel
#'
#' @import methods
odeModel <- setClass(
#name of Class
  "odeModel",
  slots = c(
    func = "function",
    times = "numeric",
    parms = "numeric",
    input = "data.frame",
    measFunc = "function",
    y = "numeric",
    meas = "data.frame",
    sd = "data.frame",
    custom = 'logical',
    nnStates = 'numeric',
    nnTollerance = 'numeric',
    resetValue = "numeric"
  ),

  prototype = list(
    func = function(x) { },
    times = numeric(),
    parms = numeric(),
    input = data.frame(matrix(numeric(0), ncol = 0)),
    measFunc = function(x) { },
    y = numeric(0),
    meas = data.frame(matrix(numeric(0), ncol = 0)),
    sd = data.frame(matrix(numeric(0), ncol = 0)),
    custom = FALSE,
    nnStates = numeric(),
    nnTollerance = numeric(),
    resetValue = numeric()
  ),
  validity = function(object) {
    # check inputs of matrix slot

    if (sum(object@times) == 0) {
      return("You have to specify the times on which the equation should be evaluated. A solution can only be calculated if the a intervall or specific timesteps are given. Set the 'times'' parameter.")
    }

    if (length(object@y) != 0 && object@custom == FALSE && sum(colSums(object@meas)) != 0) {
      m <- matrix(rep(0, length(object@y)), ncol = length(object@y))
      if (is.null(object@measFunc(m)) == FALSE) {
        testMeas <- object@measFunc(m)
        if (ncol(testMeas) != (ncol(object@meas) - 1)) {
          return("The returned results of the measurement function does not have the same
                   dimensions as the given measurements")
        }
      }
    }

    return(TRUE)
  }
)

setMethod('initialize', "odeModel", function(.Object, ...) {
  .Object <- callNextMethod()
  return(.Object)
})

checkMatrix <- function(argMatrix) {
  if (sum(argMatrix) == 0) {
    argName <- toString(deparse(substitute(argMatrix)))
    errortext <- ' has to contain values not equal to 0.'
    return(paste0(argName, errortext))
  }
}




#' Set the model equation
#' 
#' Set the model equation of the system. Has to be a function that can be used with the deSolve package
#' 
#' @param theObject an object of the class odeModel
#' @param func function describing the ode equation of the model 
#' 
#' @export
setGeneric(name = "setModelEquation",
           def = function(theObject, func) {
            standardGeneric("setModelEquation")
           }
)

#' @rdname setModelEquation
setMethod(f = "setModelEquation",
          signature = "odeModel",
          definition = function(theObject, func) {
            theObject@func <- func
            validObject(theObject)
            return(theObject)
          }
)


#' Set the model parameters
#' 
#'  a method to set the model parameters of an odeModel object. 
#' 
#' @param theObject an object of the class odeModel
#' @param parms a vector containing the parmeters of the model 
#' 
#' @export
setGeneric(name = "setParms",
           def = function(theObject, parms) {
            standardGeneric("setParms")
           }
)

#' @rdname setParms
setMethod(f = "setParms",
          signature = c("odeModel", 'numeric'),
          definition = function(theObject, parms) {
            theObject@parms <- parms
            validObject(theObject)
            return(theObject)
          }
)

#' Set the inputs of the model. 
#' 
#' @param theObject an object of the class modelClass
#' @param input function describing the ode equation of the model 
#' 
#' @export
setGeneric(name = "setInput",
           def = function(theObject, input) {
            standardGeneric("setInput")
           }
)

#' @rdname setInput 
setMethod(f = "setInput",
          signature = "odeModel",
          definition = function(theObject, input) {
            theObject@input <- input
            validObject(theObject)
            return(theObject)
          }
)




#' Set the measurement equation for the model
#' 
#' @param theObject an object of the class odeModel
#' @param measFunc measurement function of the model. Has to be a R functions.
#' @param costum costum indexing for the measurement function
#' 
#' @export
setGeneric(name = "setMeasFunc",
           def = function(theObject, measFunc, costum) {
            standardGeneric("setMeasFunc")
           }
)

#' @rdname setMeasFunc
setMethod(f = "setMeasFunc",
          signature = c('odeModel', 'function', 'missing'),
          definition = function(theObject, measFunc, costum) {
            theObject@measFunc <- measFunc
            validObject(theObject)
            return(theObject)
          }
)


#' @rdname setMeasFunc
setMethod(f = "setMeasFunc",
          signature = c('odeModel', 'function', 'logical'),
          definition = function(theObject, measFunc, costum) {
            theObject@meas <- measFunc
            theObject@custom <- costum
            validObject(theObject)
            return(theObject)
          }

)

#' Set the vector with the initial (state) values
#' 
#' @param theObject an object of the class odeModel
#' @param y vector with the initial values
#' 
#' @export
setGeneric(name = "setY",
           def = function(theObject, y) {
            standardGeneric("setY")
           }
)

#' @rdname setY
setMethod(f = "setY",
          signature = "odeModel",
          definition = function(theObject, y) {
            theObject@y <- y
            validObject(theObject)
            return(theObject)
          }
)


#' set measurements of the model
#' 
#' @param theObject an object of the class odeModel
#' @param meas measurements of the model, a matrix with measurements of the model
#' and the corresponding time values
#' 
#' @export
setGeneric(name = "setMeas",
           def = function(theObject, meas) {
            standardGeneric("setMeas")
           }
)

#' @rdname setMeas
setMethod(f = "setMeas",
          signature = 'odeModel',
          definition = function(theObject, meas) {
            theObject@meas <- meas
            validObject(theObject)
            return(theObject)
          }
)

#' Set the standard deviation of the measurements
#' 
#' @param theObject an object of the class odeModel
#' @param sd a matrix with the standard deviations of the measurements
#' 
#' @export
setGeneric(name = "setSd",
           def = function(theObject, sd) {
            standardGeneric("setSd")
           }
)

#' @rdname setSd
setMethod(f = "setSd",
          signature = "odeModel",
          definition = function(theObject, sd) {
            theObject@sd <- sd
            validObject(theObject)
            return(theObject)
          }
)


#' 
#'
#'
#'
#'


setGeneric(name = 'genCCode',
           def = function(odeModel, bden, nnStates) {
            standardGeneric('genCCode')
           }
)


setMethod(f = 'genCCode',
          signature = c('odeModel', 'logical', 'missing'),
          definition = function(odeModel, bden, nnStates) {
            createCompModel(modelFunc = odeModel@func, parameters = odeModel@parms, bden = bden)
            return(odeModel)
          }
)

setMethod(f = 'genCCode',
          signature = c('odeModel', 'logical', 'numeric'),
          definition = function(odeModel, bden, nnStates) {
            createCompModel(modelFunc = odeModel@func, parameters = odeModel@parms, bden = bden, nnStates = nnStates)
            return(odeModel)
          }
)

setMethod(f = 'genCCode',
          signature = c('odeModel', 'missing', 'numeric'),
          definition = function(odeModel, bden, nnStates) {
            createCompModel(modelFunc = odeModel@func, parameters = odeModel@parms, nnStates = nnStates)
            return(odeModel)
          }
)

# nominal solution



#' Calculate the nominal solution of the model
#' 
#' @param odeModel a object of the class ode model describing the experiment
#' 
#' @export
setGeneric(name = 'nominalSol',
           def = function(odeModel) {
            standardGeneric('nominalSol')
           }
)

#' @rdname nominalSol
setMethod(f = 'nominalSol',
          signature = c('odeModel'),
          definition = function(odeModel) {

            x0 <- odeModel@y
            ### get the times from the measurements
            # add case for missing input

            times <- odeModel@times
            if (sum(colSums(odeModel@input)) == 0) {
              input <- rep(0, length(times))
              uList = list(cbind(times, input))
            } else {
              input <- odeModel@input
              u <- apply(X = input[, -1, drop = F], MARGIN = 2, FUN = function(x) stats::approx(x = input[, 1], y = x, xout = times, rule = 2))
              uList = list(cbind(times, u[[1]]$y))
            }

            w <- matrix(rep(0, length(x0) * length(times)), ncol = length(x0))


            if (grepl("Rtools", Sys.getenv('PATH')) || (.Platform$OS.type != "windows")) {


              ext <- .Platform$dynlib.ext
              compiledModel <- paste0('model', ext)
              
              if (.Platform$OS.type != "windows"){
                temp_compiled_model <- paste0(tempdir(),'/',compiledModel)
              } else { 
                temp_compiled_model <- paste0(tempdir(),'\\',compiledModel)
                temp_compiled_model = gsub('\\\\','/', temp_compiled_model)
              }
              if (is.loaded('derivsc')) {
                dyn.unload(temp_compiled_model)
              }

              createCompModel(modelFunc = odeModel@func, parameters = odeModel@parms, nnStates = odeModel@nnStates)
              
              
              
              if (.Platform$OS.type != "windows"){
                temp_file_path <- paste0(tempdir(),'/','model.c')
              } else {
                temp_file_path <- paste0(tempdir(),'\\','model.c')
                temp_file_path = gsub('\\\\', '/', temp_file_path)
              }
              
              # compile the C function of the system
              system(paste0("R CMD SHLIB ",temp_file_path))
              
              
              # system("R CMD SHLIB model.c")
              dyn.load(temp_compiled_model)


              wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
              wList <- lapply(wSplit, FUN = function(x) cbind(times, x))
              forcings <- c(uList, wList)

              if (sum(odeModel@nnStates) == 0) {

                resOde <- deSolve::ode(y = odeModel@y, times = times, func = "derivsc",
                                         parms = odeModel@parms, dllname = "model", initforc = "forcc",
                                         forcings = forcings, initfunc = "parmsc")

              } else {

                # myEvent <- eval(parse(text = createEvent(rootStates = odeModel@nnStates, tollerance = 0., value = 0.0001)))
                # RootFunc <- eval(parse(text = createRoot(rootStates = nnStates)))
                # EventFunc <- eval(parse(text = createEvent(tollerance = eventTol, value = resetValue)))
                
                eventTol <- 0.0
                resetValue <- 0.0001

                myRoot <- eval(parse(text = createRoot(rootStates = nnStates)))
                myEvent <- eval(parse(text = createEvent(tollerance = eventTol, value = resetValue)))

                resOde <- deSolve::lsoda(y = odeModel@y, times = times, func = "derivsc",
                                         parms = odeModel@parms, dllname = "model", initforc = "forcc",
                                         forcings = forcings, initfunc = "parmsc", nroot = sum(odeModel@nnStates),
                                         rootfunc = "myroot", events = list(func = myEvent, root = TRUE))

              }


              dyn.unload(temp_compiled_model)

            } else {


              odeEq <- isDynElaNet(odeModel)
              odeEq <- calculateCostate(odeModel)
              createFunctions(odeModel)
              
              
              if (.Platform$OS.type != "windows"){
                temp_hidden_input_path <- paste0(tempdir(),'/','stateHiddenInput.R')
              } else {
                temp_hidden_input_path <- paste0(tempdir(),'\\','stateHiddenInput.R')
              }

              source(temp_hidden_input_path)

              hiddenInputState <- get('hiddenInputState', envir = environment())

              input$w = apply(X = w, MARGIN = 2, FUN = function(x) stats::approxfun(x = times, y = x, method = 'linear', rule = 2))
              input$u = apply(X = input, MARGIN = 2, FUN = function(x) stats::approxfun(x = times, y = x, method = 'linear', rule = 2))

              if (sum(odeModel@nnStates) == 0) {

                resOde <- deSolve::ode(y = odeModel@y,
                                       func = hiddenInputs,
                                       times = times,
                                       parms = odeModel@parms,
                                       input = input)



              } else {
                
                eventTol <- 0.0
                resetValue <- 0.0001

                # myRoot <- eval(parse(text = createRoot(rootStates = odeModel@nnStates)))
                # myEvent <- eval(parse(text = createEvent(rootStates = odeModel@nnStates, odeModel@nnTollerance)))
                myRoot <- eval(parse(text = createRoot(rootStates = nnStates)))
                myEvent <- eval(parse(text = createEvent(tollerance = eventTol, value = resetValue)))

                resOde <- deSolve::ode(y = odeModel@y,
                                        times = time,
                                        func = hiddenInputState,
                                        parms = odeModel@params,
                                        input = input,
                                        events = list(func = myEvent, root = TRUE),
                                        rootfun = myRoot)
              }


            }



            return(resOde)
          }
)








