### a S4 class to save the important information of a experiment

#' A class to store the important information of an experiment. 
#' 
#' The slots are used to store the important information of an experiment. The class is used to create object for the
#' two algorithms implemented in seeds. Methods are implemented to easy calculate the nominal solution of the model and
#' change the details of the saved model.
#' 
#' The numerical solutions are calculated using the \pkg{deSolve} - package. 
#' 
#' @slot func A funtion containing the ode-equations of the model. For syntax look at the given examples of the \pkg{deSolve} package.
#' @slot parms the parameters of the model
#' @slot input matrix containing the inputs with the time points
#' @slot measFunc function that converts the output of the ode solution
#' @slot y initial (state) values of the ODE system, has to be a vector
#' @slot meas matrix with the (experimental) measurements of the system
#' @slot sd optional standard deviations of the measurements, is used by the algorithms as weights in the costfunction
#' 
#' @export odeModel
#' @exportClass odeModel


odeModel <- setClass(
  #name of Class
  "odeModel",
  slots = c(
    func = "function",
    parms = "numeric",
    input = "data.frame",
    measFunc = "function",
    y = "numeric",
    meas = "data.frame",
    sd = "data.frame"
  ),
  
  prototype = list(
    func = function(x) {},
    parms = numeric(),
    input = data.frame(matrix(numeric(0), ncol = 2)),
    measFunc =  function(x) {},
    y =  numeric(0),
    meas =  data.frame(matrix(numeric(0), ncol = 2)),
    sd = data.frame(matrix(numeric(0), ncol = 2))
  ),
  validity = function(object) {
    # check inputs of matrix slots
    meas <- object@meas
    measBool <- checkMatrix(meas)
    if(!is.null(measBool)) {
      return(measBool)
    }
    
    # if(ncol(object@meas) !=  ncol(object@sd)) {
    #   return("For every measurement the standard deviation has to be given.")
    # }
    
    # if(nrow(object@meas) != 1 && nrow(sd) == 1) {
    #   object@sd = NULL
    # }
    
    if(length(object@y) != 0) {
      m <- matrix(rep(0,length(object@y)),ncol = length(object@y))
      testMeas <- do.call(cbind,object@measFunc(m))

      if(ncol(testMeas) != (ncol(object@meas)-1)){
        return("The returned results of the measurement function does not have the same
               dimensions as the given measurements")
      }
    }
    
    return(TRUE)
  }
)

setMethod('initialize', "odeModel", function(.Object,...){
  .Object <- callNextMethod()
  return(.Object)
})

checkMatrix <- function(argMatrix) {
  if(sum(argMatrix)==0) {
    argName <- toString(deparse(substitute(argMatrix)))
    errortext <- ' has to contain values not equal to 0.'
    return(paste0(argName,errortext))
  }
}


setGeneric(name="setModelEquation",
           def = function(theObject,func)
           {
             standardGeneric("setModelEquation")
           }
)

#' Set the model equation
#' 
#' Set the model equation of the system. Has to be a function that can be used with the deSolve package
#' @param theObject an object of the class modelClass
#' @param func function describing the ode equation of the model 
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setModelEquation",
          signature = "odeModel",
          definition = function(theObject,func)
          {
            theObject@func <- func
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setParms",
           def = function(theObject,parms)
           {
             standardGeneric("setParms")
           }
)

#' set the model parameters 
#' 
#' @param theObject an object of the class modelClass
#' @param parms a vector containing the parmeters of the model 
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setParms",
          signature = "odeModel",
          definition = function(theObject,parms)
          {
            theObject@parms <- parms
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setInput",
           def = function(theObject,input)
           {
             standardGeneric("setInput")
           }
)

#' Set the inputs of the model. 
#' 
#' @param theObject an object of the class modelClass
#' @param func function describing the ode equation of the model 
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setInput",
          signature = "odeModel",
          definition = function(theObject,input)
          {
            theObject@input  <- input
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setMeasFunc",
           def = function(theObject,measFunc)
           {
             standardGeneric("setMeasFunc")
           }
)


#' Set the measurement equation for the model
#' 
#' @param theObject an object of the class modelClass
#' @param measFunc measurement function of the model. Has to be a R functions.
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setMeasFunc",
          signature = "odeModel",
          definition = function(theObject,measFunc)
          {
            theObject@measFunc  <- measFunc
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setY",
           def = function(theObject,y)
           {
             standardGeneric("setY")
           }
)


#' Set the vector with the initial (state) values
#' 
#' @param theObject an object of the class modelClass
#' @param y vector with the initial values
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setY",
          signature = "odeModel",
          definition = function(theObject,y)
          {
            theObject@y  <- y
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setMeas",
           def = function(theObject,meas)
           {
             standardGeneric("setMeas")
           }
)

#' set measurements of the model
#' 
#' @param theObject an object of the class modelClass
#' @param meas measurements of the model, a matrix with measurements of the model
#' and the corresponding time values
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setMeas",
          signature = "odeModel",
          definition = function(theObject,meas)
          {
            theObject@meas <- meas
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setSd",
           def = function(theObject,sd)
           {
             standardGeneric("setSd")
           }
)

#' Set the standard deviation of the measurements
#' 
#' @param theObject an object of the class modelClass
#' @param sd a matrix with the standard deviations of the measurements
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = "setSd",
          signature = "odeModel",
          definition = function(theObject,sd)
          {
            theObject@sd <- sd
            validObject(theObject)
            return(theObject)
          }
)


setGeneric(name = 'genCCode',
           def = function(odeModel,bden,logTrans)
           {
             standardGeneric('genCCode')
           }
)

setMethod(f = 'genCCode',
          signature = c('odeModel','missing','missing'),
          definition =  function(odeModel,bden,logTrans){
            createCompModel(modelFunc = odeModel@func,parameters = odeModel@parms)
            return(odeModel)
          }
)

setMethod(f = 'genCCode',
          signature = c('odeModel','logical','missing'),
          definition =  function(odeModel,bden,logTrans){
            createCompModel(modelFunc = odeModel@func,parameters = odeModel@parms, bden = bden)
            return(odeModel)
          }
)

setMethod(f = 'genCCode',
          signature = c('odeModel','logical','numeric'),
          definition =  function(odeModel,bden,logTrans){
            createCompModel(modelFunc = odeModel@func,parameters = odeModel@parms, bden = bden, logTransfVar = logTrans)
            return(odeModel)
          }
)

setMethod(f = 'genCCode',
          signature = c('odeModel','missing','numeric'),
          definition =  function(odeModel,bden,logTrans){
            createCompModel(modelFunc = odeModel@func,parameters = odeModel@parms, logTransfVar = logTrans)
            return(odeModel)
          }
)

# nominal solution
setGeneric(name = 'nominalSol',
           def = function(odeModel, logTrans)
           {
             standardGeneric('nominalSol')
           }
)


#' Calculate the nominal solution of the model
#' 
#' @param odeModel a object of the class ode model describing the experiment
#' @param logTrans a vector indicating which of the state vector componenets should be log transformed
#' 
#' @rdname modelClass-methods
#' 
#' @export
setMethod(f = 'nominalSol',
          signature = c('odeModel','missing'),
          definition =  function(odeModel,logTrans){
            createCompModel(modelFunc = odeModel@func,parameters = odeModel@parms, logTransfVar = logTrans)
            
            times <- odeModel@input[,1]
            input <- odeModel@input
            x0 <- odeModel@y
            colnames(input) <- rep('',ncol(input))
            w <- matrix(rep(0,length(x0)*length(times)), ncol = length(x0))
            wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
            
            u <- apply(X = input[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approx(x = input[,1], y = x, xout = times, rule = 2))

            uList = list(cbind(times,u[[1]]$y))
            wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
            forcings <- c(uList, wList)

            ext <- .Platform$dynlib.ext
            compiledModel <- paste0('model',ext)
          
            if(is.loaded('derivsc')){
              dyn.unload(compiledModel)
            }

            system("R CMD SHLIB model.c")
            dyn.load(compiledModel)
            
            resOde <- deSolve::ode(y = odeModel@y, times = times, func = "derivsc",
                                       parms = odeModel@parms, dllname = "model", initforc="forcc",
                                       forcings = forcings, initfunc = "parmsc")
            
            dyn.unload(compiledModel)
            
            return(resOde)
          }
)

setMethod(f = 'nominalSol',
          signature = c('odeModel','numeric'),
          definition =  function(odeModel,logTrans){
            createCompModel(modelFunc = odeModel@func,parameters = odeModel@parms, logTransfVar = logTrans)
            
            times <- odeModel@input[,1]
            input <- odeModel@input
            colnames(input) <- rep('',ncol(input))
            w <- matrix(rep(0,length(x0)*length(times)), ncol = length(y))
            wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
            
            u <- apply(X = input[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approx(x = input[,1], y = x, xout = times, rule = 2))
            
            uList = list(cbind(times,u[[1]]$y))
            wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
            forcings <- c(uList, wList)
            
            ext <- .Platform$dynlib.ext
            compiledModel <- paste0('model',ext)
            
            if(is.loaded('derivsc')){
              dyn.unload(compiledModel)
            }
            
            if(min(logTrans) >0 && max(logTrans) <= length(odeModel@y)){
              x0[logTrans] = log(x0[logTrans])
            } else {
              stop(paste0('The given values of logTrans have to be between ',0,' and ',length(odeModel@y)))
            }
            
            system("R CMD SHLIB model.c")
            dyn.load(compiledModel)
            
            out <- deSolve::ode(y = odeModel@y, times = times, func = "derivsc",
                                       parms = odeModel@parms, dllname = "model", initforc="forcc",
                                       forcings = forcings, initfunc = "parmsc")
            
            dyn.unload(compiledModel)
            
            
            
            out[,logTrans+1] = exp(out[,logTrans+1])
            
            return(out)
          }
)








