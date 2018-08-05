### a S4 class to save the important information of a experiment
odeModel <- setClass(
  #name of Class
  "odeModel",
  slots = c(
    func = "function",
    parms = "numeric",
    input = "matrix",
    measFunc = "function",
    y = "numeric",
    meas = "matrix",
    sd = "matrix"
  ),
  
  prototype = list(
    func = function(x) {},
    parms = numeric(),
    input = matrix(numeric(0), ncol = 2),
    measFunc =  function(x) {},
    y =  numeric(0),
    meas =  matrix(numeric(0), ncol = 2),
    sd = matrix(numeric(0), ncol = 2)
  ),
  validity = function(object) {
    # check inputs of matrix slots
    meas <- object@meas
    measBool <- checkMatrix(meas)
    if(!is.null(measBool)) {
      return(measBool)
    }
    
    if(ncol(object@meas) !=  ncol(object@sd)) {
      return("For every measurement the standard deviation has to be given.")
    }
    
    # if(nrow(object@meas) != 1 && nrow(sd) == 1) {
    #   object@sd = NULL
    # }
    
    if(length(object@y) != 0) {
      m <- matrix(rep(0,length(object@y)),ncol = length(object@y))
      testMeas <- object@measFunc(m)
      if(ncol(testMeas) != ncol(object@meas)){
        return("The returned results of the measurement function does not have the same
               dimensions as the given measurements")
      }
    }
    
    return(TRUE)
  }
)

checkMatrix <- function(argMatrix) {
  if(sum(argMatrix)==0) {
    argName <- toString(deparse(substitute(argMatrix)))
    errortext <- ' has to contain values not equal to 0.'
    return(paste0(argName,errortext))
  }
}

#set function 
setGeneric(name="setModelEquation",
           def = function(theObject,func)
           {
             standardGeneric("setModelEquation")
           }
)

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

setMethod(f = "setSd",
          signature = "odeModel",
          definition = function(theObject,sd)
          {
            theObject@sd <- sd
            validObject(theObject)
            return(theObject)
          }
)









