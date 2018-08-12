### a S4 class to save the important information of a experiment
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

setMethod(f = 'nominalSol',
          signature = c('odeModel','missing'),
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

            system("R CMD SHLIB model.c")
            dyn.load(compiledModel)
            
            solJakStat <- deSolve::ode(y = odeModel@y, times = times, func = "derivsc",
                                       parms = odeModel@parms, dllname = "model", initforc="forcc",
                                       forcings = forcings, initfunc = "parmsc")
            
            dyn.unload(compiledModel)
            
            return(solJakStat)
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








