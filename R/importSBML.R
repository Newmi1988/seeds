#'  Import SBML Models using the Bioconductor package 'rsbml'
#'
#'  A simple function for importing sbml models from a extensive markup language file.
#'  
#' @param  filename name of the import file. Should be located in the working directory.
#' @param  times     timestep at which the function should be evaluated
#' @param  y measurements of the model
#'  
#' @return returns a odeModel-Object
#'
#' @export importSBML
#'
importSBML <- function(filename, times, y) {

  if(!require('rsbml',character.only = TRUE)) {
    cat('Please install rsbml from the Bioconducture reposotory')
  } else {
    requireNamespace("rsbml")
    model <- rsbml::rsbml_read(filename = filename, dom = TRUE)
    
    # return(model)
    
    states <- model@model@species
    parameter <- model@model@parameters
    # print(parameter)
    
    # measurements
    rules <- model@model@rules
    reactions <- model@model@reactions
    reacList <- list()
    for (i in 1:length(reactions)) {
      reacList[[i]] <- gsub(pattern = "expression",replacement = '',deparse(model@model@reactions[[i]]@kineticLaw@math, width.cutoff = 300))
    }

    stoichM <- rsbml::stoichiometryMatrix(object = model@model)
    react <- c()
    combieReact <- function(reactStrs, stMatrix) {
      for (i in 1:nrow(stMatrix)) {
        m <- which(stMatrix[i,] !=0 )
        if(length(m)>0) {
          react <- c(react,paste0(stMatrix[i,m],'*',reactStrs[m], collapse = ' + '))
        }
      }
      return(react)
    }

    react <- combieReact(reacList,stoichM)
    
    meas <- c()
    if ( length(rules) != 0) {
      for (i in 1:length(rules)) {
        meas[i] <- gsub(pattern = "expression", replacement = '', x = rules[[i]]@math)
      }
      
    }
    
    reformatEqs <- function(reactions, states,  measureRules){
      
      
      xStates <- paste('x',1:length(states), sep = '')
      for (i in 1:length(states)) {
        regState <- paste0('\\b',states[i],'\\b')
        reactions = unlist(lapply(X = reactions,FUN = function(x) gsub(pattern = regState, replacement = xStates[i], x = x)))
        if (length(measureRules) != 0) {
          measureRules = unlist(lapply(X = measureRules,FUN = function(x) gsub(pattern = regState, replacement = xStates[i], x = x)))
        }
      }
      if ( length(measureRules) != 0 ) {
        res = list('reac' = reactions, 'meas'= measureRules) 
      } else {
        res = list('reac' = reactions)
      }
      return(res)
    }
    
    
    reactNames = rownames(stoichM[rowSums(stoichM)!=0,])
    eqList <- reformatEqs(reactions = react, states = reactNames, measureRules= meas)
    
    # format the parameter and initial vector into names vectors
    if ( length(parameter) != 0 ) {
      v <- c()
      n <- c()
      for (i in 1:length(parameter)){
        v[i] <- parameter[[i]]@value
        n[i] <- parameter[[i]]@name
      }
      namedParaVec <- v
      names(namedParaVec) <- tolower(n)
    }
    
    initToPara <- function(model,namedParaV){
      const <- c()
      for (i in 1:length(model@model@species)){
        const[i] = model@model@species[[i]]@constant
      }
      constSpec <- model@model@species[[which(const,TRUE)]]  
      
      
      nV <- c()
      namesV <- c()
      for (i in 1:length(constSpec)){
        nV[i] = constSpec@initialAmount
        namesV[i] = constSpec@id
      }
      names(nV) <- tolower(namesV)
      
      namedParaV = c(namedParaV, nV)
      return(namedParaV)
    }
   
    if ( length(parameter) != 0) {
      namedParaVec = initToPara(model, namedParaVec)
    } 
    
    initVec <- model@model@species
    print(length(initVec))
    
    if( length(initVec) != 0) {
    
      v <- c()
      n <- c()
      for (i in 1:length(initVec)){
        v[i] <- initVec[[i]]@initialAmount
      }
      initState <- v
      # names(initState) <- n
      initState = initState[rowSums(stoichM)!=0] 
      print(eqList)
      eqFuncList = writeDummy(eqList)
    } else {
      test <- model@model@reactions[[1]]
      split <- strsplit(split = 'where')
      
      print(split)
    }
      
    # print(model@model@reactions)
      
    model <- odeModel(func = eqFuncList$reac, parms = namedParaVec, measFunc = eqFuncList$meas, y = initState, times = times, meas = y)
  }
  unloadNamespace('rsbml')
  
  return(model)
}


