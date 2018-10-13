importSBML <- function(modelStr) {
  if(!require('rsbml',character.only = TRUE)) {
    cat('Please install rsbml from the Bioconducture reposotory')
  } else {
    model <- rsbml_read(filename = modelStr, dom = TRUE)
    
    states <- model@model@species
    parameter <- model@model@parameters
    
    # measurements
    rules <- model@model@rules
    reactions <- model@model@reactions
     
    reacList <- list()
    for (i in 1:length(reactions)) {
      reacList[[i]] <- gsub(pattern = "expression",replacement = '',deparse(model@model@reactions[[i]]@kineticLaw@math, width.cutoff = 300))
    }

    stoichM <- stoichiometryMatrix(object = model@model)
    
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
    for (i in 1:length(rules)) {
      meas[i] <- gsub(pattern = "expression", replacement = '', x = rules[[i]]@math)
    }
    
    reformatEqs <- function(reactions, states,  measureRules){
      
      
      xStates <- paste('x',1:length(states), sep = '')
      for (i in 1:length(states)) {
        regState <- paste0('\\b',states[i],'\\b')
        reactions = unlist(lapply(X = reactions,FUN = function(x) gsub(pattern = regState, replacement = xStates[i], x = x)))
        measureRules = unlist(lapply(X = measureRules,FUN = function(x) gsub(pattern = regState, replacement = xStates[i], x = x)))
      }
      res = list('reac' = reactions, 'meas'= measureRules) 
      return(res)
    }
    
    
    reactNames = rownames(stoichM[rowSums(stoichM)!=0,])
    
    eqList <- reformatEqs(reactions = react, states = reactNames, measureRules= meas)
    
    # format the parameter and initial vector into names vectors
    v <- c()
    n <- c()
    for (i in 1:length(parameter)){
      v[i] <- parameter[[i]]@value
      n[i] <- parameter[[i]]@name
    }
    namedParaVec <- v
    names(namedParaVec) <- tolower(n)
    
    initVec <- model@model@species
    v <- c()
    n <- c()
    for (i in 1:length(initVec)){
      v[i] <- initVec[[i]]@initialAmount
      n[i] <- initVec[[i]]@name
    }
    initState <- v
    names(initState) <- n
    initState = initState[rowSums(stoichM)!=0] 
    eqFuncList = writeDummy(eqList)
    
    model <- odeModel(func = eqFuncList$reac, parms = namedParaVec, measFunc = eqFuncList$meas, y = initState)
  }
  
  return(model)
}


