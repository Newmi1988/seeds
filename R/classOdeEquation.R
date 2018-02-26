#' An S4 class used to handle formating ODE-Equation and calculate the needed functions for the SEEDS-algorithm
#'
#' @slot modelStr a vector of strings describing the ODE
#' @slot measureStr a vector of strings describing the
#' @slot origEq a vector of strings containing the original model function
#' @slot measureFunction a vector of strings containing the original measurement function
#' @slot costateEq a vector of strings describing the costate equation
#' @slot JhT a matrix of strings describing the jacobian matrix of the measurement function
#' @slot jacobian a matrix of strings representing the jacobian matrix model equations
#' @slot costFunction a string containing the cost function
#' @slot hamiltonian a string representing the Hamilton function of the model


odeEq <- setClass(
  #name of Class
  "odeEquations",
  slots = c(
    modelStr = "character",
    measureStr = "character",
    origEq = "character",
    measureFunction = "character",
    costateEq = "character",
    invSimMatrix = "matrix",
    JhT = "matrix",
    jacobian = "matrix",
    costFunction = "character",
    hamiltonian = "character",
    dynamicElasticNet = "logical"
  ),

  prototype = list(
    modelStr = character(0),
    measureStr = character(0),
    origEq = character(0),
    measureFunction =  character(0),
    costateEq = character(0),
    costFunction = character(0),
    jacobian = matrix(list(),nrow = 2, ncol = 2),
    hamiltonian = character(0),
    dynamicElasticNet = FALSE
  ),

  validity = function(object) {
    if(!(is.character(object@origEq) && is.character(object@costateEq) && is.character(object@hamiltonian) && is.matrix(object@jacobian))) {
      return("The Equations have to be characters and the Jacobian a matrix.")
    }
    return(TRUE)
  }
)

setGeneric(name="setCostateEq",
           def = function(theObject,costate)
           {
             standardGeneric("setCostateEq")
           }
)

setMethod(f = "setCostateEq",
          signature = "odeEquations",
          definition = function(theObject,costate)
          {
            theObject@costateEq <- costate
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setOrigEq",
           def = function(theObject,origEq)
           {
             standardGeneric("setOrigEq")
           }
)

setMethod(f = "setOrigEq",
          signature = "odeEquations",
          definition = function(theObject,origEq)
          {
            theObject@origEq <- origEq
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setJacobian",
           def = function(theObject,jacobian)
           {
             standardGeneric("setJacobian")
           }
)

setMethod(f = "setJacobian",
          signature = "odeEquations",
          definition = function(theObject,jacobian)
          {
            theObject@jacobian <- jacobian
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name="setHamiltonian",
           def = function(theObject,hamiltonian)
           {
             standardGeneric("setHamiltonian")
           }
)

setMethod(f = "setHamiltonian",
          signature = "odeEquations",
          definition = function(theObject,hamiltonian)
          {
            theObject@hamiltonian <- hamiltonian
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "calculateCostate",
           def = function(theObject) {
             standardGeneric("calculateCostate")
           }
)

setMethod(f = "calculateCostate",
          signature = "odeEquations",
          definition = function(theObject)
          {
            tempList <- symbolicDiff(theObject)
            theObject@costateEq <- tempList$costate
            theObject@invSimMatrix <- tempList$invSimMatrix
            theObject@JhT <- tempList$JhT
            theObject@jacobian <- tempList$jacobian
            theObject@origEq <- tempList$origEq
            theObject@hamiltonian <- tempList$Hamilton
            return(theObject)
          }
)

setGeneric(name = "createModelEqClass",
           def = function(theObject,theModel) {
             standardGeneric("createModelEqClass")
           }
)

setMethod(f = "createModelEqClass",
          signature = "odeEquations",
          definition = function(theObject,theModel)
          {
            theObject@modelStr <- deparse(theModel,width.cutoff = 500)
            theObject@origEq <- getEquations(theModel)
            return(theObject)
          }

)

setGeneric(name = "setCostFunc",
           def = function(theObject,costFunction){
             standardGeneric("setCostFunc")
           }
)

setMethod(f = "setCostFunc",
          signature = "odeEquations",
          definition = function(theObject,costFunction)
          {
            theObject@costFunction <- getEquations(costFunction)
            theObject@dynamicElasticNet <- FALSE
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "setMeassureFunc",
           def = function(theObject,meassureFunc) {
             standardGeneric("setMeassureFunc")
           }
)

setMethod(f = "setMeassureFunc",
          signature = "odeEquations",
          definition = function(theObject,meassureFunc)
          {
            theObject@measureStr <- deparse(meassureFunc,width.cutoff = 500)
            theObject@measureFunction <- getEquations(meassureFunc)
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "isDynElaNet",
           def = function(theObject) {
             standardGeneric("isDynElaNet")
           }
)

setMethod(f = "isDynElaNet",
          signature = "odeEquations",
          definition = function(theObject)
          {
            theObject@dynamicElasticNet <- TRUE
            return(theObject)
          }
)