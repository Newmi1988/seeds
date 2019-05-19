#' A S4 class used to handle formating ODE-Equation and calculate the needed functions for the seeds-algorithm
#'
#' @slot modelStr a vector of strings describing the ODE
#' @slot measureStr a vector of strings representing the equation of the measurement function
#' @slot origEq a vector of strings containing the original model function
#' @slot measureFunction a vector of strings containing the original measurement function
#' @slot costateEq a vector of strings describing the costate equation
#' @slot JhT a matrix of strings describing the jacobian matrix of the measurement function
#' @slot jacobian a matrix of strings representing the jacobian matrix model equations
#' @slot costFunction a string containing the cost function
#' @slot hamiltonian a string representing the Hamilton function of the model
#' @slot dynamicElasticNet boolean that indicates if the system equation should be calculated for the dynamic elastic net
#' @slot parameters parameters of the model
#' @slot cond a slot to save conditionals in equations, which are used for formating the c files
#' @slot logInd a slot containing a vector which indicates state variables that should be log-transformed
#' 
#' @export odeEquations
#' @exportClass odeEquations
odeEquations <- setClass(
  #name of Class
  "odeEquations",
  slots = c(
    modelStr = "character",
    measureStr = "character",
    origEq = "character",
    measureFunction = "character",
    costateEq = "character",
    JhT = "matrix",
    jacobian = "matrix",
    costFunction = "character",
    hamiltonian = "character",
    dynamicElasticNet = "logical",
    parameters = "character",
    cond = "list",
    logInd = "numeric"
  ),

  prototype = list(
    modelStr = character(0),
    measureStr = character(0),
    origEq = character(0),
    measureFunction =  character(0),
    costateEq = character(0),
    costFunction = character(0),
    JhT = matrix(list(),nrow = 2, ncol = 2),
    jacobian = matrix(list(),nrow = 2, ncol = 2),
    costFunction = character(0),
    hamiltonian = character(0),
    dynamicElasticNet = FALSE,
    parameters = character(0),
    cond = list(),
    logInd = numeric(0)
  ),
  validity = function(object) {
    # no validation needed
  }
)

setMethod('initialize', "odeEquations", function(.Object,...){
            .Object <- callNextMethod()
            return(.Object)
          })


# setCostateEq
# assign the costate the calculated costate function
# to the odeEquationObject 
# the costate equation is used to calculate the optimal control
setGeneric(name="setCostateEq",
           def = function(odeEquationObject,costate)
           {
             standardGeneric("setCostateEq")
           }
)

setMethod(f = "setCostateEq",
          signature = "odeEquations",
          definition = function(odeEquationObject,costate)
          {
            odeEquationObject@costateEq <- costate

            return(odeEquationObject)
          }
)



# setOrigEq
# assign the original equation to the odeEquationObject
setGeneric(name="setOrigEq",
           def = function(odeEquationObject,origEq)
           {
             standardGeneric("setOrigEq")
           }
)

setMethod(f = "setOrigEq",
          signature = "odeEquations",
          definition = function(odeEquationObject,origEq)
          {
            odeEquationObject@origEq <- origEq

            return(odeEquationObject)
          }
)



# setJacobian
# internal function
# assign a jacobian function to the odeEquation object
setGeneric(name="setJacobian",
           def = function(odeEquationObject,jacobian)
           {
             standardGeneric("setJacobian")
           }
)

setMethod(f = "setJacobian",
          signature = "odeEquations",
          definition = function(odeEquationObject,jacobian)
          {
            odeEquationObject@jacobian <- jacobian

            return(odeEquationObject)
          }
)



# set the hamiltonian
# should not be exposed directly be the user
setGeneric(name="setHamiltonian",
           def = function(odeEquationObject,hamiltonian)
           {
             standardGeneric("setHamiltonian")
           }
)

setMethod(f = "setHamiltonian",
          signature = "odeEquations",
          definition = function(odeEquationObject,hamiltonian)
          {
            odeEquationObject@hamiltonian <- hamiltonian

            return(odeEquationObject)
          }
)



# Calculate Costate
# calculates the costate function from the given ode-model function
setGeneric(name = "calculateCostate",
           def = function(odeEquationObject) {
             standardGeneric("calculateCostate")
           }
)

setMethod(f = "calculateCostate",
          signature = "odeEquations",
          definition = function(odeEquationObject)
          {
            tempList <- symbolicDiff(odeEquationObject)
            odeEquationObject@costateEq <- tempList$costate
            odeEquationObject@JhT <- tempList$JhT
            odeEquationObject@jacobian <- tempList$jacobian
            odeEquationObject@origEq <- tempList$origEq
            odeEquationObject@hamiltonian <- tempList$Hamilton
            return(odeEquationObject)
          }
)



# assigns a model to the odeEquationClass
# theModel - has to be a R-functions that can be solved with deSolve
setGeneric(name = "createModelEqClass",
           def = function(odeEquationObject,theModel) {
             standardGeneric("createModelEqClass")
           }
)

setMethod(f = "createModelEqClass",
          signature = "odeEquations",
          definition = function(odeEquationObject,theModel)
          {
            odeEquationObject@modelStr <- deparse(theModel,width.cutoff = 500)
            tempList <- getEquations(theModel)
            odeEquationObject@origEq <- tempList$strM
            odeEquationObject@parameters <- tempList$strP
            odeEquationObject@cond <- tempList$cond
            return(odeEquationObject)
          }

)



# set costum costfunction
setGeneric(name = "setCostFunc",
           def = function(odeEquationObject,costFunction){
             standardGeneric("setCostFunc")
           }
)

setMethod(f = "setCostFunc",
          signature = "odeEquations",
          definition = function(odeEquationObject,costFunction)
          {
            odeEquationObject@costFunction <- getEquations(costFunction)$strM
            odeEquationObject@dynamicElasticNet <- FALSE

            return(odeEquationObject)
          }
)



# setMeassureFunc
# class method for setting the measurement function
# the functions has to be a standard r-function
setGeneric(name = "setMeassureFunc",
           def = function(odeEquationObject,meassureFunc) {
             standardGeneric("setMeassureFunc")
           }
)

setMethod(f = "setMeassureFunc",
          signature = "odeEquations",
          definition = function(odeEquationObject,meassureFunc)
          {
            odeEquationObject@measureStr <- deparse(meassureFunc,width.cutoff = 500)
            odeEquationObject@measureFunction <- getEquations(meassureFunc)$strM
            return(odeEquationObject)
          }
)



# isDynElaNet
# Indicate if the elastic net method should be use
# if not set the optimizer will use no hidden inputs
# the solutions are only the optimal control
setGeneric(name = "isDynElaNet",
           def = function(odeEquationObject) {
             standardGeneric("isDynElaNet")
           }
)

setMethod(f = "isDynElaNet",
          signature = "odeEquations",
          definition = function(odeEquationObject)
          {
            odeEquationObject@dynamicElasticNet <- TRUE
            return(odeEquationObject)
          }
)




# Function setLogTransfInd
# Set indeces for components of that shoudl be log transformed
setGeneric(name = "setLogTransInd",
          def = function(odeEquationObject,logIndVec) {
            standardGeneric("setLogTransInd")
          }
)

setMethod(f = "setLogTransInd",
          signature = "odeEquations",
          definition = function(odeEquationObject,logIndVec)
          {
            odeEquationObject@logInd <- logIndVec
            return(odeEquationObject)
          }
)


