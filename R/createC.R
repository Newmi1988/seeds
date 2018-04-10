createCFile <- function(parameters, inputs,Eq, bden){
  
  if(missing(bden)){
    bden <- FALSE
  }
  
  StringC <- '#include <R.h>'
  StringC = append(StringC,'#include <math.h>')

  
  # define parameters and inputs
  if(bden){
    paraStr <- paste0('static double parms[',length(parameters)+inputs,'];')
  } else {
    paraStr <- paste0('static double parms[',length(parameters),'];')
  }


  inpStr <- paste0('static double forc[', inputs+1,'];')

  
  
  StringC = append(StringC, values = c('',paraStr,inpStr))
  
  #define parameters
  if(is.null(names(parameters))){
    para <- paste0('#define ',Eq@parameters,' parms[',1:length(Eq@parameters)-1,']')
  } else {
    para <- paste0('#define ',names(parameters),' parms[',1:length(parameters)-1,']')
  }
  
  if(bden){
    bdent0 = paste0("#define t0 parms[",length(parameters),"]")
    bdenwt0 = paste0("#define w",1:(inputs-1),'t0 parms[',(length(parameters)+1):(length(parameters)+inputs-1),']')
    bdenPara = append(bdent0,bdenwt0)
    para = append(para,bdenPara)
  }
  
  #define inputs
    defInputU <- paste0('#define u forc[',0,']')
    defInputs <- paste0('#define w',1:(inputs-1),' forc[',1:(inputs-1),']')
    StringC = append(StringC, values = c('',para,'',defInputU,defInputs))

  
  tStr <- rep("",5)
  # format the functions
  tStr[1] = "void parmsc(void (* odeparms)(int *, double *))"
  tStr[2] = "{"
  if(bden){
    tStr[3] = paste0("\tint N=",length(parameters)+inputs,';')
  } else {
    tStr[3] = paste0("\tint N=",length(parameters),';')
  }
  tStr[4] = "\todeparms(&N, parms);"
  tStr[5] = "}"
  
  StringC = append(StringC, values = c('',tStr))
  
  tStr[1] = "void forcc(void (* odeforcs)(int *, double *))"
  tStr[3] = paste0("\tint N=",inputs,';')
  tStr[4] = "    odeforcs(&N, forc);"

  StringC = append(StringC, values = c('',tStr))
  
  # the ode function
  startStr <- "void derivsc(int *neq, double *t, double *x, double *dx, double *yout, int *ip)\n{"
  if(bden){
    conBden = "\tif(*t == t0){"
    t0Bden = paste0('\t\tw',1:(inputs-1),' = w',1:(inputs-1),'t0;')
    t0Bden = append(t0Bden,"\t}")
    conBden = append(conBden,t0Bden)
    startStr = append(startStr,conBden)
  }
  
  
  eqC <- gsub(pattern = "(d*[x])([0-9]*)", replacement = "\\1[\\2]" , Eq@origEq)
  eqC = gsub(pattern = "(\\[[0-9]*)", replacement = "\\1 -1", eqC)
  eqC = unlist(lapply(X = eqC, FUN = Deriv::Simplify))
  eqC = gsub(pattern = "(x*\\[[0-9]*\\])\\^[0-9]*", replacement = "\\1*\\1", eqC)
  eqC = gsub(pattern = "(t)([^a-z])", replacement = "*\\1\\2", eqC)

  eqC = paste0("\t",eqC,"+w",1:length(eqC),";")

  StringC = append(StringC, values = c(startStr,eqC,'}'))
  
  
  writeFileC <- function(string){
    file.create('model.c')
    fileC <- file('model.c')
    writeLines(string,fileC)
    close(fileC)
  }
  
  writeFileC(StringC)
} 

#' Create compilable c-code of a model
#' 
#' Writes a c file that can be compiled for faster solution with the \code{\link[deSolve]{ode}} solver.
#' The file created is formated to be used with the dynamic elastic net. A hidden input is 
#' added to every component of the state vector.
#' 
#' @note On the usage of compiled code in conjunction with \pkg{deSolve} take a look into the vignette 'R
#' Package deSolve, Writing Code in Compiled Languages' of the package.
#' 
#' @param modelFunc a R-function that can be solved with deSolve. External input of the system should
#'                  be declared with 'u'. To ensure that the function is working use the most general
#'                  state-space representation.
#'
#' @param parameters a vector describing the parameters of the system. If names are missing the function
#'                   tries to extract the declared parameters from the model function. 
#'                   
#' @param bden a boolean that indicates if the c-file is used for the mcmc algorithm
#'                   
#' @return None
#'                   
#' @export
createCompModel <- function(modelFunc, parameters, bden){
  odeEq <- new("odeEquations")
  odeEq <- createModelEqClass(odeEq,modelFunc)
  
  numInputs = length(odeEq@origEq)+1
  createCFile(parameters = parameters,inputs = numInputs, odeEq, bden)
}
