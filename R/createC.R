createCFile <- function(parameters, inputs, Eq, extInp){
  
  StringC <- '#include <R.h>'

  
  # define parameters and inputs
  paraStr <- paste0('static double parms[',length(parameters),'];')
  if(!is.null(extInp)) {
    inpStr <- paste0('static double forc[', inputs+1,'];')
  } else {
    inpStr <- paste0('static double forc[', inputs,'];')
  }

  
  
  StringC = append(StringC, values = c('',paraStr,inpStr))
  
  #define parameters
  if(is.null(names(parameters))){
    para <- paste0('#define p',1:length(parameters)-1,' parms[',1:length(parameters)-1,']')
  } else {
    para <- paste0('#define ',names(parameters),' parms[',1:length(parameters)-1,']')
  }
  
  #define inputs
  if(!is.null(extInp)){
    defInputU <- paste0('#define u forc[',0,']')
    defInputs <- paste0('#define import',1:(inputs),' forc[',1:(inputs),']')
    StringC = append(StringC, values = c('',para,'',defInputU,defInputs))
  } else {
    defInputs <- paste0('#define import',1:(inputs),' forc[',1:(inputs)-1,']')
    StringC = append(StringC, values = c('',para,'',defInputs))
  }
  
  tStr <- rep("",5)
  # format the functions
  tStr[1] = "void parmsc(void (* odeparms)(int *, double *))"
  tStr[2] = "{"
  tStr[3] = paste0("\tint N=",length(parameters),';')
  tStr[4] = "\todeparms(&N, parms);"
  tStr[5] = "}"
  
  StringC = append(StringC, values = c('',tStr))
  
  tStr[1] = "void forcc(void (* odeforcs)(int *, double *))"
  tStr[3] = paste0("\tint N=",inputs,';')
  tStr[4] = "    odeforcs(&N, forc);"

  StringC = append(StringC, values = c('',tStr))
  
  # the ode function
  startStr <- "void derivsc(int *neq, double *t, double *x, double *dx, double *yout, int *ip)\n{"
  
  eqC <- gsub(pattern = "(d*[x])([0-9]*)", replacement = "\\1[\\2]" ,Eq)
  eqC = gsub(pattern = "(x*\\[[0-9]*\\])\\^[0-9]*", replacement = "\\1*\\1", eqC)
  eqC = gsub(pattern = "(\\[[0-9]*)", replacement = "\\1 -1", eqC)

  # eqC = gsub(pattern = "u*\\s*", replacement = "import0", eqC)
  eqC = paste0("\t",eqC,"+import",1:length(eqC),";")
  
  
  StringC = append(StringC, values = c(startStr,eqC,'}'))
  
  
  writeFileC <- function(string){
    file.create('model.c')
    fileC <- file('model.c')
    writeLines(string,fileC)
    close(fileC)
  }
  
  writeFileC(StringC)
} 
