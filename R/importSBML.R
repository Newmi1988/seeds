importSBML <- function(modelStr) {
  if(!require('rsbml',character.only = TRUE)) {
    cat('Please install rsbml from the Bioconducture reposotory')
  } else {
    model <- rsbml_read(filename = modelStr, dom = TRUE)
    
    species <- model@model@species
    parameter <- model@model@parameters
    
    # print(paraVec)
    
    # measurements
    rules <- model@model@rules
    reactions <- model@model@reactions
    
    reacList <- list()
    for (i in 1:length(reactions)) {
      reacList[[i]] <- gsub(pattern = "expression",replacement = '',deparse(model@model@reactions[[i]]@kineticLaw@math, width.cutoff = 300))
    }

    stoichM <- stoichiometryMatrix(object = model@model)
    print(stoichM)

    react <- c()
    combieReact <- function(reactStrs, stMatrix) {
      for (i in 1:nrow(stMatrix)) {
        m <- which(stMatrix[i,] !=0 )
        if(length(m)>0) {
          react <- c(react,paste0(stMatrix[i,m],'*',reactStrs[m], collapse = ' + '))
        }
      }
      print(react)
    }

    combieReact(reacList,stoichM)
    
  }
  
  meas <- c()
  for (i in 1:length(rules)) {
    meas[[i]] <- gsub(pattern = "expression", replacement = '', x = rules[[i]]@math)
  }
  
  print(meas)
  
  return(model)
}


