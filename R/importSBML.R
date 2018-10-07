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

    # every reaction is a class of type reaction
    reactions <- model@model@reactions
    
    reacList <- list()
    for (i in 1:length(reactions)) {
      reacList[[i]] <- gsub(pattern = "expression",replacement = '',deparse(model@model@reactions[[i]]@kineticLaw@math, width.cutoff = 300))
    }

    stoichM <- stoichiometryMatrix(object = model@model)
    print(stoichM)

    combieReact <- function(reactStrs, stMatrix) {
      for (i in 1:nrow(stMatrix)) {
        m <- which(stMatrix[i,] !=0 )
        if(length(m)>0) {
          test <- paste0(stMatrix[i,m],'*',reactStrs[m], collapse = ' + ')
          print(test)
        }

      }
    }
    
    combieReact(reacList,stoichM)
    
  }
  
  return(model)
}


