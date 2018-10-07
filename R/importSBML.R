importSBML <- function(modelStr) {
  if(!require('SBMLR',character.only = TRUE)) {
    cat('Please install rsbml from the Bioconducture reposotory')
  } else {
    model <- rsbml_read(filename = modelStr, dom = TRUE)
    
    species <- model@model@species
    parameter <- model@model@parameters
    
    print(paraVec)
    
    # measurements
    rules <- model@model@rules

    # every reaction is a class of type reaction
    reactions <- model@model@reactions
    
    for (i in 1:length(reactions)) {
      print(model@model@reactions[[i]]@kineticLaw)
    }
    
    stoichM <- stoichiometryMatrix(object = model@model)
    print(stoichM)
  }
  
  return(model)
}


