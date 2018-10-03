importSBML <- function(model) {
  if(!require('SBMLR',character.only = TRUE)) {
    cat('Please install SBMLR from the Bioconducture Reposotory')
  } else {
    
    model <- SBMLR::readSBML('BIOMD0000000545_url.xml')
    sm <- summary(model)
    reactions <- sm$reactions[,2]
    
    states = sm$S0
    # states[sm$BC == FALSE] = X
    # if(sm$nRules > 0){
    #   for (j in sm$nRules) {
    #     print(model$rules[[j]]$law)
    #     states[model$rules[[2]]$idOutput] =  model$rules[[j]]$law(states[model$rule[[j]]$inputs])
    #   }
    # }
    
    v = rep(0, sm$nReactions)
    xp = rep(0, sm$nStates)
    
    react <- sm$reactions$Laws
    # print(react)
    
    incMatrix <- sm$incid
    print(incMatrix)
    
    # function to combine the reactions and the incidents
    combReact <- function(incMatrix,reactions){

      for (i in 1:nrow(incMatrix)){
        ind <- which(incMatrix[i,]!=0)
        if(length(ind)!=0) {
          x <- incMatrix[i,ind]
          
          print(reactions[ind])
          
          print(x)
        }

      }
    }
    
    combReact(incMatrix, react)
    
    
    
  }
}


