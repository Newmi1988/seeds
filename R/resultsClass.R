
plot.results <- function(obj) {
  plot.new()
  numMeas <- dim(obj$allData$measurements)[2]
  dimPlot <- ceiling(sqrt(numMeas+1))
  par(mfrow=c(dimPlot,dimPlot))
  barplot(obj$auc, col = 'red')
  inputs <- obj$hiddenInputs[,2:ncol(obj$hiddenInputs)]
  input = inputs[,colSums(inputs)!=0]
  matplot(x = obj$hiddenInputs[,1], y = input, type = 'l', col = 'red', xlab = 'w', ylab = '')
  # plot the measurements against estimates
  for( i in 2:numMeas) {
    yLab <- paste0('y',as.character(i-1))
    yMax <- max(max(obj$estimatedMeasurements[,i]),max(obj$allData$measurements[,i]))
    yMin <- min(min(obj$estimatedMeasurements[,i]),min(obj$allData$measurements[,i]))
    matplot(x = obj$allData$measurements[,1],y = obj$allData$measurements[,i], type = 'l', col = 'black', xlab = 't', ylab = yLab, ylim = c(yMin, yMax))
    par(new=T)
    matplot(x = obj$estimatedMeasurements[,1], y = obj$estimatedMeasurements[,i], type='l', col = 'red', xlab = 't', ylab = yLab, ylim = c(yMin, yMax))
  }

}

results <- setRefClass("results",
                       fields = list(modelFunction = "character", 
                                     measureFunction = "character",
                                     hiddenInputs = "matrix",
                                     auc = "matrix",
                                     estimatedStates = "matrix",
                                     estimatedMeasurements = "data.frame",
                                     optimalSolution = "numeric",
                                     allData = "list"
                                     ),
                         methods = list(
                           show = function() {
                              dimW <- sum(colSums(.self$hiddenInputs)>0)
                              cat(paste0('estimated needed hidden inputs needed ', num2str(dimW), ' at given points.\n'))
                              print(.self$hiddenInputs)
                           }
                         )
                         )