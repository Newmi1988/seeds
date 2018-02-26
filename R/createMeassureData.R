createMeassureData <- function(output,noise) {
  times = output[,1]
  output = output[,-1]
  
  m = nrow(output)
  n = ncol(output)
  
  outputSave = as.matrix(output)
  
  outputMeasurement <- list()
  A <- matrix(rep(0,m*n), nrow = m, ncol = n)

  for(i in 1:m) {
    for(j in 1:n) {
      A[,j] = output[,j]#+rnorm(length(output[,j]),mean = mean(output[,j]),sd = mean(output[,j])*noise)
    }
    A[which(A < 0)] = outputSave[which(outputSave < 0)]
    outputMeasurement[[i]] = A
  }

  it <- length(outputMeasurement)
  c = 0
  vectorConvert <- matrix(rep(0,m*n*it), nrow = it)
  for (i in 1:m) {
    for (j in 1:n) {
      c = c+1
      for (k in 1:it) {
        vectorConvert[k,c] = outputMeasurement[[k]][i,j]
      }
    }
  }
  
  xMean <- colMeans(vectorConvert)
  output = as.data.frame(cbind(times,matrix(xMean, ncol = n, byrow = TRUE)))
  names(output)[-1] = paste0(rep("y",ncol(output)-1),1:(ncol(output)-1))
  

  xSTD <- apply(vectorConvert, 2, sd)
  STD = as.data.frame(matrix(xSTD, nrow = m, byrow = TRUE))
  names(STD) <- paste0(rep("sd",ncol(STD)),rep("y",ncol(STD)) ,1:ncol(STD))
  
  
  res <- cbind(output,STD)
  
  return(res)
}
