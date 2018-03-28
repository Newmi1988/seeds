resStates <- res$estimatedStates

resHidden <- matrix(res$hiddenInputs)
resMeas <- matrix(res$estimatedMeasurements)
resData <- matrix(res$allData$measurements)
resDataStd <- matrix(cbind(uvbData[,1],uvbData[,7:11]))

resultsSeeds(stateEstimates = resStates, hiddenInputEstimates = resHidden, 
             outputEstimates = resMeas, Data = resData, DataError = resDataStd)

