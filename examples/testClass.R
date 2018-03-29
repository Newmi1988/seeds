devtools::load_all()
graphics.off()

resStates <- as.data.frame(res$estimatedStates)
colnames(resStates)[1] <- "t"
resHidden <- as.data.frame(res$hiddenInputs)
colnames(resHidden)[1] <- "t"
resMeas <- as.data.frame(res$estimatedMeasurements)
colnames(resMeas)[1] <- "t"
resData <- as.data.frame(res$allData$measurements)
resDataStd <- cbind(t=uvbData[,1],uvbData[,7:11])
colnames(resDataStd) <- c("t",paste0('y',1:(ncol(resDataStd)-1)))

resObj <- resultsSeeds(stateEstimates = resStates, hiddenInputEstimates = resHidden, 
             outputEstimates = resMeas, Data = resData, DataError = resDataStd)

library(ggplot2)
library(dplyr)
library(tidyr)


plot(resObj)

