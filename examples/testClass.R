devtools::load_all()
graphics.off()

### Lade Beispieldaten um die Plot-Funktion zu überprüfen
load('uvbRes.RData')
resStates <- as.data.frame(res$estimatedStates)
colnames(resStates)[1] <- "t"
resHidden <- as.data.frame(res$hiddenInputs)
colnames(resHidden)[1] <- "t"
resMeas <- as.data.frame(res$estimatedMeasurements)
colnames(resMeas)[1] <- "t"
resData <- as.data.frame(res$allData$measurements)
resDataStd <- cbind(t=uvbData[,1],uvbData[,7:11])
colnames(resDataStd) <- c("t",paste0('y',1:(ncol(resDataStd)-1)))

# just some data for testing
hInputLow <- cbind(t=resHidden[,1],resHidden[,2:ncol(resHidden)]-10)
hInputUp <- cbind(t=resHidden[,1],resHidden[,2:ncol(resHidden)]+10)

stateLow <- cbind(t=resStates[,1],resStates[,2:ncol(resStates)]-1)
stateUp <- cbind(t=resStates[,1],resStates[,2:ncol(resStates)]+1)

# falls keine Konfidenzintervalle berechnet wurden einfach dummy mit gleicher Struktur und NaN übergeben
# für testen ein und auskommentieren
# stateLow[,2:ncol(stateLow)] = NaN
# stateUp[,2:ncol(stateUp)] = NaN


### Benutzung der 'resultsSeeds'-Klasse
resObj <- resultsSeeds(stateEstimates = resStates, 
                       stateUnscertainLower = stateLow,
                       stateUnscertainUpper = stateUp,
                      hiddenInputEstimates = resHidden, 
                      hiddenInputUncertainLower = hInputLow,
                      hiddenInputUncertainUpper = hInputUp,
                      outputEstimates = resMeas,
                      Data = resData,
                      DataError = resDataStd)




plot(resObj)

