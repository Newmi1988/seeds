# For testing use http://biomodels.caltech.edu/BIOMD0000000545#Files

xmlfile <- 'BIOMD0000000545_url.xml'
uvb <- importSBML(xmlfile)

y <- uvbData[,1:3]
uvb <- setMeas(theObject = uvb, y)

# Plot the nominal solution
plot(nominalSol(odeModel = uvb))

res <- sgdn(odeModel = uvb, alphaStep = 400, alpha2 = 0.0001, plotEstimates = TRUE)
plot(res)
