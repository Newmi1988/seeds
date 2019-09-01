modelStr <- 'BIOMD0000000545_url.xml'
uvb <- importSBML(modelStr)

y <- uvbData[,1:3]
uvb <- setMeas(theObject = uvb, y)

# Plot the nominal solution
# plot(nominalSol(odeModel = uvb))

res <- greedyApproach(odeModel = uvb, alphaStep = 400, alpha2 = 0.0001, plotEstimates = TRUE)
plot(res)
