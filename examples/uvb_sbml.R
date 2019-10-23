# For testing use http://biomodels.caltech.edu/BIOMD0000000545#Files
rm(list = ls())
devtools::load_all()

t <- uvbData[,1]
y <- uvbData[,1:3]
uvb <- importSBML("BIOMD0000000545_url.xml", times = t, meas = y)


# Plot the nominal solution
# plot(nominalSol(odeModel = uvb))

res <- sgdn(odeModel = uvb, alphaStep = 400, alpha2 = 0.0001, plotEstimates = TRUE, measData = y)
