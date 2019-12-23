# For testing use http://biomodels.caltech.edu/BIOMD0000000545#Files
t <- uvbData[,1]
y <- uvbData[,1:3]
uvb <- importSBML("BIOMD0000000545_url.xml", times = t, meas = y)

# Plot the nominal solution
nominalSol(odeModel = uvb)