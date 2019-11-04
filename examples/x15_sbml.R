t <- seq(0,2,0.1)
model <- importSBML("BIOMD0000000015.xml", times = t)
plot(nominalSol(model))
