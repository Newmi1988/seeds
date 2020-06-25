library(deSolve)
data(uvbModel)
results <- DEN(odeModel = uvbModel, alphaStep = 500, alpha2 = 0.0001,
                 epsilon = 0.2, plotEstimates = TRUE)