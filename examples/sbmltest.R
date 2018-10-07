rm(list=ls())
devtools::load_all()
modelStr <- 'BIOMD0000000545_url.xml'
model <- importSBML(modelStr)
