rm(list = ls())
devtools::load_all()

model_obj <- rsbml::rsbml_read(filename = "BIOMD0000000015.xml", dom = TRUE)

t <- c(0,60)
model <- importSBML("BIOMD0000000015.xml", times = t)
