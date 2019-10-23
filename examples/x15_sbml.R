devtools::load_all()

t <- c(0,60)
model <- importSBML("BIOMD0000000015.xml", times = t, y = y)

# model -> model -> reactions -> kineticLaw -> parameters
model@kineticLaw@parameters

       