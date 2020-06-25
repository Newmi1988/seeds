library(devtools)
data(bden_uvb)
results <- BDEN(odeModel          = Model,
                 lambda            = .001,
                 beta_init         = c(1,1,1,1,1),
                 numbertrialsstep  = 15,
                 numbertrialseps   = 2000,
                 numbertrialinner  = 10)