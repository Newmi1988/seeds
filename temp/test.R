

source('./temp/Init_BDEN.R')



DATA <- Data_Model()


##################################################################################




BDEN(observation_time = as.matrix(DATA$observations[["time"]]),
observations     = DATA$observations,
initialvalues    = DATA$X_0,
parameters       = DATA$parameters,
inputData        = as.matrix(DATA$inputData),
numberstates     = DATA$N,
std=DATA$variance,
settings= SETTINGS,
mcmc_component = MCMC_component,
loglikelihood_func =LOGLIKELIHOOD_func,
gibbs_update = GIBBS_update,
ode_sol=ode_solv,

numbertrialsstep = 8,
numbertrialseps  = 500,
numbertrialinner = 15,
lambda           = .001)
