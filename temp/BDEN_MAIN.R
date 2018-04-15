
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


BDEN <- function(observation_time,
                 observations,
                 initialvalues,
                 parameters,
                 inputData,
                 numberstates,
                 std,
                 settings,
                 model,
                 mcmc_component,
                 loglikelihood_func,
                 gibbs_update,
                 ode_sol,
                 model,
                 measFunc, 
                 
                 numbertrialsstep = 8,
                 numbertrialseps  = 500,
                 numbertrialinner = 15,
                 lambda           = .001){


  #### Function welche die C files erstellt. ####
  
  
  createCompModel(modelFunc = model, parameters = parameter, bden = TRUE)
  ext <- .Platform$dynlib.ext #erkenne plattform um die shared libary unter Linux und Windoes laden zu können
  
  # name des shared library objects
  compiledModel <- paste0('model',ext)
  
  # falls die shared library schon geladen ist muss diese wieder "entladen" werden,
  # damit der compiler diese neu schreiben kann (falls Änderungen gemacht wurden am Modell usw.)
  if(is.loaded('derivsc')){
    dyn.unload(compiledModel)
  }
  # kompilieren des modells
  system("R CMD SHLIB model.c")
  # laden der shared lib.
  dyn.load(compiledModel)
  

  ##################################################################################
  X_MODEL        <- ode_sol(observation_time,initialvalues,parameters,inputData,matrix(rep(0,2*4),2))
  
  X_ERROR        <- cbind(abs(abs(observations["STAT5"])-abs(X_MODEL["x1"])),abs(abs(observations["STAT5ptot_cyt"])-parameters["s1"]*abs(X_MODEL["x2"]+2*X_MODEL["x3"])),abs(abs(observations["STAT5p_cyt"])-parameters["s2"]*abs(X_MODEL["x1"]+X_MODEL["x2"]+2*X_MODEL["x3"])))
  
  X_ERROR_SMOOTH <- cbind(predict(smooth.spline(observation_time,X_ERROR[["STAT5"]]),observation_time)[["y"]],predict(smooth.spline(observation_time,X_ERROR[["STAT5p_cyt"]]),observation_time)[["y"]],predict(smooth.spline(observation_time,X_ERROR[["STAT5ptot_cyt"]]),observation_time)[["y"]])
  
  GRADIENT       <- cbind(abs(diff(X_ERROR[,1])/diff(observation_time)),abs(diff(X_ERROR[,2])/diff(observation_time)),abs(diff(X_ERROR[,3])/diff(observation_time)))
  
##################################################################################
  COUNTER = 1
  COUNTER2 = 1
  
MCMC_SET                  <- list()
GIBBS_PAR_IT              <- list()
EPSILON_IT                <- list()
VAR                       <- list()
MCMC_SET$STEP_trials      <- numbertrialsstep
MCMC_SET$EPS_step_size          <- numbertrialseps
MCMC_SET$EPS_step_size_inner    <- numbertrialinner
MCMC_SET$BURNIN           <- round(MCMC_SET$EPS_step_size/3)
MCMC_SET$BURNIN_inner     <- round(MCMC_SET$EPS_step_size_inner/3)
MCMC_TRACK                <- vector("list",dim(X_MODEL)[1])        

EPSILON          <- matrix(0, nrow = dim(X_MODEL)[1]  , ncol = numberstates)
SIGMA            <- vector("list",dim(X_MODEL)[1]) 
S                <- mean(GRADIENT)*2#tuning para
for (i in 1:length(SIGMA)){
  SIGMA[[i]]      <- max(abs(diff(GRADIENT)))
}
BETA_LAMBDA      <- lambda
GIBBS_PAR        <- SETTINGS(std,DATA$N,BETA_LAMBDA)
EPS_TIME         <- observation_time
YINIT            <- initialvalues
EPSILON_IT$CONT  <- matrix(0, nrow = MCMC_SET$STEP_trials  , ncol = numberstates)
EPSILON_IT$ACT   <- matrix(0, nrow = 2  , ncol = numberstates)
##################################################################################


for (STEP in 2:length(EPS_TIME)){
 
  EPSILON_IT$Y0            <- YINIT
  GIBBS_PAR_IT$TAU         <- GIBBS_PAR$TAU
  GIBBS_PAR_IT$LAMBDA1     <- GIBBS_PAR$LAMBDA1
  GIBBS_PAR_IT$LAMBDA2     <- GIBBS_PAR$LAMBDA2
  EPSILON_IT$ACT[1,]       <- EPSILON[STEP-1,]
  EPSILON_IT$CONT[1,]      <- EPSILON[STEP-1,]
  VAR$SIGMA                <- SIGMA[[STEP]]  
  
  for (TRIALS in 2:MCMC_SET$STEP_trials){
    print(TRIALS)
    VAR$DIAG                    <- diag((GIBBS_PAR_IT$TAU+GIBBS_PAR_IT$LAMBDA2)^-1)
   
     MCMC_RESULTS                 <- mcmc_component(loglikelihood_func, MCMC_SET$EPS_step_size, MCMC_SET$EPS_step_size_inner, EPSILON_IT$CONT[TRIALS-1,],S,
                                     STEP,observations,EPSILON_IT$Y0,inputData,parameters,EPSILON_IT$ACT,VAR$SIGMA,VAR$DIAG,GIBBS_PAR,numberstates,MCMC_SET$BURNIN_inner,measFunc)
  
     MCMC_RESULT_THIN            <- coda::mcmc(MCMC_RESULTS, start = MCMC_SET$BURNIN,end=dim(MCMC_RESULTS)[1],thin=10)
    
    EPSILON_IT$CONT[TRIALS,]    <- colMeans(MCMC_RESULT_THIN)
    EPSILON_IT$ACT[2,]          <- EPSILON_IT$CONT[TRIALS,] 
    
    G_U                         <- gibbs_update(VAR$DIAG,EPSILON_IT$ACT,
                                                GIBBS_PAR$R,GIBBS_PAR$ROH,SIGMA[[1]],DATA$N,VAR$SIGMA,
                                                GIBBS_PAR_IT$LAMBDA2,GIBBS_PAR_IT$LAMBDA1,GIBBS_PAR_IT$TAU)  
    VAR$SIGMA                   <- G_U$SIGMA 
    GIBBS_PAR_IT$LAMBDA2        <- G_U$LAMBDA2 
    GIBBS_PAR_IT$LAMBDA1        <- G_U$LAMBDA1 
    GIBBS_PAR_IT$TAU            <- G_U$TAU 
  }
  EPSILON[STEP,]              <- colMeans(EPSILON_IT$CONT[-1,])
  SOLUTION                    <- ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILON[c(STEP-1,STEP),])
  
  SIGMA[[STEP]]               <- VAR$SIGMA
}

}

