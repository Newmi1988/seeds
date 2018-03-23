
source('./temp/Init_BDEN.R')
COUNTER = 1
COUNTER2 = 1

##################################################################################
DATA <- Data_Model()
inputF <- approxfun(x = DATA$inputData[,1], y = DATA$inputData[,2], method = 'linear', rule=2)
##################################################################################
X_MODEL        <- ode_solv(DATA$observations[["time"]],DATA$X_0,DATA$parameters,inputF,matrix(rep(0,2*DATA$N),2),0,nom_ode,DATA$name)
X_ERROR        <- cbind(abs(abs(DATA$observations["STAT5"])-abs(X_MODEL["x1"])),abs(abs(DATA$observations["STAT5ptot_cyt"])-DATA$parameters["s1"]*abs(X_MODEL["x2"]+2*X_MODEL["x3"])),abs(abs(DATA$observations["STAT5p_cyt"])-DATA$parameters["s2"]*abs(X_MODEL["x1"]+X_MODEL["x2"]+2*X_MODEL["x3"])))
X_ERROR_SMOOTH <- cbind(predict(smooth.spline(X_MODEL[["time"]],X_ERROR[["STAT5"]]),X_MODEL[["time"]])[["y"]],predict(smooth.spline(X_MODEL[["time"]],X_ERROR[["STAT5p_cyt"]]),X_MODEL[["time"]])[["y"]],predict(smooth.spline(X_MODEL[["time"]],X_ERROR[["STAT5ptot_cyt"]]),X_MODEL[["time"]])[["y"]])
GRADIENT       <- cbind(abs(diff(X_ERROR[,1])/diff(X_MODEL[["time"]])),abs(diff(X_ERROR[,2])/diff(X_MODEL[["time"]])),abs(diff(X_ERROR[,3])/diff(X_MODEL[["time"]])))
##################################################################################
MCMC_SET                  <- list()
GIBBS_PAR_IT              <- list()
EPSILON_IT                <- list()
VAR                       <- list()
MCMC_SET$STEP_trials      <- 8#15
MCMC_SET$EPS_step_size          <- 100*DATA$N
MCMC_SET$EPS_step_size_inner    <- 15
MCMC_SET$BURNIN           <- round(MCMC_SET$EPS_step_size/3)
MCMC_SET$BURNIN_inner     <- round(MCMC_SET$EPS_step_size_inner/3)
MCMC_TRACK                <- vector("list",dim(X_MODEL)[1])        

EPSILON          <- matrix(0, nrow = dim(X_MODEL)[1]  , ncol = DATA$N)
SIGMA            <- vector("list",dim(X_MODEL)[1]) 
S                <- mean(GRADIENT)*2#tuning para
for (i in 1:length(SIGMA)){
  SIGMA[[i]]      <- max(abs(diff(GRADIENT)))
}
BETA_LAMBDA      <- .001
GIBBS_PAR        <- SETTINGS(DATA$variance,DATA$N,BETA_LAMBDA,DATA$name)
EPS_TIME         <- DATA$observations$time
YINIT            <- DATA$X_0
EPSILON_IT$CONT  <- matrix(0, nrow = MCMC_SET$STEP_trials  , ncol = DATA$N)
EPSILON_IT$ACT   <- matrix(0, nrow = 2  , ncol = DATA$N)
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
   
     MCMC_RESULTS                 <-MCMC_component(LOGLIKELIHOOD_func, MCMC_SET$EPS_step_size, MCMC_SET$EPS_step_size_inner, EPSILON_IT$CONT[TRIALS-1,],S,
                                    STEP,DATA$observations,EPSILON_IT$Y0,DATA$parameters,EPSILON_IT$ACT,inputF,VAR$SIGMA,VAR$DIAG,GIBBS_PAR,DATA$N,MCMC_SET$BURNIN_inner,nom_ode,DATA$name)
  
     MCMC_RESULT_THIN            <- coda::mcmc(MCMC_RESULTS, start = MCMC_SET$BURNIN,end=dim(MCMC_RESULTS)[1],thin=10)
    
    EPSILON_IT$CONT[TRIALS,]    <- colMeans(MCMC_RESULT_THIN)
    EPSILON_IT$ACT[2,]          <- EPSILON_IT$CONT[TRIALS,] 
    
    G_U                         <- GIBBS_update(VAR$DIAG,EPSILON_IT$ACT,
                                                GIBBS_PAR$R,GIBBS_PAR$ROH,SIGMA[[1]],DATA$N,VAR$SIGMA,
                                                GIBBS_PAR_IT$LAMBDA2,GIBBS_PAR_IT$LAMBDA1,GIBBS_PAR_IT$TAU)  
    VAR$SIGMA                   <- G_U$SIGMA 
    GIBBS_PAR_IT$LAMBDA2        <- G_U$LAMBDA2 
    GIBBS_PAR_IT$LAMBDA1        <- G_U$LAMBDA1 
    GIBBS_PAR_IT$TAU            <- G_U$TAU 
  }
  EPSILON[STEP,]              <- colMeans(EPSILON_IT$CONT[-1,])
  SOLUTION                    <- ode_solv(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,DATA$parameters,inputF,EPSILON[c(STEP-1,STEP),],EPS_TIME[c(STEP-1,STEP)],nom_ode,DATA$name)
  YINIT                       <- setNames(SOLUTION[2,seq(2,dim(SOLUTION)[2],1)],c("x1_0","x2_0","x3_0","x4_0"))
  SIGMA[[STEP]]               <- VAR$SIGMA
}



