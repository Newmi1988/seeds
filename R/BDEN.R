#' Bayesian Dynamic Elastic Net
#'
#' Full Bayesian algortihm to detect hidden inputs in ODE based models.The algortihm 
#' is an extention of the Dynamic Elastic Net algorithm (Engelhardt et al. 2016) inpiered by the Elastic-Net Regression.
#'
#' @param observation_time     observed time points
#' @param observations         observed state dynamics e.g. protein concentrations
#' @param initialvalues        initial values
#' @param parameters           model parameters 
#' @param inputData            discrete input function e.g. stimuli
#' @param numberstates         number of modeled states
#' @param std                  standard error
#' @param settings             initial model specific settings
#' @param model                ODE system
#' @param mcmc_component       used sampling algorithm
#' @param loglikelihood_func   used likelihood function
#' @param gibbs_update         used gibbs algorithm
#' @param ode_sol              used ode solver
#' @param measFunc             link function to match observations with modeled states
#' @param numbertrialsstep     number of sampels per timepoint
#' @param numbertrialseps      number of samples per mcmc step
#' @param numbertrialinner     number of inner samples
#' @param lambda               inital shrinkage parameter
#' @param Grad_correct         used for intial mcmc step size calculation 
#' @param alpha                mcmc tuning paramter
#' @param beta                 mcmc tunig parameter
#'
#' @return                     returns a results-object with default plot function
#'
#' @example /examples/exampleBDEN.R
#' 
#' @export
#' 
#' 

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
                 measFunc, 
                 
                 numbertrialsstep = 15,
                 numbertrialseps  = 500*4,
                 numbertrialinner = 30,
                 lambda           = .001,
                 Grad_correct     = 0,
                 alpha            = c(1,1,1,1  ),
                 beta             = c(1,1,1,0.1)){
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  createCompModel(modelFunc = model, parameters = parameters, bden = TRUE)
  ext <- .Platform$dynlib.ext 
  
  compiledModel <- paste0('model',ext)
  
  
  if(is.loaded('derivsc')){
    dyn.unload(compiledModel)
  }
  
  system("R CMD SHLIB model.c")
  
  dyn.load(compiledModel)
  
  
  PROGRESS <- R.utils::ProgressBar(max=numbertrialsstep*numbertrialseps*numbertrialinner, ticks= numbertrialseps , stepLength=1, newlineWhenDone=FALSE)
  
  ##################################################################################
  X_MODEL        <- ode_sol(observation_time,initialvalues,parameters,inputData,matrix(rep(0,2*4),2))
  
  
  X_ERROR        <- matrix(0,length(observation_time),numberstates)
  
  for (i in 1:length(observation_time)){
    X_ERROR [i,] <- abs(abs(as.numeric(observations[i,-1]))-abs(as.numeric(sapply(1:4,measFunc,y=X_MODEL[i,],parameter=parameters[5:6],USE.NAMES = TRUE))))
  }
  
  #  X_ERROR_SMOOTH   <-        predict(smooth.spline(observation_time,X_ERROR[,1]),observation_time)[["y"]]
  #  for (i in 2:(numberstates-1)){
  #   X_ERROR_SMOOTH <- cbind(X_ERROR_SMOOTH,predict(smooth.spline(observation_time,X_ERROR[,i]),observation_time)[["y"]])
  # }
  
  GRADIENT        <- matrix(0,length(observation_time)-1,numberstates-1)
  
  for (i in 1:(numberstates-(1+Grad_correct))){   
    GRADIENT[,i]       <- abs(diff(X_ERROR[,i])/diff(observation_time))
  }
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
  MCMC_SET$BURNIN           <- round(MCMC_SET$EPS_step_size*(3/5))
  MCMC_SET$BURNIN_inner     <- round(MCMC_SET$EPS_step_size_inner/3)
  MCMC_TRACK                <- vector("list",dim(X_MODEL)[1])        
  
  EPSILON          <- matrix(0, nrow = dim(X_MODEL)[1]  , ncol = numberstates)
  EPSILONLOW       <- matrix(0, nrow = dim(X_MODEL)[1]  , ncol = numberstates)
  EPSILONUP        <- matrix(0, nrow = dim(X_MODEL)[1]  , ncol = numberstates)
  SIGMA            <- vector("list",dim(X_MODEL)[1]) 
  S                <- mean(GRADIENT)*2
  for (i in 1:length(SIGMA)){
    SIGMA[[i]]      <- max(abs(diff(GRADIENT)))*0.1
  }
  BETA_LAMBDA      <- lambda
  GIBBS_PAR        <- SETTINGS(std,numberstates,BETA_LAMBDA,alpha,beta)
  EPS_TIME         <- observation_time
  YINIT            <- initialvalues
  EPSILON_IT$CONT  <- matrix(0, nrow = MCMC_SET$STEP_trials  , ncol = numberstates)
  EPSILON_IT$ACT   <- matrix(0, nrow = 2  , ncol = numberstates)
  SOLUTION         <- YINIT
  SOLUTIONLOW      <- YINIT*0
  SOLUTIONUP       <- YINIT*0
  EPSILON_IT$NEW   <- YINIT
  print('S')
  print(S)
  print('SIGMA')
  print(SIGMA)
  print(  GIBBS_PAR)
  ##################################################################################
  
  
  for (STEP in 2:length(EPS_TIME)){
    print(STEP)
    
    EPSILON_IT$Y0            <- EPSILON_IT$NEW
    GIBBS_PAR_IT$TAU         <- GIBBS_PAR$TAU
    GIBBS_PAR_IT$LAMBDA1     <- GIBBS_PAR$LAMBDA1
    GIBBS_PAR_IT$LAMBDA2     <- GIBBS_PAR$LAMBDA2
    EPSILON_IT$ACT[1,]       <- EPSILON[STEP-1,]
    EPSILON_IT$CONT[1,]      <- EPSILON[STEP-1,]
    VAR$SIGMA                <- SIGMA[[STEP]]  
    
    for (TRIALS in 2:MCMC_SET$STEP_trials){
      R.utils::increase(PROGRESS)
      print(TRIALS)
      print('DIAGONAL')
      VAR$DIAG                    <- diag((GIBBS_PAR_IT$TAU+GIBBS_PAR_IT$LAMBDA2)^-1)
      print(VAR$DIAG)
      
      
      MCMC_RESULTS                 <- mcmc_component(loglikelihood_func, MCMC_SET$EPS_step_size, MCMC_SET$EPS_step_size_inner, EPSILON_IT$CONT[TRIALS-1,],S,
                                                     STEP,observations,EPSILON_IT$Y0,inputData,parameters,EPSILON_IT$ACT,VAR$SIGMA,VAR$DIAG,GIBBS_PAR,numberstates,MCMC_SET$BURNIN_inner,measFunc)
      
      MCMC_RESULT_THIN            <- coda::mcmc(MCMC_RESULTS, start = MCMC_SET$BURNIN,end=dim(MCMC_RESULTS)[1],thin=20)
      
      
      
      EPSILON_IT$CONT[TRIALS,]    <- colMeans(MCMC_RESULT_THIN[-1,])
      EPSILON_IT$ACT[2,]          <- EPSILON_IT$CONT[TRIALS,] 
      
      G_U                         <- gibbs_update(VAR$DIAG,EPSILON_IT$ACT,
                                                  GIBBS_PAR$R,GIBBS_PAR$ROH,SIGMA[[1]],numberstates,VAR$SIGMA,
                                                  GIBBS_PAR_IT$LAMBDA2,GIBBS_PAR_IT$LAMBDA1,GIBBS_PAR_IT$TAU)  
      VAR$SIGMA                   <- G_U$SIGMA 
      GIBBS_PAR_IT$LAMBDA2        <- G_U$LAMBDA2 
      GIBBS_PAR_IT$LAMBDA1        <- G_U$LAMBDA1 
      GIBBS_PAR_IT$TAU            <- G_U$TAU 
      print(paste0('SIGMA: ', VAR$SIGMA ))

    }
    
    EPSILON[STEP,]              <- matrixStats::colQuantiles(EPSILON_IT$CONT[-1,], probs = 0.5, na.rm = TRUE)
    EPSILONLOW[STEP,]           <- matrixStats::colQuantiles(EPSILON_IT$CONT[-1,], probs = 0.05, na.rm = TRUE) 
    EPSILONUP[STEP,]            <- matrixStats::colQuantiles(EPSILON_IT$CONT[-1,], probs = 0.95, na.rm = TRUE)
    
    
    
    UP                          <- as.numeric(tail(ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILONUP[c(STEP-1,STEP),]),1))
    
    LOW                         <- as.numeric(tail(ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILONLOW[c(STEP-1,STEP),]),1))
    
    print(EPSILON[c(STEP-1,STEP),])
    EPSILON_IT$NEW              <- as.numeric(tail(ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILON[c(STEP-1,STEP),]),1))
    
    
    
    SOLUTION                    <- rbind(SOLUTION,EPSILON_IT$Y0)
    SOLUTIONLOW                 <- rbind(SOLUTIONLOW,LOW)
    SOLUTIONUP                  <- rbind(SOLUTIONUP ,UP)
    
    SIGMA[[STEP]]               <- VAR$SIGMA
  }
  
  
  X_OUTPUT <- matrix(0,length(observation_time),numberstates)
  for (i in 1:length(observation_time)){
    X_OUTPUT [i,] <- as.numeric(sapply(1:4,measFunc,y=SOLUTION[i,],parameter=parameters[5:6],USE.NAMES = TRUE))
  }
  
  print(EPSILONLOW)
  print(EPSILONUP)
  hiddenInpUnsclower <- EPSILONLOW
  colnames(hiddenInpUnsclower)[1] <- "t"
  
  hiddenInpUnscupper <- EPSILONUP
  colnames(hiddenInpUnscupper)[1] <- "t"
  
  stateUnscertainlower <- SOLUTIONLOW 
  colnames(stateUnscertainlower)[1] <- "t"
  
  stateUnscertainupper <- SOLUTIONUP 
  colnames(stateUnscertainupper)[1] <- "t"
  
  states <- as.data.frame(cbind(observation_time,SOLUTION))
  colnames(states)[1] <- "t"
  
  hiddenInp <- as.data.frame(cbind(observation_time,EPSILON))
  colnames(hiddenInp)[1] <- "t"
  
  outputMeas <- as.data.frame(cbind(observation_time,X_OUTPUT))
  colnames(outputMeas)[1] <- "t"
  
  nomStates <- as.data.frame(cbind(observation_time,X_MODEL))
  colnames(nomStates)[1] = "t"
  
  dataError <- std
  colnames(dataError) <- c("t",paste0('y',1:(ncol(sd)-1)))
  
  measData <- observations
  colnames(measData) <- c("t",paste0('y',1:(ncol(measData[,-1]))))
  
  res <- seeds::resultsSeeds(stateNominal = nomStates,
                             stateEstimates = states,
                             stateUnscertainLower =  stateUnscertainlower,
                             stateUnscertainUpper =  stateUnscertainupper,
                             hiddenInputEstimates = hiddenInp,
                             hiddenInputUncertainLower = hiddenInpUnsclower,
                             hiddenInputUncertainUpper = hiddenInpUnscupper,
                             outputEstimates = outputMeas,
                             Data = measData,
                             DataError = dataError
  )
  
  print(res)
  return(res)
  
}





