#' Bayesian Dynamic Elastic Net
#'
#' Full Bayesian algortihm to detect hidden inputs in ODE based models.The algortihm 
#' is an extension of the Dynamic Elastic Net algorithm (Engelhardt et al. 2016) inspired by the Elastic-Net Regression.
#' 
#' Ordinary differential equations (ODEs) are a popular approach to quantitatively model molecular networks based on biological knowledge. 
#' However, such knowledge is typically restricted. Wrongly modelled biological mechanisms as well as relevant external influence factors 
#' that are not included into the model likely manifest in major discrepancies between model predictions and experimental data. 
#' Finding the exact reasons for such observed discrepancies can be quite challenging in practice. 
#' In order to address this issue we suggest a Bayesian approach to estimate hidden influences in ODE based models. 
#' The method can distinguish between exogenous and endogenous hidden influences. Thus, we can detect wrongly specified as well as missed 
#' molecular interactions in the model. 
#' The BDEN as a new and fully probabilistic approach, supports the modeller in an algorithmic manner to identify possible sources of errors 
#' in ODE based models on the basis of experimental data.  THE BDEN does not require pre-specified hyper-parameters. 
#' BDEN thus provides a systematic Bayesian computational method to identify target nodes and reconstruct the corresponding 
#' error signal including detection of missing and wrong molecular interactions within the assumed model. 
#' The method works for ODE based systems even with uncertain knowledge and noisy data. 
#' In contrast to approaches based on point estimates the Bayesian framework incorporates the given uncertainty and circumvents 
#' numerical pitfalls which frequently arise from optimization methods (Engelhardt et al. 2017).
#' 
#'
#' @param observation_time     observed time points
#' @param observations         observed state dynamics e.g. protein concentrations
#' @param initialvalues        initial values of the system
#' @param parameters           model parameters estimates
#' @param inputData            discrete input function e.g. stimuli
#' @param numberstates         number of system states
#' @param sd                  standard error of the observed stat dynamics (per time point)
#' @param settings             initial model specific settings (autmaticly calculated based on the nominal model and data)
#' @param model                ODE system
#' @param mcmc_component       used sampling algorithm
#' @param loglikelihood_func   used likelihood function
#' @param gibbs_update         used gibbs algorithm
#' @param ode_sol              used ode solver
#' @param measFunc             link function to match observations with modeled states. Takes states, index of observed variable and involved parameters. Returns the estimated observation at the given index.
#' @param LogTransform         use the log transformed ODE system 
#' @param numbertrialsstep     number of gibbs updates per timepoint. This should be at least 10. Values have direct influnce on the runtime. 
#' @param numbertrialseps      number of samples per mcmc step. This should be greater than numberStates*500.Values have direct influnce on the runtime.
#' @param numbertrialinner     number of inner samples. This should be greater 15 to guarantee a reasonable exploration of the sample space. Values have direct influnce on the runtime.
#' @param lambda               inital shrinkage parameter.
#' @param Grad_correct         used for intial mcmc step size calculation 
#' @param alpha                mcmc tuning paramter (weigthing of observed states)
#' @param beta                 mcmc tunig parameter (weigthing of observed states)
#'
#' @return                     returns a results-object with default plot function
#'
#' @example /examples/exampleBDEN.R
#' 
#' @export
#' 
#' 



BDEN <- function(measData,
                 x0,
                 parameters,
                 systemInput,
                 sd,
                 settings,
                 mcmc_component,
                 loglikelihood_func,
                 gibbs_update,
                 ode_sol,
                 measFunc, 
                 modelFunc,
                 LogTransform     = FALSE,
                 numbertrialsstep = 15,
                 numbertrialseps  = 500*4,
                 numbertrialinner = 10,
                 lambda           = .001,
                 Grad_correct     = 0,
                 alpha            = c(1,1,1,1),
                 beta_init        = c(1,1,1,0.1)){
  
  
  
  
  observation_time <- measData[,1]
  observations     <- measData

  
  
   initialvalues    <- x0
   inputData        <- as.matrix(systemInput)
   model            <- modelFunc
   numberstates     <- length(x0)
  
  
  if(LogTransform) {createCompModel(modelFunc = model, parameters = parameters, bden = TRUE,logTransVar=1:numberstates)}
  else{createCompModel(modelFunc = model, parameters = parameters, bden = TRUE)}
  ext <- .Platform$dynlib.ext 
  
  compiledModel <- paste0('model',ext)
  
  
  if(is.loaded('derivsc')){
    dyn.unload(compiledModel)
  }
  
  system("R CMD SHLIB model.c")
  
  dyn.load(compiledModel)
  
  
  PROGRESS <- R.utils::ProgressBar(max=numbertrialsstep*numbertrialseps*numbertrialinner, ticks= numbertrialseps , stepLength=1, newlineWhenDone=FALSE)
  
  ##################################################################################
  X_MODEL        <- ode_sol(observation_time,initialvalues,parameters,inputData,matrix(rep(0,2*4),2),LogTransform)

  print('Nominal State Dynamics')
  print(X_MODEL)
  print('#################################')
  
  X_ERROR        <- matrix(0,length(observation_time),numberstates)
  
  for (i in 1:length(observation_time)){
    X_ERROR [i,] <- abs(abs(as.numeric(observations[i,-1]))-abs(as.numeric(sapply(1:4,measFunc,y=X_MODEL[i,],parameter=parameters[5:6],USE.NAMES = TRUE))))
  }

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
  S                <- mean(GRADIENT)*10
  for (i in 1:length(SIGMA)){
    SIGMA[[i]]      <- max(abs(diff(GRADIENT)))*0.5
  }
  BETA_LAMBDA      <- lambda
  print(numberstates)
  print(BETA_LAMBDA)
  print(beta_init)
  GIBBS_PAR        <- SETTINGS(sd,numberstates,BETA_LAMBDA,alpha,beta_init)
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
  print(GIBBS_PAR)
  print
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
      print('TAU')
      print(GIBBS_PAR_IT$TAU)
      print('LAMBDA')
      print(GIBBS_PAR_IT$LAMBDA2)
      print('DIFF EPSILON')
      print(EPSILON_IT$ACT[2,]-EPSILON_IT$ACT[1,])
      print(paste0('SIGMA'))
      print(VAR$SIGMA)     
      print('#####################################################################')
      
      MCMC_RESULTS                 <- mcmc_component(loglikelihood_func, MCMC_SET$EPS_step_size, MCMC_SET$EPS_step_size_inner, EPSILON_IT$CONT[TRIALS-1,],S,
                                                     STEP,observations,EPSILON_IT$Y0,inputData,parameters,EPSILON_IT$ACT,VAR$SIGMA,VAR$DIAG,GIBBS_PAR,numberstates,MCMC_SET$BURNIN_inner,measFunc,LogTransform)
      
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
      
      

    }
    
    EPSILON[STEP,]              <- matrixStats::colQuantiles(EPSILON_IT$CONT[-1,], probs = 0.5, na.rm = TRUE)
    EPSILONLOW[STEP,]           <- matrixStats::colQuantiles(EPSILON_IT$CONT[-1,], probs = 0.05, na.rm = TRUE) 
    EPSILONUP[STEP,]            <- matrixStats::colQuantiles(EPSILON_IT$CONT[-1,], probs = 0.95, na.rm = TRUE)
    
    
    
    UP                          <- as.numeric(tail(ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILONUP[c(STEP-1,STEP),],LogTransform),1))
    
    LOW                         <- as.numeric(tail(ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILONLOW[c(STEP-1,STEP),],LogTransform),1))
    


     EPSILON_IT$NEW              <- as.numeric(tail(ode_sol(EPS_TIME[c(STEP-1,STEP)],EPSILON_IT$Y0,parameters,inputData,EPSILON[c(STEP-1,STEP),],LogTransform),1))
    
    
    SOLUTION                    <- rbind(SOLUTION,EPSILON_IT$Y0)
    SOLUTIONLOW                 <- rbind(SOLUTIONLOW,LOW)
    SOLUTIONUP                  <- rbind(SOLUTIONUP ,UP)
    
    SIGMA[[STEP]]               <- VAR$SIGMA
  }
  


  

  hiddenInpUnsclower <- as.data.frame(cbind(observation_time,EPSILONLOW))
  colnames(hiddenInpUnsclower) = c('t','w1','w2','w3','w4')
  
  hiddenInpUnscupper <- as.data.frame(cbind(observation_time,EPSILONUP))
  colnames(hiddenInpUnscupper) = c('t','w1','w2','w3','w4')
  
  stateUnscertainlower <- as.data.frame(cbind(observation_time,SOLUTIONLOW))
  colnames(stateUnscertainlower)= c('t','x1','x2','x3','x4')
  
  
  stateUnscertainupper <- as.data.frame(cbind(observation_time,SOLUTIONUP))
  colnames(stateUnscertainupper)= c('t','x1','x2','x3','x4')
  
  states <- as.data.frame(cbind(observation_time,SOLUTION))
  colnames(states)= c('t','x1','x2','x3','x4')
  
  hiddenInp <- as.data.frame(cbind(observation_time,EPSILON))
  colnames(hiddenInp)= c('t','w1','w2','w3','w4')
  
  
  X_OUTPUT <- matrix(0,length(observation_time),numberstates)
  for (i in 1:length(observation_time)){
    X_OUTPUT [i,] <- as.numeric(sapply(1:4,measFunc,y=SOLUTION[i,],parameter=parameters[5:6],USE.NAMES = TRUE))
  }
  
  outputMeas <- as.data.frame(cbind(observation_time,X_OUTPUT))
  colnames(outputMeas) = c('t','y1','y2','y3','y4')
  
  nomStates <- as.data.frame(cbind(observation_time,X_MODEL))
  colnames(nomStates) = c('t','x1','x2','x3','x4')

  dataError <- sd
  colnames(dataError) <- c("t",paste0('s',1:(ncol(sd)-1)))

  measData <- observations
  colnames(measData) <- c("t",paste0('y',1:(ncol(measData[,-1]))))
  
  res <- seeds::resultsSeeds(stateNominal = nomStates,
                             stateEstimates = states,
                             stateUnscertainLower =  stateUnscertainlower,
                             stateUnscertainUpper =  stateUnscertainupper,
                             hiddenInputEstimates = hiddenInp,
                             hiddenInputUncertainLower = hiddenInpUnsclower,
                             hiddenInputUncertainUpper = hiddenInpUnscupper,
                             outputEstimatesUncLower = NA,
                             outputEstimatesUncUpper = NA,
                             outputEstimates = outputMeas,
                             Data = measData,
                             DataError = dataError
  )
  

  return(res)
  
}




