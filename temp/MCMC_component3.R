MCMC_component <- function(LOGLIKELIHOOD_func, STEP_SIZE, STEP_SIZE_INNER , EPSILON, JUMP_SCALE,
                           STEP,OBSERVATIONS,Y0,INPUTDATA,PARAMETER,EPSILON_ACT,SIGMA,DIAG,GIBBS_par, N, BURNIN){

  number_species = N
  epsilon_container=matrix(rep(0,number_species*(STEP_SIZE+1)*(STEP_SIZE_INNER+1)),((STEP_SIZE+1)*(STEP_SIZE_INNER+1)))

    
    rbind(rep(0,number_species),rep(0,number_species));
    
    
    epsilon_container[2,]=EPSILON_ACT[2,]
    
  
    JUMP_SCALE=JUMP_SCALE/2

for (ii in 2:STEP_SIZE+1){

  eps0=epsilon_container[ii-1,]

  eps1= matrix(rep(eps0,STEP_SIZE_INNER),ncol=number_species,nrow=STEP_SIZE_INNER,byrow=TRUE)

  k = sample(1:number_species,1) 
  FLAG = 0;

  MU_jump = EPSILON_ACT[1,k]
  
  eps1[,k] = rnorm(STEP_SIZE_INNER,MU_jump,JUMP_SCALE)


{

        ratio_new <- sapply(eps1[,k],LOGLIKELIHOOD_func,Step=STEP,OBSERVATIONS=OBSERVATIONS,x_0=Y0,parameters=PARAMETER,EPS_inner=EPSILON_ACT[1,],
                        D=DIAG,GIBBS_PAR=GIBBS_par,k=k,MU_JUMP=MU_jump,SIGMA_JUMP=JUMP_SCALE,eps_new=eps1[1,],INPUT=INPUTDATA)
    
        ratio_old <- LOGLIKELIHOOD_func(epsilon_container[ii-1,k],STEP,OBSERVATIONS,Y0,PARAMETER,EPSILON_ACT[1,],INPUTDATA,DIAG,GIBBS_par,k,MU_jump,JUMP_SCALE,eps_new=epsilon_container[ii-1,])

        ratio     <- exp(ratio_new-ratio_old-dnorm(eps1[,k],MU_jump,JUMP_SCALE,log=TRUE)+dnorm(epsilon_container[ii-1,k],MU_jump,JUMP_SCALE,log=TRUE))

        
     if (sum(is.na(ratio))==STEP_SIZE_INNER){print(paste0('Internal Error; BREAK;',ii))}
    else{
      if (0.95 < max((ratio),na.rm=T)){

        eps_i = eps1[which(ratio==max((ratio),na.rm=T)),];
    FLAG = 1;

      }
    }
  }


  if (FLAG==1){epsilon_container[ii,] = eps_i

  }
else{epsilon_container[ii,] = epsilon_container[ii-1,]} 


}
  return(epsilon_container)
  
} 
