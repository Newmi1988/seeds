if(is.function(input)) {
  inputFun <- input
} else {
  inputFun <- stats::approxfun()
}


Output

list
$state estimates[t x1 x2 x...]
$state uncentaintylower[t x1 x2 x...]
$state uncentaintyupper[t x1 x2 x...]
$hidden input estimates[t w1 w2 w3 ...]
$hidden input uncentaintylower[t w1 w2 w...]
$hidden input uncentaintyupper[t w1 w2 w...]
$output estimates[t y1 y2 ...]
$Data[t y1 y2 ...]
$DataError

plot <- ggplot(OUTPUT$state%>%geather(),aes(x=t,y=))+
  geom_line()+
  geom_ribbon(data=OUTPUT$stateuc%>%geather()%>%full_join(
              OUTPUT$stateuc%>%geather()),
              aes=(y_min=lower,y_max=upper))+
  scale_color_discrete()
                
        
plot <- ggplot(OUTPUT$estimates%>%geather(),aes(x=t,y=))+
  geom_line()+
  geom_errorbar(DataError%>%geather()!JOIN!)
  geom_point(Data%>%geather(),aes(x=t,y=))+
  geom_ribbon(data=OUTPUT$estucuc%>%geather()%>%full_join(
    OUTPUT$estuc%>%geather()),
    aes=(y_min=lower,y_max=upper))+
  scale_color_discrete()+
  facet_wrap(~.label)


[$likelihood trace
 $gibbs parameter trace]
  
  
  
  ggplot(A%>%gather(Label,estimate,2:3),aes(x=time,y=estimate,colour=Label))+geom_line()+theme_bw()+  facet_wrap(~Label)

  
  Data
  $ode
  $measurments
  $stimuli
  $standarderror
  $inputfunction
  
  ParameterALGoRITHM <=>
    MCMC_SET$STEP_trials      <- 8#15
  MCMC_SET$EPS_step_size          <- 100*DATA$N
  MCMC_SET$EPS_step_size_inner    <- 15
  MCMC_SET$BURNIN           <- round(MCMC_SET$EPS_step_size/3)
  MCMC_SET$BURNIN_inner     <- round(MCMC_SET$EPS_step_size_inner/3)
  MCMC_TRACK                <- vector("list",dim(X_MODEL)[1])        
  
  SETTINGS MODEL Specific
  
DATA,BETA_LAMBDA,MCMC_SET$STEP_trials,MCMC_SET$EPS_step_size,MCMC_SET$EPS_step_size_inner,nom_ode,objective} 
  
S default??!
  SIGMA default!?!
