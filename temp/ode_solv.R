ode_solv <- function(TIME,x_0,parameter,inputFUNCTION,w_estimate,STEP,nom_ode,NAME='JAKSTAT'){

  if (NAME == 'JAKSTAT') {
    state <- c(x1 = x_0[['x1_0']], x2 = x_0[['x2_0']], x3 = x_0[['x3_0']], x4 = x_0[['x4_0']])

   #  state <- paste0('[',state[1],',',state[2],',',state[3],',',state[4],']')
   #  w = paste0('[',w_estimate[1,1],',',w_estimate[1,2],',',w_estimate[1,3],',',w_estimate[1,4],';',w_estimate[2,1],',',w_estimate[2,2],',',w_estimate[2,3],',',w_estimate[2,4],']')
   #  x_0 = state
   #  if (STEP>0) {
   #  tn=TIME[STEP]
   #  t0=TIME[STEP-1]
   #  }else{
   #    t0=0
   #   tn=60}
   # 
   # print(paste0('[T,X]=dummyODEJS(',w,',',t0,',',tn,',',x_0,',',STEP,');'))
   #  evaluate(matlab, paste0('[T,X]=dummyODEJS(',w,',',t0,',',tn,',',x_0,',',STEP,');'))
   #  res <- getVariable(matlab, c("T","X"))
   #  res_bind <- cbind(res$T,res$X)
   #  colnames(res_bind) <- c('time','x1','x2','x3','x4')
   #  res_bind <- as.data.frame(res_bind)
   #  return(res_bind)
   #  
   #  code = c(paste0('[T,X]=dummyODEJS(',w,',',t0,',',tn,',',x_0,',',STEP,');'),
   #                 'OUTPUT = [T,X]',
   #                  "save('/home/benjamin/output_eval_Matlab.csv', 'OUTPUT', '-ascii')")
   #  res = run_matlab_code(code)
   #  output = as.data.frame(read.csv("/home/benjamin/output_eval_Matlab.csv",sep='',header=F))
   #  names(output) = c('time','x1','x2','x3','x4')
   #  return(output)

#vode
   #sol        <- R.utils::evalWithTimeout(ode(y = state, method = 'ode45',times = TIME, func = nom_ode, parms = parameter, inputD = inputData, w=w_estimate,TIME_W =TIME[1], name = NAME, rtol = 1e-6, atol = 1e-10),
     #                                  timeout = 3.5, onTimeout = "silent")   

   sol        <- ode(y = state,times = TIME, func = nom_ode, parms = parameter, input = inputFUNCTION, w=w_estimate,TIME_W =TIME[1], name = NAME)
   
   if (!is.null(sol)) sol[sol> -0.00001&sol<0] <-0
    
   
    
   if (is.null(sol)|((sum(sol< 0)!=0))){
    return(NA)}
    
    sol        <- sol[,1:5]
 
    return(as.data.frame(sol[,1:5]))
  
}}