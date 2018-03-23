nom_ode  <- function(t,x,parameter,inputD,w,TIME_W,name = 'JAKSTAT'){
  if (name == 'JAKSTAT'){ 
    dx     <- rep(0,4)
    if (t == TIME_W[1]){POINT <- 1}
    if (t > TIME_W[1]) {POINT <- 2}
    u <- inputF(t)
    dx[1] =  - parameter[[1]] * u * x[1]                                 + w[POINT,1];
    dx[2] =  + parameter[[1]] * u * x[1] - parameter[[2]] * x[2]*x[2]    + w[POINT,2];
    dx[3] =  + parameter[[2]] * x[2] * x[2] * 0.5  - parameter[[3]] * x[3]                  + w[POINT,3];
    dx[4] =  + parameter[[3]] * x[3]  + w[POINT,4]; 

    
    list(dx)
    

  }
}


#objective
objective  <- function(index,y,parameter,NAME = 'JAKSTAT'){
  

  if (NAME == 'JAKSTAT'){ 
    if (index == 1){
      return(y[1])
    }
    if (index == 2){
      return(parameter[1] * (y[2] + 2 * y[3]))
    }
    if (index == 3){
      return(parameter[2] * (y[1] + y[2] + 2 * y[3]))
    }
    if (index == 4){
      return(y[4])
    }
  }
}