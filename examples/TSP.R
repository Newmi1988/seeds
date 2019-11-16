# This script generates data to test the greedy-approach dynamic elastic net

parameters = c(  v1=1,
                   vi1=1,
                   ki1=1,
                   ni1=2,
                   ka1=1,
                   na1=2,
                   k1=1,
                   v2=1,
                   ki2=1,
                   ni2=2,
                   ka2=1,
                   na2=2,
                   k2=1,
                   v3=1,
                   ki3=1,
                   ni3=2,
                   ka3=1,
                   na3=2,
                   k3=1,
                   v4=0.1,
                   k4=1,
                   k4k=0.1,
                   v5=0.1,
                   k5=1,
                   k5k=0.1,
                   v6=0.1,
                   k6=1,
                   k6k=0.1,
                   kcat1=1,
                   km1=1,
                   km2=1,
                   kcat2=1,
                   km3=1,
                   km4=1,
                   kcat3=1,
                   km5=1,
                   km6=1,
                   p=0.05,
                   s=0.1)

x0 = c(0.66667,0.57254,0.41758,0.4,0.36409,0.29457,1.419,0.93464)


Model <- function(t,x,parameters) {
  with (as.list(parameters),{
    
    dx1 = v1/(1+(p/ki1)*ni1+(ka1/s)*na1)-k1*x[1]
    dx2 = v2/(1+(p/ki2)*ni2+(ka2/x7)*na2)-k2*x[2]
    dx3 = v3/(1+(p/ki3)*ni3+(ka3/x8)*na3)-k3*x[3] 
    dx4 = (v4*x[1])/(k4+x[1])-k4k*x[4]
    dx5 = (v5*x[2])/(k5+x[2])-k5k*x[5]
    dx6 = (v6*x[2])/(k6+x[3])-k6k*x[6]
    dx7 = (kcat1*x[4]*(1/km1)*(s-x[7]))/(1+(s/km1)+(x[7]/km2))-(kcat2*x[5]*(1/km3)*(x[7]-x[8]))/(1+(x[7]/km3)+(x[8]/km4))
    dx8 = (kcat2*x[5]*(1/km3)*(x[7]-x[8]))/(1+(x[7]/km3)+(x[8]/km4))-(kcat3*x[6]*(1/km5)*(x[8]-p))/(1+(x[8]/km5)+(p/km6))
    
    list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8))
  })
}


uvbMeasure <- function(x, parameters) {
  
  y1 = 2*x[,5] + x[,4] + x[,8]
  y2 = 2*x[,5] + 2* x[,3] + x[,1]
  y3 = x[,6]
  y4 = x[,1]
  y5 = x[,4]
  
  return(cbind(y1,y2,y3,y4,y5))
}

t.data <- seq(0, 60, length.out = 60)
state <- c( x1 = x0[1], x2 = x0[2] , x3 = x0[3], x4 = x0[4], x5 = x0[5], x6 = x0[6], x7 = x0[7], x8 = x0[8])

Perturbation <- function(t,x1,x2,x3,x4,x5,x6,x7,x8) {
  #perturbation function
  perturbation.temp <-  1 / (40 + ( t - 30 )**2)
  
  return(perturbation.temp)
}

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dx1 = v1/(1+(p/ki1)*ni1+(ka1/s)*na1)-k1*x1 + Perturbation(t,x1,x2,x3,x4,x5,x6,x7,x8)
    dx2 = v2/(1+(p/ki2)*ni2+(ka2/x7)*na2)-k2*x2
    dx3 = v3/(1+(p/ki3)*ni3+(ka3/x8)*na3)-k3*x3 
    dx4 = (v4*x1)/(k4+x1)-k4k*x4
    dx5 = (v5*x2)/(k5+x2)-k5k*x5
    dx6 = (v6*x2)/(k6+x3)-k6k*x6
    dx7 = (kcat1*x4*(1/km1)*(s-x7))/(1+(s/km1)+(x7/km2))-(kcat2*x5*(1/km3)*(x7-x8))/(1+(x7/km3)+(x8/km4))
    dx8 = (kcat2*x5*(1/km3)*(x7-x8))/(1+(x7/km3)+(x8/km4))-(kcat3*x6*(1/km5)*(x8-p))/(1+(x8/km5)+(p/km6))
    list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8))
  })
}

out <- deSolve::ode(y = state, times = t.data, func = Lorenz, parms = parameters)

# MAke plots

plot(out)
X <- out[,-1]
meas <- uvbMeasure(X)
meas = as.data.frame(cbind(t.data,meas))
y <- meas

output <- data.frame("t" = out[,1],"x1" = out[,2], "x2" = out[,3], "x3" = out[,4], "x4" = out[,5], "x5" = out[,6], "x6" = out[,7], "x7" = out[,8], "x8" = out[,9])


sd=NULL

model_class <- odeModel(func = Model, parms = parameters, times=t.data,
                    measFunc = uvbMeasure, y = x0, meas = y,custom=TRUE)

plot(nominalSol(model_class))

res <- sgdn(odeModel = model_class, alphaStep = 0.0001, alpha2 = 0.1, plotEstimates = TRUE, conjGrad = FALSE)

plot(res[[1]])

BDENode <- odeModel(func = Model, parms = parameters, times=t.data,sd=y[-1]*0.05,
                    measFunc = uvbMeasure, y = x0, meas = y,custom=TRUE)


A <- BDEN(odeModel               = BDENode,
          settings               = SETTINGS,
          mcmc_component         = MCMC_component,
          loglikelihood_func     = LOGLIKELIHOOD_func,
          gibbs_update           = GIBBS_update,
          ode_sol                = ode_solv,
          lambda            = 1,
          beta_init         = c(1,1,1,1,1),
          numbertrialsstep = 15,
          numbertrialseps  = 1000,
          numbertrialinner  = 25)
