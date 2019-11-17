devtools::load_all()


# This script generates data to test the greedy-approach dynamic elastic net
library(deSolve)


# The nominal system 
parameters = c( 0.1 , 0.1 , 0.1 , 0.1 , 
                0.2 ,  0.2 ,  0.2 ,  0.9 )
x0 = c(1,0.01,0.01,0.01)

testModel <- function(t,x,parameters) {
  a1       = parameters[1]
  a2       = parameters[2]
  a3       = parameters[3]
  a4       = parameters[4]
  b1       = parameters[5]
  b2       = parameters[6]
  b3       = parameters[7]
  b4       = parameters[8]
  
  dx1 =  b4 * x[4] - a1 * x[1]
  dx2 =  b1 * x[1] - a2 * x[2]
  dx3 =  b2 * x[2] - a3 * x[3]
  dx4 =  b3 * x[3] - a4 * x[4]
  
  return(list(c(dx1,dx2,dx3,dx4)))
}

testMessure <- function(x, parameters=c(0)) {
  
  
  y1 = x[,1] + x[,2] + 2*x[,3] 
  y2 = x[,2] + 2*x[,3]
  
  return(cbind(y1,y2))
}




# Solve the true (perturbated) system with deSolve

# set the time data and the states to use it in deSolve
t.data <- seq(0, 10, length.out = 40)
parameters.solve <- c()
state <- c( x1 = x0[1], x2 = x0[2] , x3 = x0[3], x4 = x0[4])

Perturbation <- function(t,x1,x2,x3,x4) {
  #perturbation function
  perturbation.temp <-  1 / (1 + ( t - 5 )**2)
  
  return(perturbation.temp)
}

Lorenz <- function(t, state, parameters.solve) {
  with(as.list(c(state, parameters.solve)),{
    dx1 <- parameters[8] * x4 - parameters[1] * x1 + Perturbation(t,x1,x2,x3,x4)
    dx2 <- parameters[5] * x1 - parameters[2] * x2 
    dx3 <- parameters[6] * x2 - parameters[3] * x3 
    dx4 <- parameters[7] * x3 - parameters[4] * x4 
    list(c(dx1,dx2,dx3,dx4))
  })
}

out <- ode(y = state, times = t.data, func = Lorenz, parms = parameters.solve)

# MAke plots

plot(out)

t.plot  <- out[,1]
x1.plot <- out[,2]
x2.plot <- out[,3]
x3.plot <- out[,4]
x4.plot <- out[,5]
y.plot <- rep(0, length(t.plot) )

for (i in 1:length(t.plot) ) {
  y.plot[i] <- Perturbation(t.plot[i], x1.plot[i], x2.plot[i], x3.plot[i], x4.plot[i])
}

# how the perturbation should look like
plot(t.plot, y.plot)

# export data
output <- data.frame("t" = out[,1],"x1" = out[,2], "x2" = out[,3], "x3" = out[,4], "x4" = out[,5])

stepAlpha <- 1
times <- t.data

X <- out[,-1]
meas <- as.data.frame(testMessure(X))
meas = as.data.frame(cbind(times,meas))

y <- meas

optW <- c(1,1,1,1)
STD <- NULL


results <- sgdn(alphaStep = 1, alpha2 = 0.05, 
                x0 = x0, optW = optW, 
                measFunc= testMessure,  measData = y,
                parameters = parameters, 
                modelFunc = testModel, plotEstimates = T)

plot(results[[2]])

BDENode <- odeModel(func = testModel, parms = parameters, times=t.data,sd=y*0.05,
                    measFunc = testMessure, y = x0, meas = y,custom=TRUE)


A <- BDEN(odeModel               = BDENode,
          settings               = SETTINGS,
          mcmc_component         = MCMC_component,
          loglikelihood_func     = LOGLIKELIHOOD_func,
          gibbs_update           = GIBBS_update,
          ode_sol                = ode_solv,
          lambda            = .0001,
          beta_init         = c(1,1),
          numbertrialsstep = 15,
          numbertrialseps  = 1000,
          numbertrialinner  = 25)