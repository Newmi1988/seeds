#### clear workspace, load package, remove graphics ####

rm(list=ls())
devtools::load_all()
graphics.off()

#### setze initialen zustand, parameter und Zeitpunkte der Auswertung durch deSolve ####
# parameter des models
N = 10^0.31
parameters = 10^c("k1"=0.31, "k2"=-1, "k3"=-0.49, "k4"= 0.42, "s1"=-0.21, "s2"=-0.34)

#N = 2.4290
#parameters = c("k1"=2.4290, "k2"=975.4280, "k3"=0.1157, "k4"= 0, "s1"=10^-0.21, "s2"=10^-0.34)

x0 = c(N, 0.00001, 0.00001, 0.00001)

inputData                    <- read.table('http://jeti.uni-freiburg.de/PNAS_Swameye_Data/DATA1_hall_inp')#csv
#inputData                    <- read.table('../data/DATA1_hall_inp.txt')#csv
inputData[nrow(inputData),2] <- 0.009
colnames(inputData)          <- c('t','u')
measure                      <- read.table('http://jeti.uni-freiburg.de/PNAS_Swameye_Data/DATA1_hall')#csv
#measure                      <- read.table('../data/DATA1_hall.txt')#csv
colnames(measure)            <- c("time","STAT5p_cyt" ,"sd_STAT5p_cyt","STAT5ptot_cyt","sd_STAT5ptot_cyt")

sd                           <- cbind(measure['time'],measure['sd_STAT5p_cyt'],measure['sd_STAT5ptot_cyt'])
y                            <- cbind(measure['time'],((measure['STAT5ptot_cyt']/parameters['s2'])-(measure['STAT5p_cyt']/parameters['s1'])),measure['STAT5p_cyt'],measure['STAT5ptot_cyt'],(x0[1]-(measure['STAT5ptot_cyt']/parameters['s2']))/2/(1400/450))

# ---
# Manipulate support points 

inputApprox   <- apply(X = inputData[,-1,drop=FALSE], MARGIN = 2, FUN = function(x) stats::approx(x = inputData[,1], y = x, xout = seq(0,60,5), rule = 2))
inputData     <- (cbind(inputApprox$u$x,inputApprox$u$y))
inputData     <- as.data.frame(inputData)
colnames(inputData)          <- c('t','u')

inputApprox <- apply(X = measure[,c(2,4)], MARGIN = 2, FUN = function(x) stats::approx(x = measure[,1], y = x, xout = seq(0,60,5), rule = 2))
measur     <- (cbind(inputApprox$STAT5p_cyt$x,inputApprox$STAT5p_cyt$y,inputApprox$STAT5ptot_cyt$y))
colnames(measur)            <-  c("time","STAT5p_cyt" ,"STAT5ptot_cyt")
y                            <- cbind(measur[,'time'],((measur[,'STAT5ptot_cyt']/parameters['s2'])-(measur[,'STAT5p_cyt']/parameters['s1'])),measur[,'STAT5p_cyt'],measur[,'STAT5ptot_cyt'],(x0[1]-(measur[,'STAT5ptot_cyt']/parameters['s2']))/2/(1400/450))
                                                                                                                                                                       

y<- as.data.frame(y)
# ---


y[y<0]                       <- 0
colnames(y)                  <- c("time", "STAT5" ,"STAT5p_cyt","STAT5ptot_cyt","STAT5_n")


modelJakStat  <- function(t, x, parameters, input) {
  with (as.list(parameters),{
    
    k1 = parameters[1]
    k2 = parameters[2]
    k3 = parameters[3]
    
    
    
    u <- input$u(t)
    
    dx1 = -k1 * x[1]  * u
    dx2 = k1 *  x[1]  * u - k2 * x[2]^2
    dx3 = -k3*x[3] + 0.5*k2*x[2]*x[2]
    dx4 = k3 * x[3]
    
    
    list(c(dx1 ,dx2 ,dx3 ,dx4 ))
  })
}

objectiveJakStat  <- function(index,y,parameter){
  
  if (index == 1){
    return(y[1])
  }
  if (index == 2){
    return(parameter[5] * (y[2] + 2 * y[3]))
  }
  if (index == 3){
    return(parameter[6] * (y[1] + y[2] + 2 * y[3]))
  }
  if (index == 4){
    return(y[4])
  }
  
}

JakStatModel <- odeModel(func = modelJakStat, parms = parameters, input = inputData, 
                         measFunc = objectiveJakStat, y = x0, meas = y, sd = sd,custom=TRUE)


A <- BDEN(odeModel               = JakStatModel,
          settings               = SETTINGS,
          mcmc_component         = MCMC_component,
          loglikelihood_func     = LOGLIKELIHOOD_func,
          gibbs_update           = GIBBS_update,
          ode_sol                = ode_solv,
          lambda            = .001,
          beta_init         = c(1,1,10,1),
          numbertrialsstep = 10,
          numbertrialseps  = 80)



