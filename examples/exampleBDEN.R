#### clear workspace, load package, remove graphics ####

rm(list=ls())
devtools::load_all()
graphics.off()

#### setze initialen zustand, parameter und Zeitpunkte der Auswertung durch deSolve ####
#N = 10^0.31
N = 2.4290
x0 = c(N, 0, 0, 0)


# parameter des models
parameters = c("k1"=2.4290, "k2"=975.4280, "k3"=0.1157, "k4"= 0, "s1"=10^-0.21, "s2"=10^-0.34)
#parameters = 10^c("k1"=0.31, "k2"=-1, "k3"=-0.49, "k4"= 0.42, "s1"=-0.21, "s2"=-0.34)


evalTimes <- c( 0.0208,  0.1098,   0.2696,    0.4999,    0.8002,    1.1697,    1.6077,    2.1129,    2.6843,    3.3205,    4.0200,    4.7811,    5.6020,    6.4808,    7.4154,    8.4035,    9.4429,   10.5310, 11.6653,   12.8431,   14.0616,   15.3179,   16.6090,   17.9319,   19.2834,   20.6603,   22.0594,   23.4773,   24.9107,   26.3561,   27.8102,   29.2695,   30.7305,   32.1898,   33.6439,   35.0893, 36.5227,   37.9406,   39.3397,   40.7166,   42.0681,   43.3910,   44.6821,   45.9384,   47.1569,   48.3347,   49.4690,   50.5571,   51.5965,   52.5846,   53.5192,   54.3980,   55.2189,   55.9800, 56.6795,   57.3157,   57.8871,   58.3923,   58.8303,   59.1998,   59.5001,   59.7304,   59.8902,   59.9792)

inputData                    <- read.table('http://jeti.uni-freiburg.de/PNAS_Swameye_Data/DATA1_hall_inp')
inputData[nrow(inputData),2] <- 0.009
colnames(inputData)          <- c('t','u')
measure                      <- read.table('http://jeti.uni-freiburg.de/PNAS_Swameye_Data/DATA1_hall')
colnames(measure)            <- c("time","STAT5p_cyt" ,"sd_STAT5p_cyt","STAT5ptot_cyt","sd_STAT5ptot_cyt")

sd                           <- cbind(measure['time'],measure['sd_STAT5p_cyt'],measure['sd_STAT5ptot_cyt'])

y                            <- cbind(measure['time'],((measure['STAT5ptot_cyt']/parameters['s2'])-(measure['STAT5p_cyt']/parameters['s1'])),measure['STAT5p_cyt'],measure['STAT5ptot_cyt'],(x0[1]-(measure['STAT5ptot_cyt']/parameters['s2']))/2/(1400/450))
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


measJakStat <- function(x,index,parameter) {
  
  s1 <- 10^(-0.21)
  s2 <- 10^(-0.34)
  
  y1 = s1*(x[,2]+ 2*x[,3])
  y2 = s2*(x[,1] + x[,2] + 2*x[,3])
  
  return(list(y1,y2))
}


# in BDEN loglikelhood
objectiveJakStat  <- function(index,y,parameter){
  
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



A <- BDEN(measData           = y,
     x0                      = x0,
     parameters              = parameters,
     systemInput             = inputData,
     sd                      = sd,
     model                   = modelJakStat,
     measFunc                = objectiveJakStat,
     settings                = SETTINGS,
     mcmc_component          = MCMC_component,
     loglikelihood_func      = LOGLIKELIHOOD_func,
     gibbs_update            = GIBBS_update,
     ode_sol                 = ode_solv)









