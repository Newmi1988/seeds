rm(list=ls())
devtools::load_all()
graphics.off()

modelJakStat  <- function(t, x, parameters, input) {
  with (as.list(parameters),{
    
    k1 = parameters[1]
    k2 = parameters[2]
    k3 = parameters[3]
    k4 = parameters[4]
    s1 = parameters[5]
    s2 = parameters[6]
    
    u <- input$u(t)
    
    dx1 = -k1 * x[1]  * u
    dx2 = k1 *  x[1]  * u - k2 * x[2]*x[2]
    dx3 = -k3*x[3] + 0.5*k2*x[2]*x[2]
    dx4 = k3 * x[3]
    
    list(c(dx1 ,dx2 ,dx3 ,dx4 ))
  })
}

N = 10^0.31
x0 = c(N, 0, 0, 0)
y <- c(X = x0)
evalTimes <- c( 0.0208,  0.1098,   0.2696,    0.4999,    0.8002,    1.1697,    1.6077,    2.1129,    2.6843,    3.3205,    4.0200,    4.7811,    5.6020,    6.4808,    7.4154,    8.4035,    9.4429,   10.5310, 11.6653,   12.8431,   14.0616,   15.3179,   16.6090,   17.9319,   19.2834,   20.6603,   22.0594,   23.4773,   24.9107,   26.3561,   27.8102,   29.2695,   30.7305,   32.1898,   33.6439,   35.0893, 36.5227,   37.9406,   39.3397,   40.7166,   42.0681,   43.3910,   44.6821,   45.9384,   47.1569,   48.3347,   49.4690,   50.5571,   51.5965,   52.5846,   53.5192,   54.3980,   55.2189,   55.9800, 56.6795,   57.3157,   57.8871,   58.3923,   58.8303,   59.1998,   59.5001,   59.7304,   59.8902,   59.9792)
parameters = 10^c("k1"=0.31, "k2"=-1, "k3"=-0.49, "k4"= 0.42, "s1"=-0.21, "s2"=-0.34)


dyn.unload('model.dll')
## Function welche die C files erstellt.
createCompModel(modelFunc = modelJakStat, parameters = parameters)

ext <- .Platform$dynlib.ext
compiledModel <- paste0('model',ext)

# unload dynlib if is already loaded make rewriting enable through unloading the dll
if(is.loaded('derivsc')){
  dyn.unload(compiledModel)
}
# kompilieren des modells
system("R CMD SHLIB model.c")
# laden der dynamischen library
dyn.load('model.dll')

inputData <- read.table('http://jeti.uni-freiburg.de/PNAS_Swameye_Data/DATA1_hall_inp')
inputData[nrow(inputData),2] = 0.009

times <- inputData[,1]

## Beispiel mit keinen Input auf den ws
# w werden alle mit 0 initialisiert
w <- matrix(rep(0,length(x0)*length(times)), ncol = length(x0))
# aufteilung der hidden input Matrix in einzelene vectoren; Spalten sind die Werte von w_i zu bestimmten Zeitpunkten
wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
# die einzelnen Forcings werden als liste an deSolve::ode übergeben
# ist ein externer Input gegeben (wie bei JakStat u) so muss der als erstes in der Liste stehen
# Einträge der liste sind Matrizen

# Inputs können linear interpoliert werden, hier behalten sie jedoch die Zeitpunkte des Inputfiles
inputApprox <- apply(X = inputData[,-1, drop=F], MARGIN = 2, FUN = function(x) stats::approx(x = inputData[,1], y = x, xout = times, rule = 2))
inputApprox = list(cbind(times,inputApprox$V2$y))

wList <- lapply(wSplit, FUN = function(x) cbind(times,x))
forcings <- c(inputApprox, wList)

solNominal <- deSolve::ode(y = x0, evalTimes, func = "derivsc",
                          parms = parameters, dllname = "model", initforc="forcc",
                          forcings = forcings, initfunc = "parmsc")

solNominal
plot(solNominal)
