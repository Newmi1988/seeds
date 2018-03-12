N = 10^0.31
x0 = c(N, 0, 0, 0)
y <- c(X = x0)
times <- c( 0.0208,  0.1098,   0.2696,    0.4999,    0.8002,    1.1697,    1.6077,    2.1129,    2.6843,    3.3205,    4.0200,    4.7811,    5.6020,    6.4808,    7.4154,    8.4035,    9.4429,   10.5310, 11.6653,   12.8431,   14.0616,   15.3179,   16.6090,   17.9319,   19.2834,   20.6603,   22.0594,   23.4773,   24.9107,   26.3561,   27.8102,   29.2695,   30.7305,   32.1898,   33.6439,   35.0893, 36.5227,   37.9406,   39.3397,   40.7166,   42.0681,   43.3910,   44.6821,   45.9384,   47.1569,   48.3347,   49.4690,   50.5571,   51.5965,   52.5846,   53.5192,   54.3980,   55.2189,   55.9800, 56.6795,   57.3157,   57.8871,   58.3923,   58.8303,   59.1998,   59.5001,   59.7304,   59.8902,   59.9792)
parameters = 10^c(0.31, -1, -0.49, 0.42, -0.21, -0.34)

inputData <- jakstatInput
measure <- jakstatMeasurement


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
    dx2 = k1 *  x[1]  * u - k2 * x[2]^2
    dx3 = -k3*x[3] + 0.5*k2*x[2]^2
    dx4 = k3 * x[3]

    list(c(dx1 ,dx2 ,dx3 ,dx4 ))
  })
}

measJakStat <- function(x) {

  s1 <- 10^(-0.21)
  s2 <- 10^(-0.34)

  y1 = s1*(x[,2]+ 2*x[,3])
  y2 = s2*(x[,1] + x[,2] + 2*x[,3])

  return(list(y1,y2))
}

y <- data.frame(measure$t, measure$y1, measure$y2)
sd <- data.frame(measure$y1sd, measure$y2sd)

JakStatConst <- '2*x4+ 2*x3 + x1 + x2 == N'


results <- greedyApproach(alphaStep = 0.01, alpha2 = 0.2, Beta = 0.8,
               x0 = x0, optW = c(1,1,1,1) , times=times,
               measFunc= measJakStat,  measData = y, std = sd,
               parameters = parameters, systemInput = inputData,
               modelFunc = modelJakStat, plotEstimates = TRUE, conjGrad = FALSE, cString = JakStatConst)
