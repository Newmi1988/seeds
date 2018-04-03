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

parameters = 10^c("k1"=0.31, "k2"=-1, "k3"=-0.49, "k4"= 0.42, "s1"=-0.21, "s2"=-0.34)
# you can used a unnamed vector for the parameters if these are declared like in the function above
parameters = 10^c(0.31, -1, -0.49, 0.42, -0.21,-0.34)

## Function welche die C files erstellt.
createCompModel(modelFunc = modelJakStat, parameters = parameters)
# c file darstellen
file.edit('model.c')
