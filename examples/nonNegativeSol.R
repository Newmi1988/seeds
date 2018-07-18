## ode system with negative solutions

library(deSolve) 

# inputs 
T = 0  + 273.15     # K (Kelvin ) / tTemperature 
V = 50              # mL / Volume
A_solid = 125/V     # mmol/mL = mol/L / initial concentration of product A_solid
B = 100/V           # mol/L / initial concentration of product B

# parameters
R = 8.314           # J/(K*mol) / gas constant
expfact_sat = 2
E_a_sat = 10^3

params <- c(k1 = 0.1, # rate constant of dissolution
            k2 = 3*10^(-3), # rate constant of reaction
            A_sat = expfact_sat*exp(-E_a_sat/(R*T))) # saturation of the A_bulk into the solvent

# initial values
state <- c(A_solid = A_solid, A_bulk = 0, B = B, C = 0)

# system of differential equations
derivs <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dA_solid = -k1*(A_sat - A_bulk)
    dA_bulk = -k2*A_bulk*B + k1*(A_sat - A_bulk)
    dB = -k2*A_bulk*B
    dC = k2*A_bulk*B
    return(list(c(dA_solid, dA_bulk, dB, dC)))
  })                                                        
}

times = seq(0, 500, by = 0.01)
init <- ode(y = state, func = derivs, time = times, parms = params)

l = dim(init)[1]-1
matplot(init[,1], init[,-1], type = "l", lty = 1:1, lwd = c(2),
        col = 1:l, xlab = "time [min]", ylab = "concn [mol/L]")
legend("topright", colnames(init)[-1], col = 1:l, lwd = c(2))

### rewrite
logDeriv <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dx1 = (-k1*(A_sat - x2)) 
    dx2 = -k2*x2*x3 + k1*(A_sat - x2)
    dx3 = -k2*x2*x3
    dx4 = k2*x2*x3
    return(list(c(dx1, dx2, dx3, dx4)))
  })                                                        
}

state <- c(x1 = log(A_solid), x2 = 0, x3 = B, x4 = 0)
times = seq(0, 500, by = 0.01)
init <- as.data.frame(ode(y = state, func = logDeriv, time = times, parms = params))

indTrans <- c(0,1,0,0,0)
init[,  which(indTrans>0)] = exp(init[,  which(indTrans>0)])

l = dim(init)[1]-1
matplot(init[,1], init[,-1], type = "l", lty = 1:1, lwd = c(2),
        col = 1:l, xlab = "time [min]", ylab = "concn [mol/L]")
legend("topright", colnames(init)[-1], col = 1:l, lwd = c(2))
