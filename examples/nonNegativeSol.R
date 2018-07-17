## ode system with negative solutions

library(deSolve) 

# inputs 
T = 0  + 273.15 # K (Kelvin ) / tTemperature 
V = 50 # mL / Volume
A_solid = 125/V # mmol/mL = mol/L / initial concentration of product A_solid
B = 100/V # mol/L / initial concentration of product B

# parameters
R = 8.314 # J/(K*mol) / gas constant
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
