rm(list=ls())
devtools::load_all()
graphics.off()
# library('deSolve')
# 
# testmodell
example <- function(t,x, parms) {
  with(as.list( c(x, parms)),{
    dx1 = 1
    dx2 = a * cos(t) * x[1]
    dx3 = a * sin(t) * x[1]

    list(c(dx1,dx2, dx3))
  })
}

t           <- seq(0.1,4.1,1)
parameters  = c(a = 1, b = 1)
x0          = c(x1 = 0, x2 = 1, x3 = 1)

# indicator which states should be non negative
nnStates = c(0,1,0)

# 'createEvent' erstellt automatisch eine event-funktion, die durch eval einer aus text einer variable zugewiesen werden kann
# der zweite parameter gibt den wert an auf den die funktion gesetzt werden soll

# model  defininieren und eventfunction erstellen
model <- odeModel(func = example, parms = parameters, times = t, y = x0, nnStates = nnStates, resetValue = 0.0001)
myEvent <- eval(parse(text = createEvent(tollerance = 0., value = 0.0001)))

ext <- .Platform$dynlib.ext
compiledModel <- paste0('model',ext)


if(is.loaded('derivsc')){
  dyn.unload(compiledModel)
}
createCompModel(modelFunc = model@func,parameters = model@parms, nnStates = nnStates)
system("R CMD SHLIB model.c")
dyn.load(compiledModel)

## dummy inputs - nacher selber setzen
input <- rep(0,length(t))
uList = list(cbind(t,input)) 
w0 <- rep(0, length(t))
w1 <- rep(0, length(t))
w1[as.integer(length(w1)/2+2):length(w1)] = 3

w <- matrix(c(w0,w1,w0), ncol = length(x0))
wSplit <- split(w, rep(1:ncol(w), each = nrow(w)))
wList <- lapply(wSplit, FUN = function(x) cbind(t,x))
forcings <- c(uList, wList)


## lösung
resConstant <- deSolve::lsoda(y = model@y, times = t, func = "derivsc",
                         parms = model@parms, dllname = "model", initforc="forcc",
                         forcings = forcings, initfunc = "parmsc", nroot = sum(model@nnStates),
                         rootfunc = "myroot", events = list(func = myEvent, root = TRUE),
                         fcontrol = list( method = "constant", f=0.0, rule = 2))

plot(resConstant)
max(resConstant)

resLinear <- deSolve::lsoda(y = model@y, times = t, func = "derivsc",
                         parms = model@parms, dllname = "model", initforc="forcc",
                         forcings = forcings, initfunc = "parmsc", nroot = sum(model@nnStates),
                         rootfunc = "myroot", events = list(func = myEvent, root = TRUE),
                         fcontrol = list( method = "linear", f=0.0, rule = 2))

plot(resLinear)
max(resLinear)
