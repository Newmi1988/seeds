createCLangRoot <- function(rootStates) {
  #' A function for creating a root function that can be used with the ode solver lsoda for triggering events - C Version
  #' @param rootStates Vector of the same length as the state, ones and zeros indicated which states can start an event

  head <- '\n\nvoid myroot(int *neq, double *t, double *y, int *ng, double *gout,double *out, int *ip ){\n\n'

  stateId <- which(rootStates > 0)

  body <- paste0('\tgout[', 0:(sum(rootStates) - 1), '] = y[', stateId - 1, '];')
  formatedBody <- paste0(body, collapse = '\n')

  rootFunction <- paste0(head, formatedBody, '\n}')

  return(rootFunction)

}

createRoot <- function(rootStates) {
  #' A function for creating a root function that can be used with the ode solver lsoda for triggering events
  #' @param rootStates Vector of the same length as the state, ones and zeros indicated which states can start an event

  rootFunction <- 'function(t,x,param) x'

  return(rootFunction)
}

createEvent <- function(tolerance, value) {
  #' Create an event function that sets states that are zero to a small value
  #' @param tolerance tolerance under which a state is set to a given value 'value'
  #' @param value value the states should be set to, very small values result in a much longer computation time

  head <- 'function(t,x,param) {\n\n'

  body <- paste0('\t x[which(x <= ', tolerance, ')] = ', value)

  eventFunction <- paste0(head, body, '\n\n\treturn(x)\n}')

  return(eventFunction)

}
