#' @details Details
#' 
#' The conjugate gradient based algorithm uses a greedy algorithm to estimate a 
#' sparse controll that tries to minimize the discrepancies between a given 
#' 'nominal model given the measurements (e.g from an experiment). The algorithm
#' uses \pkg{deSolve} to calculate the hidden inputs w based on the adjoint equations
#' of the ODE-System. 
#' 
#' The adjoint equations are calculated using \pkg{Deriv}. For the usage of the 
#' algorithm please look into the examples and documentation given for the
#' functions 
#' \describe{
#'    \item{\code{\link{greedyApproach}}}{a greedy algorithm to calculate a sparse control}
#'    \item{\code{\link{dynElasticNet}}}{a conjugate gradient optimal control algorithm with L2 regulation} 
#' }
"_PACKAGE"
