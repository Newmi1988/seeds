% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resultsSeeds.R
\docType{class}
\name{resultsSeeds-class}
\alias{resultsSeeds-class}
\alias{resultsSeeds}
\title{Results Class for the Algorithms}
\value{
A object of class resultsSeeds collecting all the results of the algorithm
}
\description{
A S4 class that collects the results of the two algorithms. The class also is equipped
with functions for easily plotting and extracting the different results.
}
\section{Slots}{

\describe{
\item{\code{stateNominal}}{data.frame containing the states of the nominal model}

\item{\code{stateEstimates}}{data.frame containing the state estimates}

\item{\code{stateUnscertainLower}}{lower bound of the estimated states as calculated by the baysian method}

\item{\code{stateUnscertainUpper}}{upper bound of the estimated states as calculated by the baysian method}

\item{\code{hiddenInputEstimates}}{estimated hidden input}

\item{\code{hiddenInputUncertainLower}}{lower bounds of the estimated hidden inputs}

\item{\code{hiddenInputUncertainUpper}}{upper bounds of the estimated hidden inputs}

\item{\code{outputEstimates}}{estimated measurements resulting from the control of the hidden inputs}

\item{\code{outputEstimatesUncLower}}{lower bound of the confidence bands of the estimated output}

\item{\code{outputEstimatesUncUpper}}{upper bound of the confidence bands of the estimated output}

\item{\code{Data}}{the given measurements}

\item{\code{DataError}}{standard deviation of the given measurements}
}}

