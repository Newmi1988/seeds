#' Results Class for the Algorithms
#' 
#' A S4 class that collects the results of the two algorithms
#' 
#' @slot stateEstimates data.frame containing the state estimates
#' @slot stateUnscertainLower lower bound of the estimated states as calculated by the baysian method
#' @slot stateUnscertainUpper upper bound of the estimated states as calculated by the baysian method
#' @slot hiddenInputEstimates estimated hidden input
#' @slot hiddenInputUncertainLower lower bounds of the estimated hidden inputs
#' @slot hiddenInputUncertainUpper uppper bounds of the estimated hidden inputs
#' @slot outputEstimates estimated measurements resulting from the control of the hidden inputs
#' @slot Data the given measurements
#' @slot DataError standard deviation of the given measurements
#' 
#' @export resultsSeeds
#' @exportClass resultsSeeds
#' 
#' @import methods
#' @importFrom graphics plot
resultsSeeds <- setClass(
  'resultsSeeds',
  slots = c(
    stateEstimates = "data.frame",
    stateUnscertainLower = "data.frame",
    stateUnscertainUpper = "data.frame",
    hiddenInputEstimates = "data.frame",
    hiddenInputUncertainLower = "data.frame",
    hiddenInputUncertainUpper = "data.frame",
    outputEstimates = "data.frame",
    Data = "data.frame",
    DataError = "data.frame"
  ),
  prototype = c(
    stateEstimates = data.frame(),
    stateUnscertainLower = data.frame(),
    stateUnscertainUpper = data.frame(),
    hiddenInputEstimates = data.frame(),
    hiddenInputUncertainLower = data.frame(),
    hiddenInputUncertainUpper = data.frame(),
    outputEstimates = data.frame(),
    Data = data.frame(),
    DataError = data.frame()
  )
)


plotResultsSeeds  <- function(x,y) {
  
  seedsobj = x
  
  # added formating for plotting the states in the right order
  reformatOrder <- function(df){
    df$facet = factor(df$state, levels = as.character(unique(factor(df$state))))
    
    return(df)
  }
  
  plot1 <- ggplot2::ggplot(reformatOrder(tidyr::gather(seedsobj@stateEstimates,state, value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+ 
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(seedsobj@stateUnscertainUpper, state, value, -1), tidyr::gather(seedsobj@stateUnscertainLower, state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
    ggplot2::scale_color_discrete()+
    ggplot2::facet_wrap(~facet)
  
  
  plot2 <- ggplot2::ggplot(data=reformatOrder(tidyr::gather(seedsobj@hiddenInputEstimates,state,value,-1)), ggplot2::aes(x=t,y=value, colour="red"))+
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(seedsobj@hiddenInputUncertainUpper, state, value, -1), tidyr::gather(seedsobj@hiddenInputUncertainLower, state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
    ggplot2::scale_color_discrete()+
    ggplot2::facet_wrap(~facet)
  
  
  plot3 <- ggplot2::ggplot(data=reformatOrder(tidyr::gather(seedsobj@outputEstimates,state,value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+
    ggplot2::geom_line()+
    ggplot2::geom_errorbar(data=reformatOrder(dplyr::inner_join(tidyr::gather(seedsobj@Data, state,value, -1), tidyr::gather(seedsobj@DataError,state,value, -1), by=c("t","state"))), ggplot2::aes(x=t, ymin= value.x - value.y , ymax = value.x + value.y), inherit.aes = FALSE)+
    ggplot2::geom_point(data =reformatOrder(tidyr::gather(seedsobj@Data,state,value,-1)), ggplot2::aes(x=t, y=value), colour="black") +
    ggplot2::scale_color_discrete()+
    ggplot2::facet_wrap(~facet)
  
  return(list(plot1,plot2, plot3))
  
}

#' @export
setMethod(f = "plot",
          signature = c(x="resultsSeeds",y="missing"),
          definition = function(x,y)
          {
            plotList <- plotResultsSeeds(x,y)
            
            return(plotList)
          }
)



