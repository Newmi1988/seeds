#' Results Class for the Algorithms
#' 
#' A S4 class that collects the results of the two algorithms
#' 
#' @slot stateNominal data.frame containing the states of the nominal model
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
    stateNominal = "data.frame",
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
    stateNominal = data.frame(),
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

  if(!missing(y)){
    annoX <- y[[1]]
    annoY <- y[[2]]
    if(length(annoX) != length(names(seedsobj@stateEstimates[,-1]))){
      stop('Number of names has to be be equal to the number of states')
    }
    
    if(length(annoY) != length(names(seedsobj@outputEstimates[,-1]))){
      stop('Length of measurement annotations character vector has to be equal to the number of states')
    }
    
  } else {
    annoX <- names(seedsobj@stateEstimates[,-1])
    annoY <- names(seedsobj@outputEstimates[,-1])
  }
  
  labelX <- function(ls,val) {
    dic <- list(annoX)
    return(dic[val])
  }
  
  labelY <- function(ls,val) {
    dic <- list(annoY)
    return(dic[val])
  }

  
  # added formating for plotting the states in the right order
  reformatOrder <- function(df,annoT){
    df$facet = factor(df$state, levels = as.character(unique(factor(df$state))))
    
    return(df)
  }
  
  plot1 <- ggplot2::ggplot(reformatOrder(tidyr::gather(seedsobj@stateEstimates,state, value, -1)), ggplot2::aes(x=t, y=value, colour='red'))+ 
    ggplot2::geom_line()+
    ggplot2::geom_line(data = reformatOrder(tidyr::gather(seedsobj@stateNominal,state, value, -1)), ggplot2::aes(x=t, y=value, colour='blue'))+ 
    ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(seedsobj@stateUnscertainUpper, state, value, -1), tidyr::gather(seedsobj@stateUnscertainLower, state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
    ggplot2::labs(x= 't', y='value',color = "states")+
    ggplot2::scale_color_manual(breaks= c("red","blue"), labels = c("estimate","nominal"), values = c("blue","red"))+
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())+
    ggplot2::facet_wrap(~facet, labeller = labelX)

  
  
  plot2 <- ggplot2::ggplot(data=reformatOrder(tidyr::gather(seedsobj@hiddenInputEstimates,state,value,-1)), ggplot2::aes(x=t,y=value, colour="red"))+
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(seedsobj@hiddenInputUncertainUpper, state, value, -1), tidyr::gather(seedsobj@hiddenInputUncertainLower, state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
    ggplot2::theme(legend.position = "none",
                   strip.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())+
    ggplot2::facet_wrap(~facet)
  
  plot3 <- ggplot2::ggplot(data=reformatOrder(tidyr::gather(seedsobj@outputEstimates,state,value, -1)), ggplot2::aes(x=t, y=value, colour=state))+
    ggplot2::geom_line()+
    ggplot2::geom_errorbar(data=reformatOrder(dplyr::inner_join(tidyr::gather(seedsobj@Data, state,value, -1), tidyr::gather(seedsobj@DataError,state,value, -1), by=c("t","state"))), ggplot2::aes(x=t, ymin= value.x - value.y , ymax = value.x + value.y), inherit.aes = FALSE)+
    ggplot2::geom_point(data =reformatOrder(tidyr::gather(seedsobj@Data,state,value,-1)), ggplot2::aes(x=t, y=value), colour="black") +
    ggplot2::labs(x= 't', y='value',color = "Estimated \n measurements")+
    ggplot2::scale_color_manual(labels = labels(seedsobj@outputEstimates[,-1])[[2]], values = rep("red",length(labels(seedsobj@outputEstimates[,-1])[[2]])))+
    ggplot2::theme(legend.position = "none",
                  strip.background = ggplot2::element_blank(),
                  panel.background = ggplot2::element_blank(),
                  panel.border = ggplot2::element_blank(),
                  panel.grid.major = ggplot2::element_blank(),
                  panel.grid.minor = ggplot2::element_blank())+
    ggplot2::facet_wrap(~facet, drop = TRUE, labeller = labelY)
  
  return(list(plot1,plot2, plot3))
  
}

#' Plot method for the S4 class resultsSeeds
#' 
#' A standardized plot function to display the results of the algorithms. Both
#' algorithms should result in objects of the class resultsSeeds. The results can
#' be plotted using the \code{\link{plot}}-function.
#' 
#' @param x an object of type resultsSeeds or a list of these objects. If a list
#' is given the last entry will be plotted.
#' @param y ...
#' 
#' @rdname resultsSeeds-methods
#' 
#' @aliases plot,resultsSeeds,missing-method
#' 
#' @export
setMethod(f = "plot",
          signature = c(x="resultsSeeds",y="missing"),
          definition = function(x,y)
          {
            plotList <- plotResultsSeeds(x,y)
            
            return(plotList)
          }
)

setMethod(f = "plot",
          signature = c(x="list",y="missing"),
          definition = function(x,y)
          {
            x <- x[[length(x)]]
            plotList <- plotResultsSeeds(x,y)
            
            return(plotList)
          }
)


 
setGeneric(name = "plotAnno", function(x,stateAnno,measAnno) standardGeneric("plotAnno"))

#' Create annotated plot
#' 
#' Create a annotated plot with given state and measurement names. The plots are
#' equal to the output of the normal plot function.
#' @param x an object of type resultsSeeds which contains the results of the algorithms
#' @param stateAnno a character vector describing the names of the states
#' @param measAnno a character vector describing the names of the measurements
#' 
#' @rdname resultsSeeds-methods
#' 
#' @export
setMethod(f = "plotAnno",
          signature = "resultsSeeds",
          definition = function(x,stateAnno,measAnno)
          {
            y <- list(stateAnno,measAnno)
            plotList <- plotResultsSeeds(x,y)
            
            return(plotList)
          }
)

setMethod(f = "plotAnno",
          signature = "list",
          definition = function(x,stateAnno,measAnno)
          {
            x <- x[[length(x)]]
            y <- list(stateAnno,measAnno)
            plotList <- plotResultsSeeds(x,y)
            
            return(plotList)
          }
)



