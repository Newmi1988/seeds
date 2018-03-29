plot.resultsSeeds  <- function(obj) {
  
  reformatOrder <- function(df){
    df$facet = factor(df$state, levels = as.character(unique(factor(df$state))))
    
    return(df)
  }
  
  
  plot1 <- ggplot2::ggplot(reformatOrder(tidyr::gather(obj@stateEstimates,state, value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+ 
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(obj@stateUnscertainUpper, state, value, -1), tidyr::gather(obj@stateUnscertainLower, state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
    ggplot2::scale_color_discrete()+
    ggplot2::facet_wrap(~facet)
  # if(!(is.data.frame(obj@stateUnscertainUpper) && nrow(obj@stateUnscertainUpper)==0)){
  #   plot1 + 
  # }

  
  plot2 <- ggplot2::ggplot(data=reformatOrder(tidyr::gather(obj@hiddenInputEstimates,state,value,-1)), ggplot2::aes(x=t,y=value, colour="red"))+
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(obj@hiddenInputUncertainUpper, state, value, -1), tidyr::gather(obj@hiddenInputUncertainLower, state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
    ggplot2::scale_color_discrete()+
    ggplot2::facet_wrap(~facet)

  
  #print output
  plot3 <- ggplot2::ggplot(data=reformatOrder(tidyr::gather(obj@outputEstimates,state,value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+
    ggplot2::geom_line()+
    ggplot2::geom_errorbar(data=reformatOrder(dplyr::inner_join(tidyr::gather(obj@Data, state,value, -1), tidyr::gather(obj@DataError,state,value, -1), by=c("t","state"))), ggplot2::aes(x=t, ymin= value.x - value.y , ymax = value.x + value.y), inherit.aes = FALSE)+
    ggplot2::geom_point(data =reformatOrder(tidyr::gather(obj@Data,state,value,-1)), ggplot2::aes(x=t, y=value), colour="black") +
    ggplot2::scale_color_discrete()+
    ggplot2::facet_wrap(~facet)
  
  
  

  
  return(list(plot1,plot2, plot3))

}

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

