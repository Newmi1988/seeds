plot.resultsSeeds  <- function(obj) {
  plot1 <- ggplot(obj@stateEstimates %>% gather(state, value, -1), aes(x=t, y=value, colour="red"))+ 
    geom_line()+ 
    scale_color_discrete()+
    facet_wrap(~state)

  
  #print output
  plot2 <- ggplot(data=obj@outputEstimates %>% gather(state,value, -1), aes(x=t, y=value, colour="red"))+
    geom_line()+
    geom_errorbar(data=obj@Data %>% gather(state,value, -1) %>% inner_join( obj@DataError %>% gather(state,value, -1), by=c("t","state")), aes(x=t, ymin= value.x - value.y , ymax = value.x + value.y), inherit.aes = FALSE)+
    geom_point(data = obj@Data %>% gather(state,value,-1), aes(x=t, y=value), colour="black") +
    scale_color_discrete()+
    facet_wrap(~state)
  
  
  return(list(plot1,plot2))

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

