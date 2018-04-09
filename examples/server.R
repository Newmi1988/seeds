data <- res

shiny::shinyServer(function(input,output, session){
  
  dataInput <- shiny::reactive({
    switch (input$dataset,
                      "states" = list(data@stateEstimates, data@stateUnscertainLower, data@stateUnscertainUpper, data@stateNominal),
                      "hidden inputs" = list(data@hiddenInputEstimates, data@hiddenInputUncertainLower, data@hiddenInputUncertainUpper),
                      "measurements" = list(data@outputEstimates, data@Data, data@DataError)
          
    )
  })
  
  captionSwitch <- shiny::reactive({
    paste("estimated ",input$dataset)
  })
  
  shiny::observe({
    vchoice <- names(dataInput()[[1]])
    shiny::updateSelectInput(session,
                             "variable",
                             choices = vchoice[-1])
  })
  
  output$caption <- shiny::renderText({
    captionSwitch()
  })
  
  output$statePlot <- shiny::renderPlot({
    reformatOrder <- function(df){
      df$facet = factor(df$state, levels = as.character(unique(factor(df$state))))
      
      return(df)
    }
    
    data <- dataInput()
    
    if(input$hideZero==TRUE){
      data = lapply(data, '[', colSums(data[[1]])!=0)
    }
    
    #subset for variables
    if(length(input$variable)>0){
      data = lapply(data, function(x) subset.data.frame(x, select = c("t",input$variable)))
    }
    
    if(input$dataset=="measurements"){
      ggplot2::ggplot(data=reformatOrder(tidyr::gather(data[[1]],state,value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+
        ggplot2::geom_line()+
        ggplot2::geom_errorbar(data=reformatOrder(dplyr::inner_join(tidyr::gather(data[[2]], state,value, -1), tidyr::gather(data[[3]],state,value, -1), by=c("t","state"))), ggplot2::aes(x=t, ymin= value.x - value.y , ymax = value.x + value.y), inherit.aes = FALSE)+
        ggplot2::geom_point(data =reformatOrder(tidyr::gather(data[[2]],state,value,-1)), ggplot2::aes(x=t, y=value), colour="black") +
        ggplot2::scale_color_discrete()+
        ggplot2::theme(legend.position="none")+
        ggplot2::facet_wrap(~facet)
    } else if(input$dataset=="hidden inputs") {
      ggplot2::ggplot(reformatOrder(tidyr::gather(data[[1]],state, value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+ 
        ggplot2::geom_line()+
        ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(data[[3]], state, value, -1), tidyr::gather(data[[2]], state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
        ggplot2::scale_color_discrete()+
        ggplot2::theme(legend.position="none")+
        ggplot2::facet_wrap(~facet)
    } else{
      ggplot2::ggplot(reformatOrder(tidyr::gather(data[[1]],state, value, -1)), ggplot2::aes(x=t, y=value, colour="red"))+ 
        ggplot2::geom_line()+
        ggplot2::geom_line(data = reformatOrder(tidyr::gather(data[[4]],state, value, -1)), ggplot2::aes(x=t, y=value, colour="blue"))+ 
        ggplot2::geom_ribbon(data=reformatOrder(dplyr::inner_join(tidyr::gather(data[[3]], state, value, -1), tidyr::gather(data[[2]], state, value, -1), by = c("t","state"))), ggplot2::aes(x= t, ymin = value.y, ymax=value.x), alpha=0.2, inherit.aes = FALSE)+
        ggplot2::scale_color_discrete()+
        ggplot2::theme(legend.position="none")+
        ggplot2::facet_wrap(~facet)
    }

    }
  )
}
)
