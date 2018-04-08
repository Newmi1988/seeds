shiny::shinyUI(shiny::pageWithSidebar(
  
  shiny::headerPanel(title = "Seeds Results"),
  
  shiny::sidebarPanel(
    
    shiny::selectInput("dataset","Chose a dataset:",
                        choices = c("states","hidden inputs","measurements")),
    
    shiny::selectInput("variable","select variables:","", multiple = TRUE),
    
    shiny::checkboxInput("hideZero","Hide trajectories equal to zero", FALSE)
  ),
  
  shiny::mainPanel(
    shiny::h3(shiny::textOutput("caption")),
    
    shiny::plotOutput("statePlot")
  )
  
  )
  
)
