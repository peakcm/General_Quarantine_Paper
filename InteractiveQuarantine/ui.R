
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# http://shiny.rstudio.com
# Check out http://shiny.rstudio.com/gallery/movie-explorer.html
# Check out https://github.com/hrbrmstr/ggvis-maps

require(ggvis)

shinyUI(fluidPage(

  # Application title
  titlePanel("Interactive Supplement"),

  sidebarLayout(
    
    sidebarPanel(
      h6("This website is an interactive supplement to:"),
      h6(a("Comparing nonpharmaceutical interventions for containing emerging epidemics", href="http://www.pnas.org/content/early/2017/03/27/1616438114")),
      h6("by Corey M Peak, Lauren M Childs, Yonatan H Grad, and Caroline O Buckee"),
      h6("Please contact peak@mail.harvard.edu for comments and issues."),
      width = 3
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Figure 3A", wellPanel(p("Below is an interactive supplement to Figure 3A.")), 
                 uiOutput("plot1_ui"),
                 ggvisOutput("plot1")),
        tabPanel("Additional Model Outputs", wellPanel(p("Below are additional outputs from the transmission model.")),
                 uiOutput("plot2_ui"),
                 ggvisOutput("plot2"))
    ),
    width = 9
  )
  )
  
))
