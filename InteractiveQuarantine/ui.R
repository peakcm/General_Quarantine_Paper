
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Check out http://shiny.rstudio.com/gallery/movie-explorer.html
# Check out https://github.com/hrbrmstr/ggvis-maps
require(shiny)
require(ggvis)
require(dplyr)
require(magrittr)
require(foreign)
require(RCurl)

url <- "https://raw.githubusercontent.com/peakcm/InteractiveQuarantine/master/20170328_FigureRsRq.csv"
getURL <- getURL(url) 
data_master <- read.csv(textConnection(getURL))

data_master$disease <- factor(data_master$disease, levels = c("Ebola", "HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox"), labels = c("Ebola", "HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox"), ordered = TRUE)

shinyUI(fluidPage(

  # Application title
  titlePanel("Interactive Supplement"),

  sidebarLayout(
    
    sidebarPanel(
      h6("This website is an interactive supplement to:"),
      h6(a("Comparing nonpharmaceutical interventions for containing emerging epidemics", href="http://www.pnas.org/content/early/2017/03/27/1616438114")),
      h6("by Corey M Peak, Lauren M Childs, Yonatan H Grad, and Caroline O Buckee"),
      h6("Please contact peak@mail.harvard.edu for comments and issues."),
      # uiOutput("plot1_ui"),
      # uiOutput("plot2_ui"),
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
