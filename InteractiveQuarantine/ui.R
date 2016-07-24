
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Check out http://shiny.rstudio.com/gallery/movie-explorer.html
# Check out https://github.com/hrbrmstr/ggvis-maps
library(shiny)
library(ggvis)
library(dplyr)
library(magrittr)
library(foreign)
library(RCurl)

url <- "https://raw.githubusercontent.com/peakcm/InteractiveQuarantine/master/20151113_FigureRsRq.csv"
getURL <- getURL(url) 
data_master <- read.csv(textConnection(getURL))

data_master$disease <- factor(data_master$disease, levels = c("Ebola", "HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox"), labels = c("Ebola", "HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox"), ordered = TRUE)

shinyUI(fluidPage(

  # Application title
  titlePanel("Interactive Supplement"),

  sidebarLayout(
    
    sidebarPanel(
      h6("This website is a <draft> interactive supplement to\n
'Containing Emerging Epidemics: a Quantitative Comparison of Quarantine and Symptom Monitoring'\n
by Corey M Peak, Lauren M Childs, Yonatan H Grad, and Caroline O Buckee.\n
         Please contact peak@mail.harvard.edu for issues and comments."),
      # uiOutput("plot1_ui"),
      # uiOutput("plot2_ui"),
      width = 3
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Figure 3A", wellPanel(p("Interactive Supplement to Figure 3A")), 
                 uiOutput("plot1_ui"),
                 ggvisOutput("plot1")),
        tabPanel("Additional Model Outputs", wellPanel(p("Interactive Supplement for Intermediate Outputs from Model Simulations")),
                 uiOutput("plot2_ui"),
                 ggvisOutput("plot2"))
    ),
    width = 9
  )
  )
  
))
