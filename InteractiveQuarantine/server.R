
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

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

shinyServer(function(input, output) {

  # Plot 1
  slider <- input_slider(min=1, max=5, step=0.01, label = "Reproductive Number (+/- 0.5)")
  data_master %>% 
    ggvis(x = ~R_s, y = ~R_q, fill = ~disease, stroke = ~disease, opacity := 0.3) %>%
    filter(R_0 <= (eval(slider)+0.5)) %>%
    filter(R_0 >= (eval(slider)-0.5)) %>%
    filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
    filter(disease %in% eval(input_checkboxgroup(c("Ebola" = "Ebola", "HepatitisA" = "HepatitisA", "InfluenzaA" = "InfluenzaA", "MERS"="MERS", "Pertussis"="Pertussis", "SARS"="SARS", "Smallpox"="Smallpox"), select = "Ebola"))) %>%
#     layer_rects(x = 0, x2 = 1, y = 1, y2 = 4, opacity := 0.1, fill = "yellow") %>%
#       layer_rects(x = 1, x2 = 4, y = 0, y2 = 1, opacity := 0.1, fill = "blue") %>%
#       layer_rects(x = 0, x2 = 1, y = 0, y2 = 1, opacity := 0.1, fill = "green") %>%
    layer_points(stroke := NA) %>%
    # layer_points(fill := NA) %>%
    layer_points(x = eval(slider), y = eval(slider), fill := "grey", shape := "cross") %>%
    #   layer_text(x = eval(slider), y = eval(slider), text := "..R") %>%
    scale_numeric("x", domain = c(0, 5), nice = FALSE, label = "Symptom Monitoring") %>%
    scale_numeric("y", domain = c(0, 5), nice = FALSE, label = "Quarantine") %>%
    bind_shiny("plot1", "plot1_ui")
    
    # Plot 2
    slider2 <- input_slider(min=0, max=5, c(0,5),step = 0.1, label = "Y axis upper limit")
    data_master %>%
      ggvis(~R_0,  input_select(c("R_s", "R_q", "Abs_Benefit","Rel_Benefit", "Abs_Benefit_per_Qday", "Rel_Benefit_per_Qday"), map=as.name, label = "Outcome"), fill = ~disease, stroke = ~disease) %>%
      group_by(disease) %>% 
      filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
      filter(disease %in% eval(input_checkboxgroup(c("Ebola" = "Ebola", "HepatitisA" = "HepatitisA", "InfluenzaA" = "InfluenzaA", "MERS"="MERS", "Pertussis"="Pertussis", "SARS"="SARS", "Smallpox"="Smallpox"), select = "Ebola"))) %>%
      layer_points(opacity := 0.1) %>%
      scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
      scale_numeric("y", domain = slider2, clamp=TRUE, label = "Outcome") %>%
      bind_shiny("plot2", "plot2_ui")
    
})
