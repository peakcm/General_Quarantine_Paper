# This is the server logic for a Shiny web application.

require(shiny)
require(ggvis)
require(dplyr)
require(magrittr)
require(foreign)
require(RCurl)

url <- "https://raw.githubusercontent.com/peakcm/InteractiveQuarantine/master/20151113_FigureRsRq.csv"
getURL <- getURL(url)
data_master <- read.csv(textConnection(getURL))

data_master$disease <- factor(data_master$disease, levels = c("Ebola", "HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox"), labels = c("Ebola", "HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox"), ordered = TRUE)

shinyServer(function(input, output) {

  # Plot 1
  slider <- input_slider(min=1, max=5, step=0.01, label = "Basic Reproductive Number (+/- 0.5)")
  data_master %>% 
    ggvis(x = ~R_s, y = ~R_q, fill = ~disease, stroke = ~disease, opacity := 0.3) %>% 
    # add_legend("stroke", title = "Disease(s) Selected")  %>%
    filter(R_0 <= (eval(slider)+0.5)) %>%
    filter(R_0 >= (eval(slider)-0.5)) %>%
    filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("High Intervention Performance Setting" = "HR","Low Intervention Performance Setting" = "LR")))) %>%
    filter(disease %in% eval(input_checkboxgroup(c("Ebola" = "Ebola", "Hepatitis A" = "HepatitisA", "Influenza A" = "InfluenzaA", "MERS"="MERS", "Pertussis"="Pertussis", "SARS"="SARS", "Smallpox"="Smallpox"), select = "Ebola", label = "Select Disease(s)"))) %>%
#     layer_rects(x = 0, x2 = 1, y = 1, y2 = 4, opacity := 0.1, fill = "yellow") %>%
#       layer_rects(x = 1, x2 = 4, y = 0, y2 = 1, opacity := 0.1, fill = "blue") %>%
#       layer_rects(x = 0, x2 = 1, y = 0, y2 = 1, opacity := 0.1, fill = "green") %>%
    layer_points() %>%
    layer_points(x = eval(slider), y = eval(slider), fill := "grey", shape := "cross") %>%
    #   layer_text(x = eval(slider), y = eval(slider), text := "..R") %>%
    scale_numeric("x", domain = c(0, 5), nice = FALSE, label = "Symptom Monitoring") %>%
    scale_numeric("y", domain = c(0, 5), nice = FALSE, label = "Quarantine") %>%
    bind_shiny("plot1", "plot1_ui")
    
    # Plot 2
    slider2 <- input_slider(min=0, max=5, c(0,5),step = 0.1, label = "Define y-axis Limits")
    data_master %>%
      ggvis(~R_0,  input_select(c("Symptom Monitoring Reproductive Number" = "R_s", "Quarantine Reproductive Number" = "R_q","Absolute Difference (Rs-Rq)" = "Abs_Benefit","Relative Difference (Rs-Rq)/Rs" = "Rel_Benefit", "Absolute Difference per Quarantine Day" = "Abs_Benefit_per_Qday", "Relative Difference per Quarantine Day" = "Rel_Benefit_per_Qday"), map=as.name, label = "Select a Model Outcome"), fill = ~disease, stroke = ~disease) %>%
      group_by(disease) %>% 
      filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("High Intervention Performance Setting" = "HR","Low Intervention Performance Setting" = "LR")))) %>%
      filter(disease %in% eval(input_checkboxgroup(c("Ebola" = "Ebola", "Hepatitis A" = "HepatitisA", "Influenza A" = "InfluenzaA", "MERS"="MERS", "Pertussis"="Pertussis", "SARS"="SARS", "Smallpox"="Smallpox"), select = "Ebola", label = "Select Disease(s)"))) %>%
      layer_points(opacity := 0.1) %>%
      scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "Basic Reproductive Number") %>%
      scale_numeric("y", domain = slider2, clamp=TRUE, label = "Selected Outcome") %>%
      # add_legend("stroke", title = "Disease(s) Selected")  %>%
      bind_shiny("plot2", "plot2_ui")
    
})
