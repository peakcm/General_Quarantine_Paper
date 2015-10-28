#### Header ####
# Figure with Rs as X axis and Rq as Y axis

#### Load libraries ####
library(ggplot2)
library(ggvis)
library(magrittr)
library(dplyr)

#### Load data from case_study_generic outputs ####

# Generic
desired_root <- "20151026_Pertussis" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

# Ebola
desired_root <- "20151024_Ebola"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.Ebola <- data.hr.lr
data.hr.lr.Ebola$disease <- "Ebola"

set_1 <- c("Ebola", 2, 2.5, "HR")
names(set_1) <- c("name", "lower", "upper", "Setting")

# SARS
desired_root <- "20151022_SARS"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.SARS <- data.hr.lr
data.hr.lr.SARS$disease <- "SARS"

set_2 <- c("SARS", 2, 2.5, "HR")
names(set_2) <- c("name", "lower", "upper", "Setting")

# MERS
desired_root <- "20151025_MERS"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.MERS <- data.hr.lr
data.hr.lr.MERS$disease <- "MERS"

set_3 <- c("MERS", 2, 2.5, "HR")
names(set_3) <- c("name", "lower", "upper", "Setting")

# Hepatitis A
desired_root <- "20151026_HepatitisA"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.HepatitisA <- data.hr.lr
data.hr.lr.HepatitisA$disease <- "HepatitisA"

set_4 <- c("HepatitisA", 2, 2.5, "HR")
names(set_4) <- c("name", "lower", "upper", "Setting")

# Influenza A
desired_root <- "20151027_InfluenzaA"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.InfluenzaA <- data.hr.lr
data.hr.lr.InfluenzaA$disease <- "InfluenzaA"

set_5 <- c("InfluenzaA", 2, 2.5, "HR")
names(set_5) <- c("name", "lower", "upper", "Setting")

# Pertussis
desired_root <- "20151026_Pertussis"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.Pertussis <- data.hr.lr
data.hr.lr.Pertussis$disease <- "Pertussis"

set_6 <- c("Pertussis", 2, 2.5, "HR")
names(set_6) <- c("name", "lower", "upper", "Setting")

# Smallpox
desired_root <- "20151027_Smallpox"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.Smallpox <- data.hr.lr
data.hr.lr.Smallpox$disease <- "Smallpox"

set_6 <- c("Smallpox", 2, 2.5, "HR")
names(set_6) <- c("name", "lower", "upper", "Setting")

#### combine datasets ####
data_master <- rbind(data.hr.lr.Ebola,
                     data.hr.lr.SARS,
                     data.hr.lr.MERS,
                     data.hr.lr.HepatitisA,
                     data.hr.lr.InfluenzaA,
                     data.hr.lr.Pertussis,
                     data.hr.lr.Smallpox)

#### ggplot ####

plot <- # data
  ggplot(data_master, aes(color=disease)) +
  geom_point(data = data_master[data_master$disease == set_1["name"] &
                                  data_master$R_0 >= set_1["lower"] & 
                                  data_master$R_0 <= set_1["upper"] &
                                  data_master$Setting == set_1["Setting"],],
               aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_2["name"] &
                                  data_master$R_0 >= set_2["lower"] & 
                                  data_master$R_0 <= set_2["upper"] &
                                  data_master$Setting == set_2["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_3["name"] &
                                  data_master$R_0 >= set_3["lower"] & 
                                  data_master$R_0 <= set_3["upper"] &
                                  data_master$Setting == set_3["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_4["name"] &
                                  data_master$R_0 >= set_4["lower"] & 
                                  data_master$R_0 <= set_4["upper"] &
                                  data_master$Setting == set_4["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_5["name"] &
                                  data_master$R_0 >= set_5["lower"] & 
                                  data_master$R_0 <= set_5["upper"] &
                                  data_master$Setting == set_5["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_6["name"] &
                                  data_master$R_0 >= set_6["lower"] & 
                                  data_master$R_0 <= set_6["upper"] &
                                  data_master$Setting == set_6["Setting"],],
             aes(x=R_s, y=R_q))

plot + # accoutrements
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(0, 4) + ylim(0, 4) + 
  xlab("Effective Reproductive Number under Symptom Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  geom_vline(x=1, col="grey") + geom_hline(y=1, col="grey") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = 4, ymin = 0, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = 3, y = 0.1, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = 0.5, y = 1.5, label = "Control with\nSymptom Monitoring", col = "orange")
  
#### ggvis ####
p <- ggvis(data = data_master, x = ~R_s, y = ~R_q) 
layer_points(p)

v_line <- data.frame(cbind(c(1.0000, 1.0001), c(0, 5)))
names(v_line) <- c("R_s", "R_q")

diseases <- unique(data_master$disease)
slider <- input_slider(min=1, max=5, step=0.01, label = "Reproductive Number")
data_master %>% 
  ggvis(x = ~R_s, y = ~R_q, fill = ~disease, opacity := 0.5) %>%
  filter(R_0 <= (eval(slider)+0.5)) %>%
  filter(R_0 >= (eval(slider)-0.5)) %>%
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
#   layer_rects(x = 0, x2 = 1, y = 1, y2 = 4, opacity := 0.1, fill = "yellow") %>%
#   layer_rects(x = 1, x2 = 4, y = 0, y2 = 1, opacity := 0.1, fill = "blue") %>%
#   layer_rects(x = 0, x2 = 1, y = 0, y2 = 1, opacity := 0.1, fill = "green") %>%
  layer_points() %>%
  scale_numeric("x", domain = c(0, 5), nice = FALSE, label = "Symptom Monitoring") %>%
  scale_numeric("y", domain = c(0, 5), nice = FALSE, label = "Quarantine")

  

