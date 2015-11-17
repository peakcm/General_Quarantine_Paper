#### Header ####
# Figure with Rs as X axis and Rq as Y axis

#### Load libraries ####
library(ggplot2)
library(ggvis)
library(magrittr)
library(dplyr)
library(RColorBrewer)

#### Load Workspace ####
date_desired <- "20151104"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date_desired, "_FigureRsRq.RData", sep=""))

#### Load data from case_study_generic outputs ####
# Ebola
desired_root <- "20151024_Ebola"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.Ebola <- data.hr.lr
data.hr.lr.Ebola$disease <- "Ebola"

# SARS
desired_root <- "20151022_SARS"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.SARS <- data.hr.lr
data.hr.lr.SARS$disease <- "SARS"

# MERS
desired_root <- "20151027_MERS"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.MERS <- data.hr.lr
data.hr.lr.MERS$disease <- "MERS"

# Hepatitis A
desired_root <- "20151026_HepatitisA"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.HepatitisA <- data.hr.lr
data.hr.lr.HepatitisA$disease <- "HepatitisA"

# Influenza A
desired_root <- "20151028_InfluenzaA"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.InfluenzaA <- data.hr.lr
data.hr.lr.InfluenzaA$disease <- "InfluenzaA"

# Pertussis
desired_root <- "20151026_Pertussis"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.Pertussis <- data.hr.lr
data.hr.lr.Pertussis$disease <- "Pertussis"

# Smallpox
desired_root <- "20151028_Smallpox"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

data.hr.lr.Smallpox <- data.hr.lr
data.hr.lr.Smallpox$disease <- "Smallpox"

#### Set ranges for each disease for Plot2 ####
set_1 <- c("Ebola", 1.72, 1.94, "HR")
names(set_1) <- c("name", "lower", "upper", "Setting")

set_2 <- c("SARS", 2.2, 3.6, "HR")
names(set_2) <- c("name", "lower", "upper", "Setting")

set_3 <- c("MERS", 0.6, 1.3, "HR")
names(set_3) <- c("name", "lower", "upper", "Setting")

set_4 <- c("HepatitisA", 2, 2.5, "HR")
names(set_4) <- c("name", "lower", "upper", "Setting")

set_5 <- c("InfluenzaA", 1.28, 1.8, "HR")
names(set_5) <- c("name", "lower", "upper", "Setting")

set_6 <- c("Pertussis", 4.5, 5, "HR")
names(set_6) <- c("name", "lower", "upper", "Setting")

set_7 <- c("Smallpox", 4.5, 5, "HR")
names(set_7) <- c("name", "lower", "upper", "Setting")

#### Combine datasets ####
data_master <- rbind(data.hr.lr.Ebola,
                     data.hr.lr.Smallpox,
                     data.hr.lr.SARS,
                     data.hr.lr.MERS,
                     data.hr.lr.HepatitisA,
                     data.hr.lr.InfluenzaA,
                     data.hr.lr.Pertussis)

#### Plot settings ####
xlim.min <- 0
xlim.max <- 3
ylim.min <- 0
ylim.max <- 3

scale_colour_brewer(type="qual", palette=6)
my.cols <- brewer.pal(n = 7, name = "Set1")
my.cols <- my.cols[c(3, 7, 4, 5, 1, 6, 2)]

#### Plot 1: All same R ####
set_all <- c("All", 2.5, 3, "HR")
names(set_all) <- c("name", "lower", "upper", "Setting")

plot1 <- # data
  ggplot(data_master, aes(color=disease)) + 
  geom_point(data = data_master[data_master$disease == set_1["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_1["Setting"],],
               aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_2["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_2["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_3["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_3["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_4["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_4["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_5["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_5["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_6["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_6["Setting"],],
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_7["name"] &
                                  data_master$R_0 >= set_all["lower"] & 
                                  data_master$R_0 <= set_all["upper"] &
                                  data_master$Setting == set_7["Setting"],],
             aes(x=R_s, y=R_q)) +
  scale_color_manual(values = my.cols)

#### Plot 2: Disease-specific R ####
plot2 <- # data
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
             aes(x=R_s, y=R_q)) +
  geom_point(data = data_master[data_master$disease == set_7["name"] &
                                  data_master$R_0 >= set_7["lower"] & 
                                  data_master$R_0 <= set_7["upper"] &
                                  data_master$Setting == set_7["Setting"],],
             aes(x=R_s, y=R_q)) +
  scale_color_manual(values = my.cols)

#### Plot with accoutrements ####
plot1_pretty <- plot1 +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(xlim.min, xlim.max) + ylim(ylim.min, ylim.max) + 
  xlab("Effective Reproductive Number under Symptom Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  geom_vline(x=1, col="grey") + geom_hline(y=1, col="grey") +
  annotate("rect", xmin = xlim.min, xmax = 1, ymin = 1, ymax = ylim.max, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = xlim.max, ymin = ylim.min, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = xlim.min, xmax = 1, ymin = ylim.min, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = xlim.max - 1, y = ylim.min + 0.1, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = xlim.min + 0.5, y = 1.5, label = "Control with\nSymptom Monitoring", col = "orange")
plot1_pretty

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PlotRsRq1.pdf", sep=""))
plot(plot1_pretty)
dev.off()

plot2_pretty <- plot2 +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(xlim.min, xlim.max) + ylim(ylim.min, ylim.max) + 
  xlab("Effective Reproductive Number under Symptom Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  geom_vline(x=1, col="grey") + geom_hline(y=1, col="grey") +
  annotate("rect", xmin = xlim.min, xmax = 1, ymin = 1, ymax = ylim.max, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = xlim.max, ymin = ylim.min, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = xlim.min, xmax = 1, ymin = ylim.min, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = xlim.max - 1, y = ylim.min + 0.1, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = xlim.min + 0.5, y = 1.5, label = "Control with\nSymptom Monitoring", col = "orange")
plot2_pretty

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PlotRsRq2.pdf", sep=""))
plot(plot2_pretty)
dev.off()

#### ggvis: Interactive Supplement for this figure ####
diseases <- unique(data_master$disease)
slider <- input_slider(min=1, max=5, step=0.01, label = "Reproductive Number (+/- 0.5)")
data_master %>% 
  ggvis(x = ~R_s, y = ~R_q, fill = ~disease, stroke = ~disease, opacity := 0.5) %>%
  filter(R_0 <= (eval(slider)+0.5)) %>%
  filter(R_0 >= (eval(slider)-0.5)) %>%
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  # layer_rects(x = 0, x2 = 1, y = 1, y2 = 4, opacity := 0.1, fill = "yellow") %>%
#   layer_rects(x = 1, x2 = 4, y = 0, y2 = 1, opacity := 0.1, fill = "blue") %>%
#   layer_rects(x = 0, x2 = 1, y = 0, y2 = 1, opacity := 0.1, fill = "green") %>%
  layer_points(stroke := NA) %>%
  # layer_points(fill := NA) %>%
  layer_points(x = eval(slider), y = eval(slider), fill := "grey", shape := "cross") %>%
#   layer_text(x = eval(slider), y = eval(slider), text := "..R") %>%
  scale_numeric("x", domain = c(0, 5), nice = FALSE, label = "Symptom Monitoring") %>%
  scale_numeric("y", domain = c(0, 5), nice = FALSE, label = "Quarantine")
  

#### ggvis: model outputs across levels of R_0 for all diseases ####
# Several outcomes to select
slider2 <- input_slider(min=0, max=5, c(0,5),step = 0.1, label = "Y axis upper limit")
data_master %>%
  ggvis(~R_0,  input_select(c("R_s", "R_q", "Abs_Benefit","Rel_Benefit", "Abs_Benefit_per_Qday", "Rel_Benefit_per_Qday"), map=as.name, label = "Outcome"), fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", domain = slider2, clamp=TRUE, label = "Outcome")

#### ggvis: model outputs across levels of R_0, including lowess curves ####
# Rel_Benefit
data_master %>%
  ggvis(x = ~R_0, y = ~Rel_Benefit*100, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", domain = c(-50, 100), nice = FALSE, label = "% reduction in R by Q over SM")

# Abs_Benefit
data_master %>%
  ggvis(x = ~R_0, y = ~Abs_Benefit, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", domain = c(-1, 5), nice = FALSE, label = "Reduction in R by Q over SM")

# NNQ
data_master %>%
  ggvis(x = ~R_0, y = ~log10(NNQ), fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  filter(NNQ <= eval(input_slider(value=100, min = 50, max = 1000, step = 10, label = "Outlier Threshold"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", nice = FALSE, label = "Log 10 NNQ")

# Abs_Benefit_per_Qday
data_master %>%
  ggvis(x = ~R_0, y = ~Abs_Benefit_per_Qday, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", nice = FALSE, label = "Reduction in R by Q over SM per day of quarantine")

# Rel_Benefit_per_Qday
data_master$Rel_Benefit_per_Qday <- data_master$Rel_Benefit / data_master$obs_to_iso_q
data_master %>%
  ggvis(x = ~R_0, y = ~Rel_Benefit_per_Qday*100, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", nice = FALSE, label = "% reduction in R by Q over SM per day of quarantine")

# obs_to_iso_q
data_master %>%
  ggvis(x = ~R_0, y = ~obs_to_iso_q, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0") %>%
  scale_numeric("y", domain = c(0, 30), nice = FALSE, label = "Number of days in quarantine")

# ggplot version instead of ggvis

# look at Rel_Benefit for one disease across a range of R_0
R_0_input <- 3
sd <- 0.5
for (R_0_input in c(1,2,3,4)){
  summary <- summary(data_master[data_master$disease == "Smallpox" &
                                   data_master$R_0 <= R_0_input + sd & 
                                   data_master$R_0 >= R_0_input - sd, "Rel_Benefit"])
  cat("\nR_0_input is", R_0_input, "\n")
  print(summary)
}

ggplot(data_master[data_master$Setting == "HR",], aes(x=R_0, y=Rel_Benefit, group=disease, col=disease)) +
  geom_point(alpha = 0.4) +
  stat_smooth(size = 1.5, method = "loess") +
  theme_bw()

#### ggvis: plots for effect of R_s and R_q as compared to R_0 ####

# Compare R_s to R_0
data_master$R0_Rs <- data_master$R_0 - data_master$R_s
data_master %>%
  ggvis(x = ~R_0, y = ~R0_Rs/R_0, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0")

# Compare R_q to R_0
data_master$R0_Rq <- data_master$R_0 - data_master$R_q
data_master %>%
  ggvis(x = ~R_0, y = ~R0_Rq/R_0, fill = ~disease, stroke = ~disease) %>%
  group_by(disease) %>% 
  filter(Setting == eval(input_radiobuttons(selected = "HR",label = "Setting", choices = c("HR","LR")))) %>%
  filter(disease %in% eval(input_checkboxgroup(diseases, select = "Ebola"))) %>%
  layer_points(opacity := 0.1) %>%
  layer_smooths() %>%
  scale_numeric("x", domain = c(1, 5), nice = FALSE, label = "R_0")
  

#### Save Workspace Image ####
date <- format(Sys.time(), "%Y%m%d")
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_FigureRsRq.RData", sep=""))
write.csv(data_master, paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_FigureRsRq.csv", sep=""))
