#### Header ####
# Create a figure that has Re on the x axis and then each column is a disease. There is a triangle shape for the Rsm and a square for Rq and a circle for R0
# Connect the two with a line. The length of the line is Rsm-Rq.
# Set the color of the points and line to match the other figures
# Draw a vertical line at Re = 1 to symbolize control
# Draw a dashed horizontal line to group off the disease where 1 < Rq < Rs from Rq < 1 < Rs and also from Rq < Rs < 1

#### Load Libraries ####
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

#### Set colors ####
scale_colour_brewer(type="qual", palette=6)
my.cols <- brewer.pal(n = 7, name = "Set1")
my.cols <- my.cols[c(3, 7, 4, 5, 1, 6, 2)]

#### Load data ####
load("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/20151113_FigureRsRq.RData")

#### Set ranges for R that are the same for each disease ####
lower <- 2.5
upper <- 3

data_master_subset <- data_master[data_master$R_0 > lower & data_master$R_0 < upper & data_master$Setting == "HR",]

#### Set ranges for R for each disease ####
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

data_master_subset.1 <- data_master[data_master$disease == set_1["name"] & data_master$R_0 > set_1["lower"] & data_master$R_0 < set_1["upper"] & data_master$Setting == "HR",]
data_master_subset.2 <- data_master[data_master$disease == set_2["name"] & data_master$R_0 > set_2["lower"] & data_master$R_0 < set_2["upper"] & data_master$Setting == "HR",]
data_master_subset.3 <- data_master[data_master$disease == set_3["name"] & data_master$R_0 > set_3["lower"] & data_master$R_0 < set_3["upper"] & data_master$Setting == "HR",]
data_master_subset.4 <- data_master[data_master$disease == set_4["name"] & data_master$R_0 > set_4["lower"] & data_master$R_0 < set_4["upper"] & data_master$Setting == "HR",]
data_master_subset.5 <- data_master[data_master$disease == set_5["name"] & data_master$R_0 > set_5["lower"] & data_master$R_0 < set_5["upper"] & data_master$Setting == "HR",]
data_master_subset.6 <- data_master[data_master$disease == set_6["name"] & data_master$R_0 > set_6["lower"] & data_master$R_0 < set_6["upper"] & data_master$Setting == "HR",]
data_master_subset.7 <- data_master[data_master$disease == set_7["name"] & data_master$R_0 > set_7["lower"] & data_master$R_0 < set_7["upper"] & data_master$Setting == "HR",]

data_master_subset <- rbind(data_master_subset.1, data_master_subset.2, data_master_subset.3, data_master_subset.4, data_master_subset.5, data_master_subset.6, data_master_subset.7)

#### Create a new data frame ####
data_master_melt <- melt(data_master_subset, id.vars = c("Abs_Benefit", "Rel_Benefit", "NNQ", "obs_to_iso_q", "Abs_Benefit_per_Qday", "ks", "d_inf",
                                                         "pi_t_triangle_center", "T_lat_offset", "dispersion", "Setting", "disease", "Rel_Benefit_per_Qday", "R0_Rs", "R0_Rq"))
head(data_master_melt)
unique(data_master_melt$variable)

grouped <- group_by(data_master_melt, disease, Setting, variable)
data_master_melt_mean <- summarise(grouped, mean=mean(value), sd=sd(value))

# reorder diseases
data_master_melt_mean$disease <- factor(data_master_melt_mean$disease, levels = rev(c("Pertussis", "Smallpox", "SARS", "HepatitisA", "InfluenzaA", "Ebola", "MERS")), labels = rev(c("Pertussis", "Smallpox", "SARS", "Hepatitis A", "Influenza A", "Ebola", "MERS")), ordered = TRUE)

my.cols <- brewer.pal(n = 7, name = "Set1")
my.cols <- my.cols[c(5, 3, 4, 7, 6, 2, 1)]

# Create a data frame for drawing a line segment
data_master_melt_mean_segment <- melt(data_master_melt_mean[data_master_melt_mean$variable %in% c("R_s", "R_q"),], id.vars = c("disease"), measure.vars = c("mean"))

#### Plot ####
plot <- ggplot(data_master_melt_mean[data_master_melt_mean$Setting == "HR" & data_master_melt_mean$variable %in% c("R_0","R_s", "R_q"),]) +
  theme_bw() + 
  geom_vline(xintercept = 1, color = "darkgrey") +

  scale_shape_manual(values = c(23, 24, 21),
                     name="Intervention",
                     breaks=c("R_0", "R_s", "R_q"),
                     labels=c(  expression(paste("None (",R[0], ")")), expression(paste("Symptom Monitoring (", R[S], ")")),expression(paste("Quarantine (", R[Q], ")")))) +
  theme(legend.text.align = 0)+
  # theme(legend.position=c(.8, .25)) +
  # theme(legend.background = element_rect(color = "white")) +
  # theme(legend.text = element_text(size = 10, face = "bold")) +
  # theme(legend.title = element_blank()) +
  scale_color_manual(values = my.cols, guide = FALSE) +
  scale_fill_manual(values = my.cols, guide = FALSE) +
  theme(axis.title.y=element_blank()) +
  geom_line(data = data_master_melt_mean_segment, aes(y = disease, x = value, color = disease, fill = disease), lwd = 1) +
  annotate("rect", xmin = 0, xmax = 5, ymin = 0.5, ymax = 2.5, alpha = .08, fill = "green") +
  annotate("rect", xmin = 0, xmax = 5, ymin = 2.5, ymax = 6.5, alpha = .08, fill = "blue") +
  geom_point(aes(x = mean, y = disease, group = disease, fill = disease, shape = variable), size = 3, color = "darkgrey") +
  xlim(0,5) + 
  xlab(expression(R[0])) +
  theme(text = element_text(size=10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank())
plot

date <- format(Sys.time(), "%Y%m%d")

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PlotRsRq_horizontal_8by2-5.pdf", sep=""), width = 8, height = 2.5)
plot(plot)
dev.off()

  

#### Plot FOR POSTER ####
plot_poster <- ggplot(data_master_melt_mean[data_master_melt_mean$Setting == "HR" & data_master_melt_mean$variable %in% c("R_0","R_s", "R_q"),]) +
  theme_bw() + 
  geom_vline(xintercept = 1, color = "darkgrey") +
  
  scale_shape_manual(values = c(23, 24, 21),
                     name="Intervention",
                     breaks=c("R_0", "R_s", "R_q"),
                     labels=c(  expression(paste("None (",R[0], ")")), expression(paste("Symptom Monitoring (", R[S], ")")),expression(paste("Quarantine (", R[Q], ")")))) +
  theme(legend.text.align = 0)+
  # theme(legend.position=c(.8, .25)) +
  # theme(legend.background = element_rect(color = "white")) +
  # theme(legend.text = element_text(size = 10, face = "bold")) +
  # theme(legend.title = element_blank()) +
  scale_color_manual(values = my.cols, guide = FALSE) +
  scale_fill_manual(values = my.cols, guide = FALSE) +
  theme(axis.title.y=element_blank()) +
  geom_line(data = data_master_melt_mean_segment, aes(y = disease, x = value, color = disease, fill = disease), lwd = 1) +
  annotate("rect", xmin = 0, xmax = 5, ymin = 0.5, ymax = 2.5, alpha = .08, fill = "green") +
  annotate("rect", xmin = 0, xmax = 5, ymin = 2.5, ymax = 6.5, alpha = .08, fill = "blue") +
  geom_point(aes(x = mean, y = disease, group = disease, fill = disease, shape = variable), size = 3, color = "darkgrey") +
  xlim(0,5) + 
  xlab(expression(R[0])) +
  theme(text = element_text(size=18)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank())
plot_poster

date <- format(Sys.time(), "%Y%m%d")

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PlotRsRq_horizontal_8by2-5.pdf", sep=""), width = 10, height = 2.5)
plot(plot_poster)
dev.off()


