#### Header ####
# Code for the R_e vs R_0 plots

#### Load Libraries ####
library(ggplot2)
library(reshape)

#### Load HR and LR workspaces ####
desired_root <- "20151104_InfluenzaA" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

#### Create new dataset ####
# Set range for relevant R_0 values
R_0_relevant.min <- 1
R_0_relevant.max <- 5

data.hr$Setting <- "HR"
data.lr$Setting <- "LR"
data.hr.lr <- rbind(data.hr, data.lr)

#### Plot 1 ####
plot1 <- ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min & data.hr.lr$R_0 < R_0_relevant.max,]) +
  geom_vline(xintercept=1, col="grey") + geom_hline(yintercept=1, col="grey") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = 4, ymin = 0, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = 2.5, y = 0.9, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = 0.5, y = 3.5, label = "Control with\nSymptom Monitoring", col = "orange") +
  geom_point(aes(x=R_s, y=R_q, col = R_0, shape=Setting) ) +
  xlim(0,4) + ylim(0,4) +
  scale_colour_gradient(low="yellow", high="darkred") +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("High Resource", "Low Resource")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab("Effective Reproductive Number under Symptom Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste("Disease: ", disease))
plot1

#### Plot 2 ####
plot2 <- ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min & data.hr.lr$R_0 < R_0_relevant.max,]) +
  annotate("rect", xmin = R_0_relevant.min, xmax = R_0_relevant.max, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  geom_hline(yintercept=1, col = "grey") +
  geom_point(aes(x=R_0, y=R_s, shape = Setting), col = "darkgreen", alpha = 0.7) +
  geom_point(aes(x=R_0, y=R_q, shape = Setting), col = "blue", alpha = 0.7) +
  stat_smooth(aes(x=R_0, y=R_s, shape = Setting), method = "loess", color="darkgreen", size = 1.2) +
  stat_smooth(aes(x=R_0, y=R_q, shape = Setting), method = "loess", color="blue", size = 1.2) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("High Resource", "Low Resource")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab(expression("Basic Reproductive Number R" [0])) + ylab(expression("Effective Reproductive Number R" [e])) +
  ggtitle(paste("Disease: ", disease))
plot2

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot2.pdf", sep=""))
plot(plot2)
dev.off()

#### Plot 3 ####
point_size = 0.1
line_size = 1
plot3 <- ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min &
                                    data.hr.lr$R_0 < R_0_relevant.max &
                                    data.hr.lr$Setting == "HR",]) +
  annotate("rect", xmin = R_0_relevant.min, xmax = R_0_relevant.max, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  geom_point(aes(x=R_0, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = point_size) +
  geom_point(aes(x=R_0, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = point_size) +
  geom_point(aes(x=R_0, y=R_s), col = "darkgoldenrod", alpha = 0.5, size = point_size) +
  geom_point(aes(x=R_0, y=R_q), col = "cornflowerblue", alpha = 0.5, size = point_size) +
  stat_smooth(aes(x=R_0, y=R_hsb), method = "loess", color="mediumturquoise", size = line_size, se=FALSE, alpha = 0.5) +
  stat_smooth(aes(x=R_0, y=R_s), method = "loess", color="darkgoldenrod", size = line_size, se=FALSE, alpha = 0.5) +
  stat_smooth(aes(x=R_0, y=R_q), method = "loess", color="cornflowerblue", size = line_size, se=FALSE, alpha = 0.5) +
  geom_hline(yintercept=1, col = "grey", size = 0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression("R" [0])) + ylab(expression("R" [e])) +
  geom_text(data = NULL, x = 1.46, y = 0, label = "|", size = 2) +
  # geom_text(data = NULL, x = 1.83, y = 0, label = "|", size = 2) +
  # ggtitle(paste("Influenza A")) +
  # ggtitle(paste("Ebola")) +
  theme(text = element_text(size=8))
plot3

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot3.pdf", sep=""), width = 2, height = 2)
plot(plot3)
dev.off()


#### Plot 3 FOR POSTER ####
plot3_poster <- ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min &
                                    data.hr.lr$R_0 < R_0_relevant.max &
                                    data.hr.lr$Setting == "HR",]) +
  annotate("rect", xmin = R_0_relevant.min, xmax = R_0_relevant.max, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  geom_hline(yintercept=1, col = "grey") +
  geom_point(aes(x=R_0, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = 0.5) +
  geom_point(aes(x=R_0, y=R_hsb), col = "mediumturquoise", alpha = 0.5, size = 0.5) +
  geom_point(aes(x=R_0, y=R_s), col = "darkgoldenrod", alpha = 0.5, size = 0.5) +
  geom_point(aes(x=R_0, y=R_q), col = "cornflowerblue", alpha = 0.5, size = 0.5) +
  stat_smooth(aes(x=R_0, y=R_hsb), method = "loess", color="mediumturquoise", size = 1, se=FALSE) +
  stat_smooth(aes(x=R_0, y=R_s), method = "loess", color="darkgoldenrod", size = 1, se=FALSE) +
  stat_smooth(aes(x=R_0, y=R_q), method = "loess", color="cornflowerblue", size = 1, se=FALSE) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression("Basic Reproductive Number R" [0])) + ylab(expression("Effective Reproductive Number R" [e])) +
  geom_text(data = NULL, x = 1.46, y = 0, label = "*") +
  # geom_text(data = NULL, x = 1.83, y = 0, label = "*") +
  theme(text = element_text(size=18)) +
  ggtitle(paste("Influenza A"))
# ggtitle(paste("Ebola"))
plot3_poster

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot3_5by5.pdf", sep=""), width = 5, height = 5)
plot(plot3_poster)
dev.off()

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plots.RData", sep=""))
