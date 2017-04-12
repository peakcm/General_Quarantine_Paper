#### Header ####
# Create a figure that shows a heatmap for how gamma and Prob_CT influence the relative benefit (heat) of Q over SM

#### Load Libraries ####

#### Load Workspaces ####
desired_root <- "20151024_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.ebola <- data.prcc
data.prcc.ebola$disease <- "Ebola"

desired_root <- "20151022_SARS" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.sars <- data.prcc
data.prcc.sars$disease <- "SARS"

desired_root <- "20151026_HepatitisA" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.hepa <- data.prcc
data.prcc.hepa$disease <- "HepatitisA"

desired_root <- "20151026_Pertussis" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.pertussis <- data.prcc
data.prcc.pertussis$disease <- "Pertussis"

desired_root <- "20151027_MERS" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.mers <- data.prcc
data.prcc.mers$disease <- "MERS"

desired_root <- "20151028_InfluenzaA" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.flu <- data.prcc
data.prcc.flu$disease <- "InfluenzaA"

desired_root <- "20151028_Smallpox" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
data.prcc.smallpox <- data.prcc
data.prcc.smallpox$disease <- "Smallpox"

#### Generate Dataset ####
data.prcc_master <- rbind(data.prcc.ebola,
                          data.prcc.sars,
                          data.prcc.hepa,
                          data.prcc.pertussis,
                          data.prcc.mers,
                          data.prcc.flu,
                          data.prcc.smallpox)

#### Plot 1 ####
ggplot(data.prcc_master, aes(x=gamma, y=prob_CT)) +
  geom_point(data = data.prcc_master[data.prcc_master$Rel_Benefit > 0.0 & data.prcc_master$Rel_Benefit < 0.25,], aes(x=gamma, y=prob_CT), color = "lightgrey", alpha = 0.5, size = 4) +
  geom_point(data = data.prcc_master[data.prcc_master$Rel_Benefit > 0.25 & data.prcc_master$Rel_Benefit < 0.5,], aes(x=gamma, y=prob_CT), color = "red", alpha = 0.5, size = 4) +
  geom_point(data = data.prcc_master[data.prcc_master$Rel_Benefit > 0.5 & data.prcc_master$Rel_Benefit < 0.75,], aes(x=gamma, y=prob_CT), color = "blue", alpha = 0.5, size = 4) +
  geom_point(data = data.prcc_master[data.prcc_master$Rel_Benefit > 0.75,], aes(x=gamma, y=prob_CT), color = "green", alpha = 0.5, size = 4) +
  facet_grid(.~disease) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#### Plot 2 ####
plot2 <- ggplot(data.prcc.hepa, aes(x=gamma, y=prob_CT)) +
  geom_point(data = data.prcc.hepa[data.prcc.hepa$Rel_Benefit > 0.0 & data.prcc.hepa$Rel_Benefit < 0.25,], aes(x=gamma, y=prob_CT), color = "lightgrey", alpha = 0.5, size = 4) +
  geom_point(data = data.prcc.hepa[data.prcc.hepa$Rel_Benefit > 0.25 & data.prcc.hepa$Rel_Benefit < 0.5,], aes(x=gamma, y=prob_CT), color = "red", alpha = 0.5, size = 4) +
  geom_point(data = data.prcc.hepa[data.prcc.hepa$Rel_Benefit > 0.5 & data.prcc.hepa$Rel_Benefit < 0.75,], aes(x=gamma, y=prob_CT), color = "blue", alpha = 0.5, size = 4) +
  geom_point(data = data.prcc.hepa[data.prcc.hepa$Rel_Benefit > 0.75,], aes(x=gamma, y=prob_CT), color = "green", alpha = 0.5, size = 4) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Isolation Effectivenss") + ylab("Percent of Contacts Traced") +
  scale_x_continuous(limits = c(0.2, 1.1), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = c("0", "25%", "50%", "75", "100%")) +
  scale_y_continuous(limits = c(0.2, 1.1), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = c("0", "25%", "50%", "75", "100%")) +
  annotate("text", x = 1.05, y = 0.9, label = ">75%", col = "green", size = 5) +
  annotate("text", x = 1.05, y = 0.7, label = ">50%", col = "blue", size = 5) +
  annotate("text", x = 1.05, y = 0.35, label = ">25%", col = "red", size = 5) +
  annotate("text", x = 1.05, y = 1.05, label = "Hepatitis A\nRelative Benefit", col = "black", size = 3) +
  geom_rect(xmin = 0.19, xmax = 1.01, ymin = 0.19, ymax = 1.01, fill = NA, color = "darkgrey")
plot2

date <- format(Sys.time(), "%Y%m%d")
pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_Gamma_ProbCT_HepA.pdf", sep=""))
plot(plot2)
dev.off()
  
#### Plot 3 ####
ggplot(data.prcc_master, aes(x=gamma, y=prob_CT, shape = disease, color = Rel_Benefit, size = Rel_Benefit)) +
  geom_point() + 
  scale_color_continuous(low = "darkgrey", high = "lightgreen") +
  theme_bw()

