#### Header ####
# Create a figure that shows how Rel_Benefit and Rel_Benefit_per_Qday are different for each disease
# Grouped bar chart by each disease. One bar for each outcome

#### Load Libraries ####
library(ggplot2)

#### Load Workspace _FigureRsRq.RData ####
desired_date <- "20151113"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_date, "_FigureRsRq.RData", sep=""))

#### New dataset ####
data_bar <- data.frame(disease = unique(data_master$disease))
data_bar$Rel_Benefit <- NA
data_bar$Rel_Benefit_per_Qday <- NA

for (disease in data_bar$disease){
  data_bar[data_bar$disease == disease, "Rel_Benefit"] <- mean(data_master[data_master$disease == disease & data_master$Setting == "HR", "Rel_Benefit"])
  data_bar[data_bar$disease == disease, "Rel_Benefit_per_Qday"] <- mean(data_master[data_master$disease == disease & data_master$Setting == "HR", "Rel_Benefit_per_Qday"])
}

data_bar$disease <- factor(as.character(data_bar$disease), levels = rev(c("Pertussis", "Smallpox", "SARS", "HepatitisA", "InfluenzaA", "Ebola", "MERS")), labels = rev(c("Pertussis", "Smallpox", "SARS", "Hepatitis A", "Influenza A", "Ebola", "MERS")), ordered = TRUE)

data_bar_melt <- melt(data_bar, id="disease")

#### Plot Settings ####
scale_colour_brewer(type="qual", palette=6)
my.cols <- brewer.pal(n = 7, name = "Set1")
my.cols <- my.cols[c(3, 7, 4, 5, 1, 6, 2)]

#### Plot ####
plot1 <- ggplot(data_bar_melt, aes(disease, y = value, fill = disease, color = variable)) +
  geom_bar(stat = "identity", position="dodge", width = 0.8, size = 1) +
  annotate("text", x = 6, y = 0.6, label = "Bold Border denotes Percent Reduction\nper Day of Quarantine", col = "black", size = 4) +
  scale_fill_manual(values = my.cols) +
  scale_color_manual(values = c(NA, "black")) +
  # scale_x_discrete(breaks = c("Pertussis", "Smallpox", "SARS", "HepatitisA", "InfluenzaA", "Ebola", "MERS"), labels = c("Pertussis", "Smallpox", "SARS", "Hepatitis A", "Influenza A", "Ebola", "MERS")) +
  ylab(expression(paste(frac(R[S]-R[Q],R[S])))) +
  xlab("Disease") +
  coord_flip() +
  theme(text = element_text(size=10)) +
  theme_bw() +  guides(color=FALSE, fill = FALSE)
plot1

date <- format(Sys.time(), "%Y%m%d")
pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_RelBenefitBars.pdf", sep=""), width = 8, height = 4)
plot(plot1)
dev.off()

#### Save Workspace Image ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_RelBenefitBars.RData", sep=""))

