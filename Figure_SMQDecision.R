#### Header ####
# Code for SMQDecision plot

#### Load Libraries ####
library(ggplot2)
library(reshape2)

#### Load Workspace ####
desired_root <- "20151024_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_MR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_MR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMQDecision.RData", sep=""))

#### Create dataframe ####
R_0_relevant.min <- 1.7
R_0_relevant.max <- 2

data.hr$Setting <- "HR"
data.mr$Setting <- "MR"
data.lr$Setting <- "LR"
data.hr.mr.lr <- rbind(data.hr, data.mr, data.lr)

data.hr.mr.lr.melt <- melt(data.hr.mr.lr, id = c("Setting"))

Setting  = c("HR", "MR", "LR")
data.SMQ <- data.frame(Setting = Setting)
data.SMQ$R_s <- NA
data.SMQ$R_q <- NA
data.SMQ$R_s.var <- NA
data.SMQ$R_q.var<- NA

for (i in 1:nrow(data.SMQ)){
  q <- mean(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                       data.hr.mr.lr$R_0_input > R_0_relevant.min &
                       data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_q"])
  q.var <- var(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                            data.hr.mr.lr$R_0_input > R_0_relevant.min &
                            data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_q"])
  sm <- mean(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                            data.hr.mr.lr$R_0_input > R_0_relevant.min &
                            data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_s"])
  sm.var <- var(data.hr.mr.lr[data.hr.mr.lr$R_0_input < R_0_relevant.max &
                             data.hr.mr.lr$R_0_input > R_0_relevant.min &
                             data.hr.mr.lr$Setting == data.SMQ[i, "Setting"], "R_s"])
  data.SMQ[i, "R_s"] <- sm
  data.SMQ[i, "R_s.var"] <- sm.var
  data.SMQ[i, "R_q"] <- q
  data.SMQ[i, "R_q.var"] <- q.var
}
data.SMQ.melt <- melt(data.SMQ, id = c("Setting", "R_q.var", "R_s.var"))
data.SMQ.melt$variance <- NA
data.SMQ.melt[data.SMQ.melt$variable == "R_q", "variance"] <- data.SMQ.melt[data.SMQ.melt$variable == "R_q", "R_q.var"]
data.SMQ.melt[data.SMQ.melt$variable == "R_s", "variance"] <- data.SMQ.melt[data.SMQ.melt$variable == "R_s", "R_s.var"]
data.SMQ.melt <- data.SMQ.melt[,c("Setting", "variable", "value", "variance")]

data.SMQ.melt$variable <- factor(data.SMQ.melt$variable, levels = c("R_q", "R_s"), ordered = TRUE)

#### Plot 1 ####
plot1 <- ggplot(data.SMQ.melt, aes(Setting, y = value, fill = variable)) +
  theme_bw() +
  theme(legend.justification=c(0,0), legend.position=c(0.1,0.85)) +
  geom_hline(y=1, color = "lightgreen", size=2) +
  geom_bar(stat = "identity", position="dodge", width = 0.8, size = 1) +
  geom_errorbar(aes(ymin = value - 1.96*sqrt(variance), ymax = value + 1.96*sqrt(variance)), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values = c("cornflowerblue", "darkgoldenrod"), labels = c("Quarantine", "Symptom Monitoring")) +
  scale_x_discrete(limits = rev(c("LR", "MR", "HR")), labels = rev(c("Low Resource", "Mid Resource", "High Resource")), name = "") +
  ylab(expression(paste(R[e]))) +
  guides(fill = guide_legend(title = element_blank()))
plot1

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotSMQ1.pdf", sep=""), height = 7, width = 4)
plot(plot1)
dev.off()

#### Calculate RBQD and NQD ####
data.hr.mr.lr$NQD <- 1/data.hr.mr.lr$Abs_Benefit_per_Qday
data.hr.mr.lr$Rel_Benefit_per_Qday <- data.hr.mr.lr$Rel_Benefit / data.hr.mr.lr$obs_to_iso_q

summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "NQD"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "NQD"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "NQD"])

summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit_per_Qday"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit_per_Qday"])
summary(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit_per_Qday"])

#### Calculate NQD50 ####
# NQDX = [ (X*D1) + ((1-X)*D0) ] / (Rs – Rq)

D1_HR <- data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "obs_to_iso_q"]
D1_MR <- data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "obs_to_iso_q"]
D1_LR <- data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "obs_to_iso_q"]

D0 <- 21 # for ebola

x <- 0.5
data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "NQD50"] <- (D1_HR + (1/x - 1)*(D0)) / (data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Abs_Benefit"])
data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "NQD50"] <- (D1_MR + (1/x - 1)*(D0)) / (data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Abs_Benefit"])
data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "NQD50"] <- (D1_LR + (1/x - 1)*(D0)) / (data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Abs_Benefit"])

median(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "NQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "NQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "NQD50"])

#### Calculate RBQD50 ####
# Rel_Benefit_per_QdayX = [ (Rs – Rq)/Rs ] / [ (X*D1) + ((1-X)*D0) ]

data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "RBQD50"] <- data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit"] / (D1_HR + (1/x - 1)*(D0))
data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "RBQD50"] <- data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit"] / (D1_MR + (1/x - 1)*(D0))
data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "RBQD50"] <- data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit"] / (D1_LR + (1/x - 1)*(D0))

median(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "RBQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "RBQD50"])
median(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "RBQD50"])


#### Explore Data ####
layout(cbind(c(1, 2, 3), c(4,5,6)))
hist(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit"], xlim = c(0, max(data.hr.mr.lr$Rel_Benefit)), main = "Relative Benefit in HR setting")
hist(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit"], xlim = c(0, max(data.hr.mr.lr$Rel_Benefit)), main = "Relative Benefit in MR setting")
hist(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit"], xlim = c(0, max(data.hr.mr.lr$Rel_Benefit)), main = "Relative Benefit in LR setting")
hist(D1_HR, xlim = c(0, max(D1_HR)), main = "Days in Q for HR setting")
hist(D1_MR, xlim = c(0, max(D1_HR)), main = "Days in Q for MR setting")
hist(D1_LR, xlim = c(0, max(D1_HR)), main = "Days in Q for LR setting")

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_SMQDecision.RData", sep=""))

