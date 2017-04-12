#### Header ####
# Figure to plot the Rel_Benefit and the Rel_Benefit_per_Qday (Y) with respect the intervention metrics like gamma and risk profiling

#### Load Libraries ####
library(ggplot2)

#### Load Workspaces ####
desired_root <- "20151024_Ebola"

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMQDecision.RData", sep=""))

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_IntEffect.RData", sep=""))

gbr_cols <- c("#00BA38", "#619CFF", "#F8766D")

#### Vary proportion of contacts who are infected (X) ####
data.hr.mr.lr$Rel_Benefit_positive <- data.hr.mr.lr$Rel_Benefit
data.hr.mr.lr[data.hr.mr.lr$Rel_Benefit < 0, "Rel_Benefit_positive"] <- 0

sum(data.hr.mr.lr$Rel_Benefit_per_Qday != (data.hr.mr.lr$Rel_Benefit / data.hr.mr.lr$obs_to_iso_q))

rand_rp <- 1-runif(n=nrow(data.hr.mr.lr)/3, min = 0, max = 1)^3

data.hr.mr.lr$riskprofile <- rep(rand_rp, times = 3)

maxInc <- 21
data.hr.mr.lr$obs_to_iso_q_rp <- (data.hr.mr.lr$obs_to_iso_q + (1/data.hr.mr.lr$riskprofile - 1)*maxInc)
data.hr.mr.lr$RBQDrp <- data.hr.mr.lr$Rel_Benefit / (data.hr.mr.lr$obs_to_iso_q + (1/data.hr.mr.lr$riskprofile - 1)*maxInc)
data.hr.mr.lr$RBQDrp_positive <- data.hr.mr.lr$Rel_Benefit_positive / (data.hr.mr.lr$obs_to_iso_q + (1/data.hr.mr.lr$riskprofile - 1)*maxInc)

data.hr.mr.lr$Setting <- factor(data.hr.mr.lr$Setting, levels = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), ordered = TRUE)

# Rel_Benefit
ggplot(data.hr.mr.lr, aes(x=riskprofile, y=Rel_Benefit*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Fraction of Contacts Truly Infected (", P[inf], ")", sep=""))) +
  ylab(expression(paste("Relative Difference ", frac(R[S]-R[Q],R[S]), sep=""))) +
  # scale_y_continuous(breaks = seq(0, 20, by=5), labels = c("0%", "5%", "10%", "15%", "20%")) +
  # scale_y_log10(breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1), labels = c("0.005%", "0.01%", "0.05%", "0.1%", "0.5%", "1%")) +
  # scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr, aes(x=riskprofile, y=Rel_Benefit_positive*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") + 
  xlab(expression(paste("Fraction of Contacts Truly Infected (", P[inf], ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence ", frac(R[S]-R[Q],R[S]), sep=""))) +
  # scale_y_continuous(breaks = seq(0, 20, by=5), labels = c("0%", "5%", "10%", "15%", "20%")) +
  # scale_y_log10(breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1), labels = c("0.005%", "0.01%", "0.05%", "0.1%", "0.5%", "1%")) +
  # scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) +
  geom_point(alpha = 0.2)

# Rel_Benefit_per_Qday
ggplot(data.hr.mr.lr, aes(x=riskprofile, y=RBQDrp*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +  
  xlab(expression(paste("Fraction of Contacts Truly Infected (", P[inf], ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence per Quarantine Day  ", (frac(R[S]-R[Q],R[S]))/d[Q], sep=""))) +
  # scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by=1), labels = c("0%", "1%", "2%", "3%", "4%", "5%", "6%", "7%", "8%", "9%", "10%")) +
  # scale_y_log10(breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1), labels = c("0.005%", "0.01%", "0.05%", "0.1%", "0.5%", "1%")) +
  # scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr, aes(x=riskprofile, y=RBQDrp_positive*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Fraction of Contacts Truly Infected (", P[inf], ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence per Quarantine Day  ", (frac(R[S]-R[Q],R[S]))/d[Q], sep=""))) +
  scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by=1), labels = c("0%", "1%", "2%", "3%", "4%", "5%", "6%", "7%", "8%", "9%", "10%")) +
  # scale_y_log10(breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1), labels = c("0.005%", "0.01%", "0.05%", "0.1%", "0.5%", "1%")) +
  # scale_y_log10() +
  # scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  geom_point(alpha = .2) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) 

#### Vary proportion of contacts who are infected (X): Analytic version ####
x.vals <- seq(0.01, 1, 0.01)
test <- data.frame(x <- x.vals)
test$NQD <- NA
test$RBQD <- NA
test$Setting <- "HR"
df <- test
test$Setting <- "MR"
df <- rbind(df, test)
test$Setting <- "LR"
df <- rbind(df, test)
names(df) <- c("x.vals", "NQD", "RBQD", "Setting")
df$Setting <- factor(df$Setting, levels = c("HR", "MR", "LR"))
df$RB <- NA

for (i in 1:nrow(df)){
  x <- df[i, "x.vals"]
  Setting <- df[i, "Setting"]
  
  if (Setting == "HR"){
    NQD <- median( ( (D1_HR + (1/x - 1)*(D0))) / data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Abs_Benefit"])
    RBQD <- median(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit"] / ( D1_HR + (1/x - 1)*(D0)))
    RB <- median(data.hr.mr.lr[data.hr.mr.lr$Setting == "HR", "Rel_Benefit"])
  }
  if (Setting == "MR"){
    NQD <- median( ( (D1_HR + (1/x - 1)*(D0))) / data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Abs_Benefit"])
    RBQD <- median(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit"] / ( (D1_HR + (1/x - 1)*(D0))))
    RB <- median(data.hr.mr.lr[data.hr.mr.lr$Setting == "MR", "Rel_Benefit"])
  }
  if (Setting == "LR"){
    NQD <- median(( (D1_HR + (1/x - 1)*(D0))) / data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Abs_Benefit"])
    RBQD <- median(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit"] / ( (D1_HR + (1/x - 1)*(D0))))
    RB <- median(data.hr.mr.lr[data.hr.mr.lr$Setting == "LR", "Rel_Benefit"])
  }
  
  df[i, "NQD"] <- NQD
  df[i, "RBQD"] <- RBQD
  df[i, "RB"] <- RB
}

plot2 <- ggplot(df, aes(x = x.vals, y = RBQD*100, color = Setting)) +
  theme_bw() +
  geom_line(size = 1.5) +
  theme(legend.justification=c(0,0), legend.position=c(0.8,0.1)) +
  # scale_y_log10(breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1), labels = c("0.005%", "0.01%", "0.05%", "0.1%", "0.5%", "1%")) +
  # scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  ggtitle("Relative Benefit of Quarantine Over\nSymptom Monitoring per Quarantine Day") +
  xlab("Proportion of contacts who are truly infected") + ylab("Percent Reduction per Quarantine Day") +
  scale_color_manual(c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), values = c("green", "blue", "red"), name = "Setting") +
  theme(panel.grid.minor = element_blank())
plot2

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotSMQ2.pdf", sep=""), height = 7, width = 7)
plot(plot2)
dev.off()

plot3 <- ggplot(df, aes(x = x.vals, y = NQD, color = Setting)) +
  theme_bw() +
  geom_line(size = 1.5) +
  theme(legend.justification=c(0,0), legend.position=c(0.8,0.7)) +
  scale_y_log10(breaks = c(50, 100, 500, 1000, 5000, 10000), labels = c("50", "100", "500", "1,000", "5,000", "10,000")) +
  scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  ggtitle("Number of Quarantine Days to Avert\nOne Case over Symptom Monitoring") +
  xlab("Proportion of contacts who are truly infected") + ylab("Days") +
  scale_color_manual(c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), values = c("green", "blue", "red"), name = "Setting") +
  theme(panel.grid.minor = element_blank())
plot3

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotSMQ3.pdf", sep=""), height = 7, width = 7)
plot(plot3)
dev.off()

plot4 <- ggplot(df, aes(x = x.vals, y = RB*100, color = Setting)) +
  theme_bw() +
  geom_line(size = 1.5) +
  theme(legend.justification=c(0,0), legend.position=c(0.8,0.1)) +
  # scale_y_log10(breaks = c(0.005, 0.01, 0.05, 0.1, 0.5, 1), labels = c("0.005%", "0.01%", "0.05%", "0.1%", "0.5%", "1%")) +
  # scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
  ylim(0, 20) +
  ggtitle("Relative Benefit of Quarantine Over\nSymptom Monitoring") +
  xlab("Proportion of contacts who are truly infected") + ylab("Percent Reduction in Re") +
  scale_color_manual(c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), values = c("green", "blue", "red"), name = "Setting") +
  theme(panel.grid.minor = element_blank())
plot4

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotSMQ4.pdf", sep=""), height = 7, width = 7)
plot(plot4)
dev.off()

#### Effect and efficiency of Q over R_0 is highest for HR then MR then LR ####
data.SMQ$Qday <- c(median(D1_HR), median(D1_MR), median(D1_LR))
data.SMQ$R_0 <- 1.8
data.SMQ$Q_effect <- data.SMQ$R_0 - data.SMQ$R_q
data.SMQ$Q_efficiency <- data.SMQ$Q_effect / data.SMQ$Qday

data.SMQ$S_effect <- data.SMQ$R_0 - data.SMQ$R_s
data.SMQ$S_efficiency <- data.SMQ$S_effect / data.SMQ$Qday

data.SMQ

#### Run model with varying Gamma in HR, MR, and LR ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))

# High Resource everything else
# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 100
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.hr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  gamma = runif(n = times, min = 0, max = 1),
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1.72, max = 1.94)) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  gamma <- params.set[i, "gamma"]
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == "s"){
      data.hr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.hr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.hr[,"Abs_Benefit"] <- data.hr[,"R_s"] - data.hr[,"R_q"]
data.hr[,"Rel_Benefit"] <- data.hr[,"Abs_Benefit"] / data.hr[,"R_s"]
data.hr[,"NNQ"] <- 1 / data.hr[,"Abs_Benefit"]
data.hr[data.hr$NNQ < 1,"NNQ"] <- 1
data.hr[data.hr$NNQ > 9999,"NNQ"] <- 9999
data.hr[data.hr$NNQ == Inf,"NNQ"] <- 9999
data.hr[,"Abs_Benefit_per_Qday"] <- data.hr[,"Abs_Benefit"] / data.hr[,"obs_to_iso_q"]
data.hr$d_inf <- params.set[,"d_inf"]
data.hr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.hr$R_0_input <- params.set[,"R_0"]
data.hr$T_lat_offset <- params.set[,"T_lat_offset"]
data.hr$dispersion <- params.set[,"dispersion"]
data.hr$gamma <- params.set[,"gamma"]

## Middle Resource ##
# Interventions
background_intervention <- "u"

prob_CT <- 0.75

gamma <- 0.75

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.mr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.mr) <- names

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  gamma <- params.set[i,"gamma"]
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == "s"){
      data.mr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.mr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.mr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.mr[,"Abs_Benefit"] <- data.mr[,"R_s"] - data.mr[,"R_q"]
data.mr[,"Rel_Benefit"] <- data.mr[,"Abs_Benefit"] / data.mr[,"R_s"]
data.mr[,"NNQ"] <- 1 / data.mr[,"Abs_Benefit"]
data.mr[data.mr$NNQ < 1,"NNQ"] <- 1
data.mr[data.mr$NNQ > 9999,"NNQ"] <- 9999
data.mr[data.mr$NNQ == Inf,"NNQ"] <- 9999
data.mr[,"Abs_Benefit_per_Qday"] <- data.mr[,"Abs_Benefit"] / data.mr[,"obs_to_iso_q"]
data.mr$d_inf <- params.set[,"d_inf"]
data.mr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.mr$R_0_input <- params.set[,"R_0"]
data.mr$T_lat_offset <- params.set[,"T_lat_offset"]
data.mr$dispersion <- params.set[,"dispersion"]
data.mr$gamma <- params.set[,"gamma"]

## Low Resource ##

# Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

parms_epsilon = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.lr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.lr) <- names

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  gamma <- params.set[i,"gamma"]
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == "s"){
      data.lr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.lr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.lr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.lr[,"Abs_Benefit"] <- data.lr[,"R_s"] - data.lr[,"R_q"]
data.lr[,"Rel_Benefit"] <- data.lr[,"Abs_Benefit"] / data.lr[,"R_s"]
data.lr[,"NNQ"] <- 1 / data.lr[,"Abs_Benefit"]
data.lr[data.lr$NNQ < 1,"NNQ"] <- 1
data.lr[data.lr$NNQ > 9999,"NNQ"] <- 9999
data.lr[data.lr$NNQ == Inf,"NNQ"] <- 9999
data.lr[,"Abs_Benefit_per_Qday"] <- data.lr[,"Abs_Benefit"] / data.lr[,"obs_to_iso_q"]
data.lr$d_inf <- params.set[,"d_inf"]
data.lr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.lr$R_0_input <- params.set[,"R_0"]
data.lr$T_lat_offset <- params.set[,"T_lat_offset"]
data.lr$dispersion <- params.set[,"dispersion"]
data.lr$gamma <- params.set[,"gamma"]

#### Combine Datasets ####
data.hr$Setting <- "HR"
data.mr$Setting <- "MR"
data.lr$Setting <- "LR"
data.hr.mr.lr2 <- rbind(data.hr, data.mr, data.lr)
data.hr.mr.lr2$Setting <- factor(data.hr.mr.lr2$Setting, levels = c("HR", "MR", "LR"))


data.hr.mr.lr2$Rel_Benefit_positive <- data.hr.mr.lr2$Rel_Benefit
data.hr.mr.lr2[data.hr.mr.lr2$Rel_Benefit_positive < 0, "Rel_Benefit_positive"] <- 0

data.hr.mr.lr2$Rel_Benefit_per_Qday <- data.hr.mr.lr2$Rel_Benefit / data.hr.mr.lr2$obs_to_iso_q
data.hr.mr.lr2$Rel_Benefit_per_Qday_positive <- data.hr.mr.lr2$Rel_Benefit_positive / data.hr.mr.lr2$obs_to_iso_q

data.hr.mr.lr2$Setting <- factor(data.hr.mr.lr$Setting, levels = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), ordered = TRUE)

#### Plot gamma by Rel_Benefit ####
ggplot(data.hr.mr.lr2, aes(x=gamma, y=Rel_Benefit*100, color = Setting, fill = Setting)) + 
  theme_bw() +  
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence  ", frac(R[S]-R[Q],R[S]), sep=""))) +
  # scale_y_continuous(breaks = seq(0, 30, by=5), labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr2, aes(x=gamma, y=Rel_Benefit_positive*100, color = Setting, fill = Setting)) + 
  theme_bw() +  
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence  ", frac(R[S]-R[Q],R[S]), sep=""))) +
  scale_y_continuous(breaks = seq(0, 30, by=5), labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

#### Plot gamma by Rel_Benefit_per_Qday ####
ggplot(data.hr.mr.lr2, aes(x=gamma, y=Rel_Benefit_per_Qday*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence per Quarantine Day  ", (frac(R[S]-R[Q],R[S]))/d[Q], sep=""))) +
  scale_y_continuous(breaks = seq(0, 3, by=.5), labels = c("0%", "0.5%", "1%", "1.5%", "2%", "2.5%", "3%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr2, aes(x=gamma, y=Rel_Benefit_per_Qday_positive*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Relative Diffrence per Quarantine Day  ", (frac(R[S]-R[Q],R[S]))/d[Q], sep=""))) +
  scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by=1), labels = c("0%", "1%", "2%", "3%", "4%", "5%", "6%", "7%", "8%", "9%", "10%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

#### Plot gamma by other things ####
ggplot(data.hr.mr.lr2[data.hr.mr.lr2$NNQ < 150,], aes(x=gamma, y=R_s, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Symptom Monitoring (", R[S], ")", sep=""))) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)+
  geom_point()

ggplot(data.hr.mr.lr2[data.hr.mr.lr2$NNQ < 150,], aes(x=gamma, y=R_q, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Quarantine (", R[Q], ")", sep="")))+    
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) +
  geom_point()

ggplot(data.hr.mr.lr2[data.hr.mr.lr2$NNQ < 150,], aes(x=gamma, y=Abs_Benefit, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Absolute Difference ", (R[S]-R[Q]), sep=""))) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr2[data.hr.mr.lr2$NNQ < 150,], aes(x=gamma, y=Abs_Benefit_per_Qday, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab(expression(paste("Isolation Effectiveness (", gamma, ")", sep=""))) +
  ylab(expression(paste("Absolute Diffrence per Quarantine Day  ", (frac(R[S]-R[Q],d[Q])), sep=""))) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr2[data.hr.mr.lr2$NNQ < 150,], aes(x=gamma, y=NNQ, color = Setting, fill = Setting)) + 
  theme_bw() + 
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab("Isolation Effectiveness") +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) + geom_point() + facet_grid(Setting~.)

ggplot(data.hr.mr.lr2[data.hr.mr.lr2$NNQ < 150,], aes(x=gamma, y=obs_to_iso_q, color = Setting, fill = Setting)) + 
  theme_bw() + 
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab("Isolation Effectiveness") +
  ylab(expression(paste("Quarantine Days  ", (d[Q]), sep=""))) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) + geom_point()

#### Save Workspace ####
root <- "20151024_Ebola"
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_IntEffect.RData", sep=""))

#### Run model with varying Prob_CT in HR, MR, and LR ####
desired_root <- "20151024_Ebola"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

# High Resource everything else
# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 100
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.hr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  prob_CT = runif(n = times, min = 0, max = 1),
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1.72, max = 1.94)) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  prob_CT <- as.numeric(params.set[i, "prob_CT"])
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c("s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == "s"){
      data.hr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.hr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.hr[,"Abs_Benefit"] <- data.hr[,"R_s"] - data.hr[,"R_q"]
data.hr[,"Rel_Benefit"] <- data.hr[,"Abs_Benefit"] / data.hr[,"R_s"]
data.hr[,"NNQ"] <- 1 / data.hr[,"Abs_Benefit"]
data.hr[data.hr$NNQ < 1,"NNQ"] <- 1
data.hr[data.hr$NNQ > 9999,"NNQ"] <- 9999
data.hr[data.hr$NNQ == Inf,"NNQ"] <- 9999
data.hr[,"Abs_Benefit_per_Qday"] <- data.hr[,"Abs_Benefit"] / data.hr[,"obs_to_iso_q"]
data.hr$d_inf <- params.set[,"d_inf"]
data.hr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.hr$R_0_input <- params.set[,"R_0"]
data.hr$T_lat_offset <- params.set[,"T_lat_offset"]
data.hr$dispersion <- params.set[,"dispersion"]
data.hr$prob_CT <- params.set[,"prob_CT"]

## Middle Resource ##
# Interventions
background_intervention <- "u"

prob_CT <- 0.75

gamma <- 0.75

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.mr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.mr) <- names

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  prob_CT <- params.set[i,"prob_CT"]
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == "s"){
      data.mr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.mr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.mr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.mr[,"Abs_Benefit"] <- data.mr[,"R_s"] - data.mr[,"R_q"]
data.mr[,"Rel_Benefit"] <- data.mr[,"Abs_Benefit"] / data.mr[,"R_s"]
data.mr[,"NNQ"] <- 1 / data.mr[,"Abs_Benefit"]
data.mr[data.mr$NNQ < 1,"NNQ"] <- 1
data.mr[data.mr$NNQ > 9999,"NNQ"] <- 9999
data.mr[data.mr$NNQ == Inf,"NNQ"] <- 9999
data.mr[,"Abs_Benefit_per_Qday"] <- data.mr[,"Abs_Benefit"] / data.mr[,"obs_to_iso_q"]
data.mr$d_inf <- params.set[,"d_inf"]
data.mr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.mr$R_0_input <- params.set[,"R_0"]
data.mr$T_lat_offset <- params.set[,"T_lat_offset"]
data.mr$dispersion <- params.set[,"dispersion"]
data.mr$prob_CT <- params.set[,"prob_CT"]

## Low Resource ##

# Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

parms_epsilon = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.lr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.lr) <- names

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  prob_CT <- params.set[i,"prob_CT"]
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_inf$parm2 <- params.set[i,"d_inf"]
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == "s"){
      data.lr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.lr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.lr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.lr[,"Abs_Benefit"] <- data.lr[,"R_s"] - data.lr[,"R_q"]
data.lr[,"Rel_Benefit"] <- data.lr[,"Abs_Benefit"] / data.lr[,"R_s"]
data.lr[,"NNQ"] <- 1 / data.lr[,"Abs_Benefit"]
data.lr[data.lr$NNQ < 1,"NNQ"] <- 1
data.lr[data.lr$NNQ > 9999,"NNQ"] <- 9999
data.lr[data.lr$NNQ == Inf,"NNQ"] <- 9999
data.lr[,"Abs_Benefit_per_Qday"] <- data.lr[,"Abs_Benefit"] / data.lr[,"obs_to_iso_q"]
data.lr$d_inf <- params.set[,"d_inf"]
data.lr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.lr$R_0_input <- params.set[,"R_0"]
data.lr$T_lat_offset <- params.set[,"T_lat_offset"]
data.lr$dispersion <- params.set[,"dispersion"]
data.lr$prob_CT <- params.set[,"prob_CT"]

#### Combine Datasets ####
data.hr$Setting <- "HR"
data.mr$Setting <- "MR"
data.lr$Setting <- "LR"
data.hr.mr.lr3 <- rbind(data.hr, data.mr, data.lr)
data.hr.mr.lr3$Setting <- factor(data.hr.mr.lr3$Setting, levels = c("HR", "MR", "LR"))

View(data.hr.mr.lr3[is.na(data.hr.mr.lr3$obs_to_iso_q)==1,])

data.hr.mr.lr3$Rel_Benefit_positive <- data.hr.mr.lr3$Rel_Benefit
data.hr.mr.lr3[data.hr.mr.lr3$Rel_Benefit_positive < 0, "Rel_Benefit_positive"] <- 0
data.hr.mr.lr3[data.hr.mr.lr3$Rel_Benefit_positive > 10, "Rel_Benefit_positive"]

data.hr.mr.lr3$Rel_Benefit_per_Qday <- data.hr.mr.lr3$Rel_Benefit / data.hr.mr.lr3$obs_to_iso_q
data.hr.mr.lr3$Rel_Benefit_per_Qday_positive <- data.hr.mr.lr3$Rel_Benefit_positive / data.hr.mr.lr3$obs_to_iso_q

data.hr.mr.lr3$Setting <- factor(data.hr.mr.lr$Setting, levels = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), ordered = TRUE)

#### Plot prob_CT by Rel_Benefit ####
ggplot(data.hr.mr.lr3, aes(x=prob_CT, y=Rel_Benefit*100, color = Setting, fill = Setting)) + 
  theme_bw() +  
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab("Fraction of Contacts Traced") +
  ylab("Percent reduction in Re by quarantine over symptom monitoring") +
  scale_y_continuous(breaks = seq(0, 30, by=5), labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr3, aes(x=prob_CT, y=Rel_Benefit_positive*100, color = Setting, fill = Setting)) + 
  theme_bw() +  
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab("Fraction of Contacts Traced") +
  ylab("Percent reduction in Re by quarantine over symptom monitoring") +
  scale_y_continuous(breaks = seq(0, 30, by=5), labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) + geom_point()

#### Plot prob_CT by Rel_Benefit_per_Qday ####
ggplot(data.hr.mr.lr3, aes(x=prob_CT, y=Rel_Benefit_per_Qday*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab("Fraction of Contacts Traced") +
  ylab("Percent reduction in Re by quarantine over symptom monitoring per quarantine day") +
  scale_y_continuous(breaks = seq(0, 3, by=.5), labels = c("0%", "0.5%", "1%", "1.5%", "2%", "2.5%", "3%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE)

ggplot(data.hr.mr.lr3, aes(x=prob_CT, y=Rel_Benefit_per_Qday_positive*100, color = Setting, fill = Setting)) + 
  theme_bw() +
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  xlab("Fraction of Contacts Traced") +
  ylab("Percent reduction in Re by quarantine over symptom monitoring per quarantine day") +
  # scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by=1), labels = c("0%", "1%", "2%", "3%", "4%", "5%", "6%", "7%", "8%", "9%", "10%")) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) + geom_point()

#### Plot prob_CT by other things
ggplot(data.hr.mr.lr3, aes(x=prob_CT, y=obs_to_iso_q, color = Setting, fill = Setting)) + 
  theme_bw() + 
  scale_color_manual(values = gbr_cols, breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  scale_fill_manual(values = alpha(gbr_cols, 0.2), breaks = c("HR", "MR", "LR"), labels = c("High Resource", "Mid Resource", "Low Resource"), name = "Setting") +
  ggtitle("") +
  xlab(expression(paste("Fraction of Contacts Traced ", (P[CT]), sep=""))) +
  ylab(expression(paste("Quarantine Days  ", (d[Q]), sep=""))) +
  geom_smooth(size = 2,alpha = 0.2) +
  geom_smooth(size = 2, se=FALSE) + geom_point()

#### Save Workspace ####
root <- "20151024_Ebola"
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_IntEffect.RData", sep=""))
