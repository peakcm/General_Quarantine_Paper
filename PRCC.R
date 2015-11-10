#### Header ####
# Extracted code for PRCC calculation and plots

#### Load Libraries ####
library(ggplot2)
library(RColorBrewer)

#### Source Functions ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

#### Load Workspaces ####
desired_root <- "20151024_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_PRCC.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))

#### Pull from SMC data with independent draws ####

# Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Initialize
n_pop = 500
num_generations <- 5
times <- 1000
names <- c("R_0_input", "R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "Rel_Benefit_per_Qday", "NQD")
data.prcc <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.prcc) <- names

# Sample from posterior distributions for each parameter independently
params.set.prcc <- cbind(
  T_lat_offset = sample(data$T_lat_offset, size = times, replace = TRUE),
  d_inf = sample(data$d_inf, size = times, replace = TRUE),
  pi_t_triangle_center = sample(data$pi_t_triangle_center, size = times, replace = TRUE) )

# Set range for other parameters to vary
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0", "dispersion")
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0
gamma.max <- 1
prob_CT.min <- 0
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 7
epsilon.min <- 0
epsilon.max <- 7
R_0.min <- 5
R_0.max <- 5
dispersion.min <- 1
dispersion.max <- 4

params.set.prcc <- cbind(params.set.prcc,
                         gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
                         prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
                         CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
                         epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
                         dispersion = lhs[,6]*(dispersion.max - dispersion.min) + dispersion.min,
                         R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min)
params.set.prcc <- data.frame(params.set.prcc)

i=1
while (i <= times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  gamma <- as.numeric(params.set.prcc[i,"gamma"])
  prob_CT <- as.numeric(params.set.prcc[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set.prcc[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set.prcc[i,"epsilon"])
  parms_pi_t$triangle_center <- as.numeric(params.set.prcc[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set.prcc[i,"T_lat_offset"])
  parms_d_inf$parm2 <- as.numeric(params.set.prcc[i,"d_inf"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set.prcc[i,"R_0"]), as.numeric(params.set.prcc[i,"R_0"]))
  dispersion <- as.numeric(params.set.prcc[i, "dispersion"])
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 200
    } else {n_pop_input <- n_pop}
    
    In_Out <- repeat_call_fcn(n_pop=n_pop_input, 
                              parms_T_inc = parms_T_inc, 
                              parms_T_lat = parms_T_lat, 
                              parms_d_inf = parms_d_inf, 
                              parms_d_symp = parms_d_symp, 
                              parms_R_0 = parms_R_0, 
                              parms_epsilon = parms_epsilon, 
                              parms_pi_t = parms_pi_t,
                              num_generations = num_generations,
                              background_intervention = background_intervention,
                              subseq_interventions = subseq_interventions,
                              gamma = gamma,
                              prob_CT = prob_CT,
                              parms_CT_delay = parms_CT_delay,
                              parms_serial_interval = parms_serial_interval,
                              dispersion = dispersion,
                              printing = printing)
    if (subseq_interventions == background_intervention){
      data.prcc[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.prcc[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.prcc[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.prcc[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.prcc[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.prcc[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
  
  if (is.na(data.prcc[i,"R_0"])==1 |
      is.na(data.prcc[i,"R_hsb"])==1 |
      is.na(data.prcc[i,"R_s"])==1 |
      is.na(data.prcc[i,"R_q"])==1){
    i=i  #re-run that set
    cat("x")
  } else {i = i+1}
  
}

# Check for missing data
if (sum(is.na(data.prcc[,1:4]))>0){cat("Something's missing")}

data.prcc$R_0_input <- params.set.prcc[,"R_0"]
data.prcc[,"Abs_Benefit"] <- data.prcc[,"R_s"] - data.prcc[,"R_q"]
data.prcc[,"NNQ"] <- 1 / data.prcc[,"Abs_Benefit"]
data.prcc[,"Abs_Benefit_per_Qday"] <- data.prcc[,"Abs_Benefit"] / data.prcc[,"obs_to_iso_q"]
data.prcc[,"NQD"] <- 1 / data.prcc[,"Abs_Benefit_per_Qday"]
data.prcc[data.prcc$NNQ < 1,"NNQ"] <- 1
data.prcc[data.prcc$NNQ > 9999,"NNQ"] <- 9999
data.prcc[data.prcc$NNQ == Inf,"NNQ"] <- 9999

data.prcc[,"Rel_Benefit"] <- data.prcc[,"Abs_Benefit"] / data.prcc[,"R_s"]
data.prcc[,"Rel_Benefit_per_Qday"] <- data.prcc[,"Rel_Benefit"] / data.prcc[,"obs_to_iso_q"]

data.prcc$gamma <- params.set.prcc[,"gamma"]
data.prcc$prob_CT <- params.set.prcc[,"prob_CT"]
data.prcc$CT_delay <- params.set.prcc[,"CT_delay"]
data.prcc$epsilon <- params.set.prcc[,"epsilon"]
data.prcc$dispersion <- params.set.prcc[,"dispersion"]
data.prcc$pi_t_triangle_center <- params.set.prcc[,"pi_t_triangle_center"]
data.prcc$T_lat_offset <- params.set.prcc[,"T_lat_offset"]
data.prcc$d_inf <- params.set.prcc[,"d_inf"]

#### Confirm Monotonicity ####
# Plot each of the covariate - outcome scatterplots
for (covariate in c("gamma","prob_CT","CT_delay","epsilon","R_0", "dispersion")){
# for (covariate in c("gamma","prob_CT","CT_delay","epsilon","R_0_input", "dispersion")){
    panel_plot_fcn(data = data.prcc, covariate = covariate)
    cat ("Press [enter] to continue")
    line <- readline()
}

# Compare ks value across deciles of covariates to make sure fit isn't horrendous
decile_plot_fcn(data.prcc, params.set.prcc)

#### Calculate PRCC for one disease ####
dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "Rel_Benefit","obs_to_iso_q","ks")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon","dispersion","pi_t_triangle_center","T_lat_offset","d_inf")
# dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "NQD", "Rel_Benefit", "Rel_Benefit_per_Qday","obs_to_iso_q","ks")
# indep_var <- c("gamma","prob_CT","CT_delay", "epsilon","dispersion","pi_t_triangle_center","T_lat_offset","d_inf", "R_0_input")
output <- prcc_fcn(input_data = data.prcc, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

#### Plot each disease ####
plot_prcc_1 <- ggplot(output, aes(x = parameter, y= coef)) +
  facet_grid(output ~ .) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  ggtitle(desired_root) +
  geom_errorbar(data = output, aes(ymin = CImin, ymax = CImax), width = 0.1)
plot_prcc_1

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot_prcc.pdf", sep=""))
plot(plot_prcc_1)
dev.off()

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PRCC.RData", sep=""))

#### Load PRCC data for multiple diseases ####

desired_roots <- c("20151022_SARS", "20151024_Ebola", "20151026_HepatitisA", "20151026_Pertussis", "20151027_MERS", "20151028_InfluenzaA", "20151028_Smallpox")

# load the first one to get a list of the headers
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_roots[1], "/", desired_roots[1], "_PRCC.RData", sep=""))

names <- c(names(data.prcc), "disease")
df.prcc <- data.frame(matrix(rep(NA, length(names)), nrow=1))
names(df.prcc) <- names

names <- c(names(output), "disease")
df.prcc.output <- data.frame(matrix(rep(NA, length(names)), nrow=1))
names(df.prcc.output) <- names

for (ind in desired_roots){
  load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", ind, "/", ind, "_PRCC.RData", sep=""))
  data.prcc$disease <- substr(ind, 10, nchar(ind))
  output$disease <- substr(ind, 10, nchar(ind))
  
  cat("\n",names(df.prcc))
  cat("\n",names(df.prcc.output))
  
  df.prcc <- rbind(df.prcc, data.prcc)
  df.prcc.output <- rbind(df.prcc.output, output)
  names(df.prcc.output)
}
df.prcc <- df.prcc[is.na(df.prcc$R_0)==0,]
df.prcc.output <- df.prcc.output[is.na(df.prcc.output$coef)==0,]

#### Calculate PRCC for many diseases together ####
dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "Rel_Benefit","obs_to_iso_q","ks")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon","dispersion","pi_t_triangle_center","T_lat_offset","d_inf")
# dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "NQD", "Rel_Benefit", "Rel_Benefit_per_Qday","obs_to_iso_q","ks")
# indep_var <- c("gamma","prob_CT","CT_delay", "epsilon","dispersion","pi_t_triangle_center","T_lat_offset","d_inf", "R_0_input")
output.all <- prcc_fcn(input_data = df.prcc, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

#### Plot all diseases ####
plot_prcc_2 <- ggplot(output.all, aes(x = parameter, y= coef)) +
  facet_grid(output ~ .) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  geom_errorbar(data = output.all, aes(ymin = CImin, ymax = CImax), width = 0.1) +
  ggtitle("All diseases")
plot_prcc_2

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_Plot_prcc_2.pdf", sep=""))
plot(plot_prcc_2)
dev.off()

#### Plot all diseases, horizontal bar chart ####
plot_prcc_3 <- ggplot(output.all, aes(x = parameter, y = coef)) +
  facet_grid(.~output) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  ggtitle("All diseases") +
  coord_flip()
plot_prcc_3

#### Plot each disease, horizontal bar chart ####
df.prcc.output.subset <- df.prcc.output[is.element(df.prcc.output$output, c("R_s", "R_q", "Abs_Benefit", "Rel_Benefit")),]
df.prcc.output.subset$parameter <- factor(df.prcc.output.subset$parameter, levels = rev(c("gamma", "prob_CT", "epsilon", "CT_delay", "T_lat_offset", "pi_t_triangle_center", "d_inf", "dispersion")), ordered = TRUE)
df.prcc.output.subset$output <- factor(df.prcc.output.subset$output, levels = c("R_s","R_q", "Abs_Benefit", "Rel_Benefit"), ordered = TRUE)
df.prcc.output.subset$disease <- factor(df.prcc.output.subset$disease, levels = rev(c("Ebola","HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox")), ordered = TRUE)

scale_colour_brewer(type="qual", palette=6)
my.cols <- brewer.pal(n = 7, name = "Set1")
my.cols <- my.cols[c(3, 7, 4, 5, 1, 6, 2)]

plot_prcc_4 <- ggplot(df.prcc.output.subset, aes(x = parameter, y = coef, fill = parameter, color = disease)) +
  facet_grid(.~output) +
  geom_bar(position = "dodge", stat="identity", width = .9) +
  ggtitle("Partial Rank Correlation Coefficient") +
  coord_flip() +
  # scale_fill_manual(values = rep(c("grey", "black"), 4)) +
  scale_fill_brewer(type = "div", palette = 6) +
  scale_color_manual(values = rev(my.cols), breaks = c("Ebola","HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox")) +
  theme_bw() +
  theme(axis.title.y=element_blank()) +
  ylab("Partial Rank Correlation Coefficient") +
  geom_hline(yintercept=0, color="darkgrey", size=1) +
  guides(color=FALSE) +
  guides(fill=FALSE)
plot_prcc_4

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_Plot_prcc_4.pdf", sep=""), width=9, height=6)
plot(plot_prcc_4)
dev.off()

 #### Save Workspace ####
date <- format(Sys.time(), "%Y%m%d")
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PRCC.RData", sep=""))

