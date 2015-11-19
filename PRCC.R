#### Header ####
# Extracted code for PRCC calculation and plots

#### Load Libraries ####
library(ggplot2)
library(RColorBrewer)
library(sensitivity)
library(reshape)
library(lhs)

#### Source Functions ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

#### Load Workspaces ####
desired_root <- "20151118_SARS" # Paste the desired root here "YYYYMMDD_DISEASE"

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/",  desired_root, "_PRCC.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))

desired_date <- "20151111"
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_date, "_PRCC.RData", sep=""))

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
times <- 200
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "Rel_Benefit_per_Qday","Rel_Benefit_per_Qday_rp", "NQD")
data.prcc <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.prcc) <- names

# Sample from posterior distributions for each parameter independently
params.set.prcc <- cbind(
  T_lat_offset = sample(data$T_lat_offset, size = times, replace = TRUE),
  d_inf = sample(data$d_inf, size = times, replace = TRUE),
  pi_t_triangle_center = sample(data$pi_t_triangle_center, size = times, replace = TRUE) )

# Set range for other parameters to vary
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","riskprofile", "R_0_mean","R_0_spread","dispersion","T_inc_stretch")
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0
gamma.max <- 1
prob_CT.min <- 0
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 7
epsilon.min <- 0
epsilon.max <- 7
riskprofile.min <- 0.01
riskprofile.max <- 1
R_0_mean.min <- 1
R_0_mean.max <- 5
R_0_spread.min <- 0
R_0_spread.max <- 0.9
dispersion.min <- 1
dispersion.max <- 4
T_inc_stretch.min <- 0.5
T_inc_stretch.max <- 1.5

params.set.prcc <- cbind(params.set.prcc,
                         gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
                         prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
                         CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
                         epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
                         riskprofile = lhs[,5]*(riskprofile.max - riskprofile.min) + riskprofile.min,
                         R_0_mean = lhs[,6]*(R_0_mean.max - R_0_mean.min) + R_0_mean.min,
                         R_0_spread = lhs[,7]*(R_0_spread.max - R_0_spread.min) + R_0_spread.min,
                         dispersion = lhs[,8]*(dispersion.max - dispersion.min) + dispersion.min,
                         T_inc_stretch = lhs[,9]*(T_inc_stretch.max - T_inc_stretch.min) + T_inc_stretch.min)
params.set.prcc <- data.frame(params.set.prcc)

source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")
i=1
while (i <= times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
  gamma <- as.numeric(params.set.prcc[i,"gamma"])
  prob_CT <- as.numeric(params.set.prcc[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set.prcc[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set.prcc[i,"epsilon"])
  riskprofile <- as.numeric(params.set.prcc[i, "riskprofile"])
    R_0_input.min <- as.numeric(params.set.prcc[i,"R_0_mean"]) - as.numeric(params.set.prcc[i,"R_0_mean"])*as.numeric(params.set.prcc[i,"R_0_spread"])
    R_0_input.max <- as.numeric(params.set.prcc[i,"R_0_mean"]) + as.numeric(params.set.prcc[i,"R_0_mean"])*as.numeric(params.set.prcc[i,"R_0_spread"])
  parms_R_0[c("parm1","parm2")] <- c( R_0_input.min, R_0_input.max)
  parms_pi_t$triangle_center <- as.numeric(params.set.prcc[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set.prcc[i,"T_lat_offset"])
  parms_d_inf$parm2 <- as.numeric(params.set.prcc[i,"d_inf"])
  dispersion <- as.numeric(params.set.prcc[i, "dispersion"])
  parms_T_inc$T_inc_stretch <- as.numeric(params.set.prcc[i,"T_inc_stretch"])
  
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

data.prcc[,"Abs_Benefit"] <- data.prcc[,"R_s"] - data.prcc[,"R_q"]
data.prcc[,"NNQ"] <- 1 / data.prcc[,"Abs_Benefit"]
data.prcc[data.prcc$obs_to_iso_q < 0.01, "obs_to_iso_q"] <- 0.01
data.prcc[,"Abs_Benefit_per_Qday"] <- data.prcc[,"Abs_Benefit"] / data.prcc[,"obs_to_iso_q"]
data.prcc[,"NQD"] <- 1 / data.prcc[,"Abs_Benefit_per_Qday"]
data.prcc[data.prcc$NNQ < 1,"NNQ"] <- 1
data.prcc[data.prcc$NNQ > 9999,"NNQ"] <- 9999
data.prcc[data.prcc$NNQ == Inf,"NNQ"] <- 9999

data.prcc[,"Rel_Benefit"] <- data.prcc[,"Abs_Benefit"] / data.prcc[,"R_s"]
data.prcc[,"Rel_Benefit_per_Qday"] <- data.prcc[,"Rel_Benefit"] / data.prcc[,"obs_to_iso_q"]

data.prcc$T_lat_offset <- params.set.prcc[,"T_lat_offset"]
data.prcc$d_inf <- params.set.prcc[,"d_inf"]
data.prcc$pi_t_triangle_center <- params.set.prcc[,"pi_t_triangle_center"]
data.prcc$gamma <- params.set.prcc[,"gamma"]
data.prcc$prob_CT <- params.set.prcc[,"prob_CT"]
data.prcc$CT_delay <- params.set.prcc[,"CT_delay"]
data.prcc$epsilon <- params.set.prcc[,"epsilon"]
data.prcc$riskprofile <- params.set.prcc[,"riskprofile"]
data.prcc[,"R_0_mean"] <- params.set.prcc[,"R_0_mean"]
data.prcc[,"R_0_spread"] <- params.set.prcc[,"R_0_spread"]
data.prcc$dispersion <- params.set.prcc[,"dispersion"]
data.prcc$T_inc_stretch <- params.set.prcc[,"T_inc_stretch"]

#### Add a Rel_Benefit_per_Qday_rp that considers Risk Profiling ####
disease <- substr(root, 10, nchar(root)) # Find upper 95 percerntile for incubation period for each disease
cat(disease)
if (disease == "Ebola"){
  T_inc_95 <- 23.80
} else if (disease == "HepatitisA"){
  T_inc_95 <- 33.30
} else if (disease == "InfluenzaA"){
  T_inc_95 <- 2.73
} else if (disease == "MERS"){
  T_inc_95 <- 12.46
} else if (disease == "Pertussis"){
  T_inc_95 <- 9.51
} else if (disease == "SARS"){
  T_inc_95 <- 10.62
} else if (disease == "Smallpox"){
  T_inc_95 <- 15.64
} 
data.prcc$Rel_Benefit_per_Qday_rp <- data.prcc$Rel_Benefit_per_Qday / ( data.prcc$obs_to_iso_q + (1/data.prcc$riskprofile - 1)*(T_inc_95))

#### Confirm Monotonicity ####
# Plot each of the covariate - outcome scatterplots
for (covariate in c("gamma", "prob_CT","CT_delay","epsilon", "riskprofile", "R_0_mean", "R_0_spread", "dispersion","T_inc_stretch", "T_lat_offset")){
    panel_plot_fcn(data = data.prcc, covariate = covariate, outputs = c("R_0", "R_s", "R_q", "Rel_Benefit", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp"))
    cat ("Press [enter] to continue")
    line <- readline()
}

# Compare ks value across deciles of covariates to make sure fit isn't horrendous
decile_plot_fcn(data.prcc, params.set.prcc)

#### Calculate PRCC for one disease ####
dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "Rel_Benefit","obs_to_iso_q","ks", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon", "riskprofile", "R_0_mean", "R_0_spread","dispersion","T_inc_stretch","pi_t_triangle_center","T_lat_offset","d_inf")
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")
output <- prcc_fcn(input_data = data.prcc, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

#### See if PRCC depends on parameter ranges ####
nrow(data.prcc[data.prcc$gamma > 0.6,]) 
output.1 <- prcc_fcn(input_data = data.prcc[data.prcc$gamma > 0.6,], dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

nrow(data.prcc[data.prcc$gamma < 0.6,])
output.2 <- prcc_fcn(input_data = data.prcc[data.prcc$gamma < 0.6,], dep_var = dep_var, indep_var = indep_var, 
                     nboot = 100, package = "sensitivity", standardize = TRUE)

nrow(data.prcc[data.prcc$epsilon < 3,])
output.3 <- prcc_fcn(input_data = data.prcc[data.prcc$epsilon < 3,], dep_var = dep_var, indep_var = indep_var, 
                     nboot = 100, package = "sensitivity", standardize = TRUE)

nrow(data.prcc[data.prcc$T_lat_offset < 3,])
output.3 <- prcc_fcn(input_data = data.prcc[data.prcc$epsilon < 3,], dep_var = dep_var, indep_var = indep_var, 
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

# pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plot_prcc.pdf", sep=""))
# plot(plot_prcc_1)
# dev.off()

#### Save Workspace ####
cat(desired_root)
# save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_PRCC.RData", sep=""))

#### Load PRCC data for multiple diseases ####
# desired_roots_list <- c("20151022_SARS", "20151024_Ebola", "20151026_HepatitisA", "20151026_Pertussis", "20151027_MERS", "20151028_InfluenzaA", "20151028_Smallpox")
desired_roots_list <- c("20151118_SARS", "20151118_Ebola", "20151118_HepatitisA", "20151118_Pertussis", "20151118_MERS", "20151118_InfluenzaA", "20151118_Smallpox")
# desired_roots_list <- c("20151022_SARS", "20151028_InfluenzaA")

# load the first one to get a list of the headers
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_roots_list[1], "/", desired_roots[1], "_PRCC.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_roots_list[1], "_PRCC.RData", sep=""))

names <- c(names(data.prcc), "disease")
length(names)
df.prcc.temporary <- data.frame(matrix(rep(NA, length(names)), nrow=1))
names(df.prcc.temporary) <- names

names <- c(names(output), "disease")
length(names)
df.prcc.output.temporary <- data.frame(matrix(rep(NA, length(names)), nrow=1))
names(df.prcc.output.temporary) <- names

for (i.temporary in 1:length(desired_roots_list)){
  desired_roots_list <- c("20151118_SARS", "20151118_Ebola", "20151118_HepatitisA", "20151118_Pertussis", "20151118_MERS", "20151118_InfluenzaA", "20151118_Smallpox")
  root_i <- desired_roots_list[i.temporary]
  # load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root_i, "/", root_i, "_PRCC.RData", sep=""))
  load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root_i, "_PRCC.RData", sep=""))
  
  data.prcc$disease <- substr(root_i, 10, nchar(root_i))
  output$disease <- substr(root_i, 10, nchar(root_i))
  
  # Add an variable and outcome for risk profiling
  if (output$disease[1] == "Ebola"){
    T_inc_95 <- 23.80
  } else if (output$disease[1] == "HepatitisA"){
    T_inc_95 <- 33.30
  } else if (output$disease[1] == "InfluenzaA"){
    T_inc_95 <- 2.73
  } else if (output$disease[1] == "MERS"){
    T_inc_95 <- 12.46
  } else if (output$disease[1] == "Pertussis"){
    T_inc_95 <- 9.51
  } else if (output$disease[1] == "SARS"){
    T_inc_95 <- 10.62
  } else if (output$disease[1] == "Smallpox"){
    T_inc_95 <- 15.64
  } 

  cat("\n",length(names(df.prcc.temporary)))
  cat("\n",length(names(df.prcc.output.temporary)))
  
  df.prcc.temporary <- rbind(df.prcc.temporary, data.prcc)
  df.prcc.output.temporary <- rbind(df.prcc.output.temporary, output)
}

df.prcc.temporary <- df.prcc.temporary[is.na(df.prcc.temporary$R_0)==0,]
df.prcc.output.temporary <- df.prcc.output.temporary[is.na(df.prcc.output.temporary$coef)==0,]

#### Save Workspace ####
cat(date)
date <- "20151119"
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PRCC.RData", sep=""))

#### Calculate PRCC for many diseases together ####
dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "NQD", "Rel_Benefit","obs_to_iso_q","ks", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon", "riskprofile", "R_0_mean", "R_0_spread","dispersion","T_inc_stretch","pi_t_triangle_center","T_lat_offset","d_inf")
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")
output.all <- prcc_fcn(input_data = df.prcc.temporary, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

# Add output.all to the df.prcc.output file so an all-diseases bar can be added to the horizontal grouped bar chart
output.all$disease <- "all"
df.prcc.output.temporary <- rbind(df.prcc.output.temporary, output.all)

#### Plot all diseases ####
plot_prcc_2 <- ggplot(output.all, aes(x = parameter, y= coef)) +
  facet_grid(output ~ .) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  geom_errorbar(data = output.all, aes(ymin = CImin, ymax = CImax), width = 0.1) +
  ggtitle("All diseases")
plot_prcc_2

cat(date)
pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_Plot_prcc_2.pdf", sep=""))
plot(plot_prcc_2)
dev.off()

# Note that obs_to_iso_q is NOT monotonic for the pooled estimate for T_lat_offset, pi_t_triangle_center, d_inf, and gamma.
for (covariate in indep_var){
  panel_plot_fcn(data = df.prcc.temporary, covariate = covariate, outputs = c("obs_to_iso_q", "Rel_Benefit", "Rel_Benefit_per_Qday"))
  cat ("Press [enter] to continue")
  line <- readline()
}

#### Plot all diseases, horizontal bar chart ####
plot_prcc_3 <- ggplot(output.all, aes(x = parameter, y = coef)) +
  facet_grid(.~output) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  ggtitle("All diseases") +
  coord_flip()
plot_prcc_3

#### Plot each disease, horizontal bar chart, subset of outputs ####
df.prcc.outut.subset <- NA
df.prcc.output.subset <- df.prcc.output.temporary[is.element(df.prcc.output.temporary$output, c("R_s", "R_q", "Rel_Benefit", "Rel_Benefit_per_Qday_rp")),]
df.prcc.output.subset$parameter <- factor(df.prcc.output.subset$parameter, levels = rev(c("gamma", "prob_CT", "riskprofile", "epsilon", "CT_delay", "R_0_mean", "R_0_spread", "T_inc_stretch", "T_lat_offset", "pi_t_triangle_center", "d_inf", "dispersion")), ordered = TRUE,
                                          labels = rev(c("Isolation\nEffectiveness","Fraction of\n Contacts Traced", "Fraction of Traced\nContacts who are Infected", "Delay from Symptom Onset\nto Isolation", "Delay in Tracing\na Contact", "R_0 Mean", "R_0 Spread", "Incubation\nPeriod", "Latent\nPeriod", "Time of Peak\nInfectiousness", "Duration of\nInfectiousness", "Super-Spreading\nDispersion Factor")))
df.prcc.output.subset$output <- factor(df.prcc.output.subset$output, levels = c("R_s","R_q", "Rel_Benefit", "Rel_Benefit_per_Qday_rp"), ordered = TRUE,
                                       labels = c("Symptom Monitoring R", "Quarantine R", "Relative Difference", "Relative Difference\nper Quarantine Day"))
df.prcc.output.subset$disease <- factor(df.prcc.output.subset$disease, levels = rev(c("Ebola","HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox", "all")), ordered = TRUE, labels = rev(c("Ebola","Hepatitis A", "Influenza A", "MERS", "Pertussis", "SARS", "Smallpox", "All Diseases")))

scale_colour_brewer(type="qual", palette=6)
my.cols <- brewer.pal(n = 7, name = "Set1")
my.cols <- c(my.cols[c(3, 7, 4, 5, 1, 6, 2)], "black")

plot_prcc_4 <- ggplot(df.prcc.output.subset, aes(x = parameter, y = coef, fill = parameter, color = disease)) +
  facet_grid(.~output) +
  geom_bar(position = "dodge", stat="identity", width = .9) +
  ggtitle("Partial Rank Correlation Coefficient") +
  coord_flip() +
  # scale_fill_manual(values = rep(c("grey", "black"), 4)) +
  scale_fill_brewer(type = "div", palette = 6) +
  scale_color_manual(values = rev(my.cols), breaks = c("Ebola","HepatitisA", "InfluenzaA", "MERS", "Pertussis", "SARS", "Smallpox", "all")) +
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

