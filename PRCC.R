#### Header ####
# Extracted code for PRCC calculation and plots

#### Load Libraries ####
library(ggplot2)
library(RColorBrewer)
library(sensitivity)
library(reshape)
library(lhs)
library(MASS)

#### Source Functions ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

#### Load Workspaces ####

# To load SMC data for a disease and then generate PRCC outputs for that disease
desired_root <- "20160507_SARS" # Paste the desired root here "YYYYMMDD_DISEASE"
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMC.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))

# To load the PRCC outputs for one disease
desired_root <- "20160507_SARS" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/",  desired_root, "_PRCC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))

# For the pooled PRCC data
desired_date <- "20160507"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_date, "_PRCC.RData", sep=""))

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
times <- 5000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks", "Rel_Benefit_per_Qday","Rel_Benefit_per_Qday_rp", "NQD")
data.prcc <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.prcc) <- names

# Sample from posterior distributions for each parameter independently
params.set.prcc <- cbind(
  T_lat_offset = sample(data$T_lat_offset, size = times, replace = TRUE),
  d_inf = sample(data$d_inf, size = times, replace = TRUE),
  pi_t_triangle_center = sample(data$pi_t_triangle_center, size = times, replace = TRUE) )

# Set range for other parameters to vary
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","riskprofile", "R_0_mean","dispersion","T_inc_stretch")
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0.01
gamma.max <- 1
prob_CT.min <- 0.01
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 7
epsilon.min <- 0
epsilon.max <- 7
riskprofile.min <- 0.01
riskprofile.max <- 1
R_0_mean.min <- 1
R_0_mean.max <- 5
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
                         dispersion = lhs[,7]*(dispersion.max - dispersion.min) + dispersion.min,
                         T_inc_stretch = lhs[,8]*(T_inc_stretch.max - T_inc_stretch.min) + T_inc_stretch.min)
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
  R_0_input.min <- as.numeric(params.set.prcc[i,"R_0_mean"])
  R_0_input.max <- as.numeric(params.set.prcc[i,"R_0_mean"])
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
data.prcc$dispersion <- params.set.prcc[,"dispersion"]
data.prcc$T_inc_stretch <- params.set.prcc[,"T_inc_stretch"]

#### Add a Abs_ and Rel_Benefit_per_Qday_rp that considers Risk Profiling ####
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
data.prcc$Abs_Benefit_per_Qday_rp <- data.prcc$Abs_Benefit_per_Qday / ( data.prcc$obs_to_iso_q + (1/data.prcc$riskprofile - 1)*(T_inc_95))

#### Explore PRCC dataset ####
apply(data.prcc, 2, summary)
nrow(data.prcc)

#### Drop PRCC dataset rows with missing values ####
nrow(data.prcc[is.na(data.prcc$obs_to_iso_q),])
nrow(data.prcc[data.prcc$Rel_Benefit_per_Qday == -Inf,])
nrow(data.prcc[data.prcc$Rel_Benefit_per_Qday == Inf,])
data.prcc[data.prcc$NQD == Inf,]

data.prcc <- data.prcc[is.na(data.prcc$obs_to_iso_q)==0 & 
                         data.prcc$Rel_Benefit_per_Qday != -Inf &
                         data.prcc$Rel_Benefit_per_Qday != Inf &
                         data.prcc$NQD != Inf,]

#### Confirm Monotonicity ####
# Plot each of the covariate - outcome scatterplots
for (covariate in c("gamma", "prob_CT","CT_delay","epsilon", "riskprofile", "R_0_mean", "dispersion","T_inc_stretch", "T_lat_offset")){
    panel_plot_fcn(data = data.prcc, covariate = covariate, outputs = c("R_0", "R_s", "R_q", "Rel_Benefit", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp"))
    cat ("Press [enter] to continue")
    line <- readline()
}

# Compare ks value across deciles of covariates to make sure fit isn't horrendous
decile_plot_fcn(data.prcc, params.set.prcc)

#### Calculate PRCC for one disease ####
dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "Rel_Benefit","obs_to_iso_q","ks", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon", "riskprofile", "R_0_mean","dispersion","T_inc_stretch","pi_t_triangle_center","T_lat_offset","d_inf")
output <- prcc_fcn(input_data = data.prcc, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

#### Add more outcomes to PRCC ####
data.prcc$Impact_SM <- data.prcc$R_0 - data.prcc$R_s
data.prcc$Impact_Q <- data.prcc$R_0 - data.prcc$R_q
data.prcc$Impact_HSB <- data.prcc$R_0 - data.prcc$R_hsb

dep_var <- c("Impact_SM", "Impact_Q", "Impact_HSB", "Abs_Benefit_per_Qday_rp")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon", "riskprofile", "R_0_mean","dispersion","T_inc_stretch","pi_t_triangle_center","T_lat_offset","d_inf")
output_append <- prcc_fcn(input_data = data.prcc, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)
output <- rbind(output, output_append)

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

#### Save Workspace and .csv's ####
# cat(desired_root)
# date <- format(Sys.time(), "%Y%m%d")
# disease <- "SARS"
(root <- paste(date, disease, sep = "_"))
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PRCC.RData", sep=""))
write.csv(data.prcc, paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_dataprcc.csv", sep=""), row.names = FALSE)
write.csv(output, paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_outputprcc.csv", sep=""), row.names = FALSE)

#### Load PRCC data for multiple diseases ####
# desired_roots_list <- c("20151022_SARS", "20151024_Ebola", "20151026_HepatitisA", "20151026_Pertussis", "20151027_MERS", "20151028_InfluenzaA", "20151028_Smallpox")
desired_roots_list <- c("20160505_Ebola", "20160506_InfluenzaA", "20160506_MERS", "20160506_Smallpox", "20160507_HepatitisA", "20160507_Pertussis", "20160507_SARS")
# desired_roots_list <- c("20151022_SARS", "20151028_InfluenzaA")

# load the first one to get a list of the headers
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_roots_list[1], "/", desired_roots[1], "_PRCC.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_roots_list[1], "_PRCC.RData", sep=""))

desired_roots_list <- c("20160505_Ebola", "20160506_InfluenzaA", "20160506_MERS", "20160506_Smallpox", "20160507_HepatitisA", "20160507_Pertussis", "20160507_SARS")

names <- c(names(data.prcc), "disease")
length(names)
df.prcc.temporary <- data.frame(matrix(rep(NA, length(names)), nrow=1))
names(df.prcc.temporary) <- names

names <- c(names(output), "disease")
length(names)
df.prcc.output.temporary <- data.frame(matrix(rep(NA, length(names)), nrow=1))
names(df.prcc.output.temporary) <- names

for (i.temporary in 1:length(desired_roots_list)){
  root_i <- desired_roots_list[i.temporary]
  
  data.prcc.temp <- read.csv(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root_i, "_dataprcc.csv", sep=""))
  output.temp <- read.csv(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root_i, "_outputprcc.csv", sep=""))
  
  data.prcc.temp$disease <- substr(root_i, 10, nchar(root_i))
  output.temp$disease <- substr(root_i, 10, nchar(root_i))
  
  cat("\n",length(names(df.prcc.temporary)))
  cat("\n",length(names(df.prcc.output.temporary)))
  
  df.prcc.temporary <- rbind(df.prcc.temporary, data.prcc.temp)
  df.prcc.output.temporary <- rbind(df.prcc.output.temporary, output.temp)
}

df.prcc.temporary <- df.prcc.temporary[is.na(df.prcc.temporary$R_0)==0,]
df.prcc.output.temporary <- df.prcc.output.temporary[is.na(df.prcc.output.temporary$coef)==0,]

#### Save Workspace with all diseases ####
(date <- format(Sys.time(), "%Y%m%d"))
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PRCC.RData", sep=""))

#### Calculate PRCC for many diseases together ####
dep_var <- c("R_0", "R_hsb","R_s", "R_q","Impact_SM", "Impact_Q", "Impact_HSB", "Abs_Benefit","Abs_Benefit_per_Qday","Abs_Benefit_per_Qday_rp", "NNQ", "NQD", "Rel_Benefit","obs_to_iso_q","ks", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon", "riskprofile", "R_0_mean","dispersion","T_inc_stretch","pi_t_triangle_center","T_lat_offset","d_inf")
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

df.prcc.temporary <- df.prcc.temporary[(df.prcc.temporary$obs_to_iso_q %in% c(0, NA))==0,]

output.all <- prcc_fcn(input_data = df.prcc.temporary, dep_var = dep_var, indep_var = indep_var, 
                   nboot = 100, package = "sensitivity", standardize = TRUE)

# Add output.all to the df.prcc.output file so an all-diseases bar can be added to the horizontal grouped bar chart
output.all$disease <- "all"
df.prcc.output.temporary <- rbind(df.prcc.output.temporary, output.all)

#### Sensitivity check - control for disease as a parameter ####
# Repeat, controlling for disease? Conclusion - it doesn't make a difference really
df.prcc.temporary.b <- df.prcc.temporary
df.prcc.temporary.b$disease <- as.numeric(as.factor(df.prcc.temporary$disease))
dep_var <- c("R_0", "R_hsb","R_s", "R_q", "Abs_Benefit","Abs_Benefit_per_Qday", "NNQ", "NQD", "Rel_Benefit","obs_to_iso_q","ks", "Rel_Benefit_per_Qday", "Rel_Benefit_per_Qday_rp")
indep_var <- c("gamma","prob_CT","CT_delay", "epsilon", "riskprofile", "R_0_mean","dispersion","T_inc_stretch","pi_t_triangle_center","T_lat_offset","d_inf", "disease")
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")
output.all.b <- prcc_fcn(input_data = df.prcc.temporary.b, dep_var = dep_var, indep_var = indep_var, 
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

cat(date)
pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_Plot_prcc_2.pdf", sep=""))
plot(plot_prcc_2)
dev.off()

plot_prcc_2.b <- ggplot(output.all.b, aes(x = parameter, y= coef)) +
  facet_grid(output ~ .) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  geom_errorbar(data = output.all, aes(ymin = CImin, ymax = CImax), width = 0.1) +
  ggtitle("All diseases")
plot_prcc_2.b

# Note that obs_to_iso_q is NOT monotonic for the pooled estimate for T_lat_offset, pi_t_triangle_center, d_inf, and gamma.
# for (covariate in indep_var){
#   panel_plot_fcn(data = df.prcc.temporary, covariate = covariate, outputs = c("obs_to_iso_q", "Rel_Benefit", "Rel_Benefit_per_Qday"))
#   cat ("Press [enter] to continue")
#   line <- readline()
# }

#### Plot all diseases, horizontal bar chart ####
plot_prcc_3 <- ggplot(output.all, aes(x = parameter, y = coef)) +
  facet_grid(.~output) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  ggtitle("All diseases") +
  coord_flip()
plot_prcc_3

#### Plot each disease, horizontal bar chart, all outputs ####
df.prcc.output.temporary.2 <- df.prcc.output.temporary[df.prcc.output.temporary$output %in% c("Impact_HSB", "Impact_SM", "Impact_Q", "Abs_Benefit", "Rel_Benefit", "Abs_Benefit_per_Qday_rp", "Rel_Benefit_per_Qday_rp", "obs_to_iso_q"),]
df.prcc.output.temporary.2$parameter <- factor(df.prcc.output.temporary.2$parameter, 
                                               levels = rev(c("d_inf", "T_lat_offset", "pi_t_triangle_center", "CT_delay", "dispersion", "riskprofile", "T_inc_stretch", "R_0_mean", "epsilon", "gamma", "prob_CT")), 
                                               ordered = TRUE)
df.prcc.output.temporary.2$output <- factor(df.prcc.output.temporary.2$output, levels = c("Impact_HSB", "Impact_SM", "Impact_Q", "Abs_Benefit", "Rel_Benefit", "Abs_Benefit_per_Qday_rp", "Rel_Benefit_per_Qday_rp", "obs_to_iso_q"), ordered = TRUE)
levels(df.prcc.output.temporary.2$output) <-  c("R[0]-R[HSB]",
                                                "R[0]-R[S]",
                                                "R[0]-R[Q]",
                                                "R[S]-R[Q]",
                                                "frac(R[S]-R[Q],R[S])",
                                                "(frac(R[S]-R[Q],d[Q]))",
                                                "(frac(R[S]-R[Q],R[S]))/d[Q]",
                                                "d[Q]")
df.prcc.output.temporary.2$disease <- factor(df.prcc.output.temporary.2$disease, levels = rev(c("Pertussis","Smallpox", "SARS", "HepatitisA", "InfluenzaA", "Ebola", "MERS", "all")), ordered = TRUE, labels = rev(c("Pertussis","Smallpox", "SARS", "Hepatitis A", "Influenza A", "Ebola", "MERS", "All Diseases")))
df.prcc.output.temporary.2 <- df.prcc.output.temporary.2[is.na(df.prcc.output.temporary.2$output)==0,]
df.prcc.output.temporary.2 <- df.prcc.output.temporary.2[df.prcc.output.temporary.2$output != "d[Q]",]

scale_colour_brewer(type="qual", palette=6)
my.cols <- brewer.pal(n = 7, name = "Set1")
# my.cols <- c(my.cols[c(3, 7, 4, 5, 1, 6, 2)])
my.cols <- c(my.cols[c(5, 3, 4, 7, 6, 2, 1)])

# Note that obs_to_iso_q is NOT monotonic for the pooled estimate for T_lat_offset, pi_t_triangle_center, d_inf, and gamma.
plot_prcc_6 <- ggplot() +
  geom_bar(data = df.prcc.output.temporary.2[df.prcc.output.temporary.2$disease == "All Diseases" &
                                               df.prcc.output.temporary.2$output != "d[Q]" &
                                               df.prcc.output.temporary.2$output != "(frac(R[S]-R[Q],R[S]))/d[Q]" &
                                               df.prcc.output.temporary.2$output != "(frac(R[S]-R[Q],d[Q]))",],
           aes(x=parameter, y=coef), stat = "identity", alpha = 0.5, fill = "grey",  color = "darkgrey") +
  facet_grid(.~output, labeller=label_parsed) +
  geom_bar(data = df.prcc.output.temporary.2[df.prcc.output.temporary.2$disease != "All Diseases",], aes(x = parameter, y = coef, fill = disease), position = "dodge", stat="identity", width = .5, alpha = 0.5) +
#   ggtitle("Partial Rank Correlation Coefficient") +
  coord_flip() +
  # scale_fill_manual(values = rep(c("grey", "black"), 4)) +
  # scale_fill_brewer(type = "div", palette = 6) +
  scale_fill_manual(values = my.cols, breaks = c("Pertussis","Smallpox", "SARS", "Hepatitis A", "Influenza A", "Ebola", "MERS", "All Diseases"), name = "Disease") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  # ylab("Partial Rank Correlation Coefficient") +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1", "", "0", "", "+1")) +
  geom_hline(yintercept=0, color="darkgrey", size=0.7) +
  guides(color=FALSE) +
  # guides(fill=FALSE) +
  scale_x_discrete(name = element_blank(),
                   labels=c(expression(paste("Fraction of\nContacts Traced ", (P[CT]), sep="")),
                            expression(paste("Isolation\nEffectiveness ", (gamma), sep="")),
                            expression(paste("Delay from Symptom\nOnset to Isolation ", (D[SM]), sep="")),
                            expression(paste("Basic Reproductive\nNumber ", (R[0]), sep="")),
                            expression(paste("Incubation\nPeriod ", (T[INC]), sep="")),
                            expression(paste("Fraction of Traced Contacts\nwho are Infected ", (P[INF]), sep="")),
                            expression(paste("Dispersion\nFactor ", (kappa), sep="")),
                            expression(paste("Delay in Tracing\na Contact ", (D[CT]), sep="")),
                            expression(paste("Time of Peak\nInfectiousness ", (tau[beta]), sep="")),
                            expression(paste("Latent Period\nOffset ", (T[OFFSET]), sep="")),
                            expression(paste("Duration of\nInfectiousness ", (d[INF]), sep=""))))
plot_prcc_6

#### Mini-plots of select variables to compliment PRCC barplots ####
plot(df.prcc.temporary[df.prcc.temporary$disease == "InfluenzaA","prob_CT"], df.prcc.temporary[df.prcc.temporary$disease == "InfluenzaA","Abs_Benefit"])

ggplot(df.prcc.temporary[df.prcc.temporary$disease == "InfluenzaA",]) +
  geom_boxplot(aes(x = prob_CT, group = cut_width(prob_CT, 0.1), y = Abs_Benefit))

desired_root <- "20151028_InfluenzaA" # Paste the desired root here "YYYYMMDD_DISEASE"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

# Vary Prob_CT
prob_CT_range <- seq(0,1,0.1)
prob_CT_data <- data.frame(prob_CT = NA, R_s = NA, R_q = NA, Abs_Benefit = NA, Abs_Benefit_per_Qday = NA)

for (prob_CT in prob_CT_range){

  out <- intervention_effect_fcn(background_intervention = "u",
                                        prob_CT = prob_CT, gamma = NA, parms_epsilon = NA, parms_CT_delay = NA,
                                        resource_level = "high",
                                        n_pop = 2000, num_generations = 5, times = 5,
                                        input_data = data,
                                        parms_T_lat, parms_d_inf, parms_pi_t, parms_R_0, dispersion, 
                                        parms_serial_interval,
                                        printing = FALSE)
  out$prob_CT <- prob_CT
  prob_CT_data <- rbind(prob_CT_data, out[,c("prob_CT", "R_s", "R_q","Abs_Benefit","Abs_Benefit_per_Qday")])
  
}

ggplot(prob_CT_data, aes(x = prob_CT)) +
  theme_classic() +
  geom_point(aes(y = Abs_Benefit), color = "darkred") +
  geom_point(aes(y = Abs_Benefit_per_Qday), color = "blue")

ggplot(prob_CT_data, aes(x = prob_CT)) + 
  theme_classic() +
  geom_boxplot(aes(group = cut_interval(prob_CT, 10), y = Abs_Benefit), color = "darkred")

ggplot(prob_CT_data, aes(x = prob_CT)) + 
  theme_classic() +
  geom_smooth(aes(y = Abs_Benefit), color = "darkred") +
  geom_smooth(aes(y = Abs_Benefit_per_Qday), color = "blue")

ggplot(prob_CT_data, aes(x = prob_CT)) + 
  theme_classic() +
  # geom_point(aes(y = R_s), color = "darkred", alpha = 0.1) +
  geom_smooth(aes(y = R_s), color = "darkred") +
  # geom_point(aes(y = R_q), color = "blue", alpha = 0.1) +
  geom_smooth(aes(y = R_q), color = "blue")
  
# Vary gamma
gamma_range <- seq(0,1,0.1)
gamma_data <- data.frame(gamma = NA, R_s = NA, R_q = NA, Abs_Benefit = NA, Abs_Benefit_per_Qday = NA)

for (gamma in gamma_range){
  
  out <- intervention_effect_fcn(background_intervention = "u",
                                 prob_CT = NA, gamma = gamma, parms_epsilon = NA, parms_CT_delay = NA,
                                 resource_level = "high",
                                 n_pop = 5000, num_generations = 5, times = 5,
                                 input_data = data,
                                 parms_T_lat, parms_d_inf, parms_pi_t, parms_R_0, dispersion, 
                                 parms_serial_interval,
                                 printing = FALSE)
  out$gamma <- gamma
  gamma_data <- rbind(gamma_data, out[,c("gamma", "R_s", "R_q","Abs_Benefit","Abs_Benefit_per_Qday")])
  
}

ggplot(gamma_data, aes(x = gamma)) +
  theme_classic() +
  geom_point(aes(y = Abs_Benefit), color = "darkred") +
  geom_point(aes(y = Abs_Benefit_per_Qday), color = "blue")

ggplot(gamma_data, aes(x = gamma)) + 
  theme_classic() +
  geom_boxplot(aes(group = cut_interval(gamma, 10), y = Abs_Benefit), color = "darkred")

ggplot(gamma_data, aes(x = gamma)) + 
  theme_classic() +
  geom_smooth(aes(y = Abs_Benefit), color = "darkred") +
  geom_smooth(aes(y = Abs_Benefit_per_Qday), color = "blue")

ggplot(gamma_data, aes(x = gamma)) + 
  theme_classic() +
  geom_smooth(aes(y = R_s), color = "darkred") +
  geom_smooth(aes(y = R_q), color = "blue")

#### circleFun for plotting on ggplot ####
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circle_0.25 <- circleFun(diameter = 0.5, npoints = 100)
circle_0.5 <- circleFun(diameter = 1, npoints = 100)

#### Make a "Target" Plot for R_SM and R_Q ####
df.prcc.output.subset_R <- df.prcc.output.subset[is.element(df.prcc.output.subset$output, c("Symptom Monitoring R", "Quarantine R")), ]
df.prcc.output.subset_R.melt <- melt(df.prcc.output.subset_R, id.vars = c("parameter", "disease", "output"))

df.prcc.output.subset_R.target <- data.frame(cbind(disease = as.character(rep(unique(df.prcc.output.subset_R.melt$disease), times = length(unique(df.prcc.output.subset_R.melt$parameter)))), parameter = as.character(rep(unique(df.prcc.output.subset_R.melt$parameter), each = length(unique(df.prcc.output.subset_R.melt$disease))))))
df.prcc.output.subset_R.target$RSM <- NA
df.prcc.output.subset_R.target$RQ <- NA

for (i in 1:nrow(df.prcc.output.subset_R.target)){
  df.prcc.output.subset_R.target[i, "RSM"] <- df.prcc.output.subset_R[as.character(df.prcc.output.subset_R$disease) == as.character(df.prcc.output.subset_R.target[i,"disease"]) & as.character(df.prcc.output.subset_R$parameter) == as.character(df.prcc.output.subset_R.target[i,"parameter"]) & as.character(df.prcc.output.subset_R$output) == "Symptom Monitoring R", "coef"]
  
  df.prcc.output.subset_R.target[i, "RQ"] <- df.prcc.output.subset_R[as.character(df.prcc.output.subset_R$disease) == as.character(df.prcc.output.subset_R.target[i,"disease"]) & as.character(df.prcc.output.subset_R$parameter) == as.character(df.prcc.output.subset_R.target[i,"parameter"]) & as.character(df.prcc.output.subset_R$output) == "Quarantine R", "coef"]
}


ggplot() +
  geom_hline(yintercept=0, color = "darkgrey") +
  geom_vline(xintercept=0, color = "darkgrey") +
  geom_path(data=circle_0.25, aes(x,y), color = "lightgrey") +
  geom_path(data=circle_0.5, aes(x,y), color = "lightgrey") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_text(data=df.prcc.output.subset_R.target[df.prcc.output.subset_R.target$disease == "All Diseases",], aes(x = RSM, y = RQ, label=parameter, color = parameter, size = ((abs(RSM)^.1)+abs(RQ)^.1)), angle=0) +
  geom_point(data=df.prcc.output.subset_R.target[df.prcc.output.subset_R.target$disease == "All Diseases",], aes(x = RSM, y = RQ, color = parameter), shape = 5) +
  xlim(-1,1) +
  ylim(-1,1) +
  guides(size = FALSE, color = FALSE) +
  coord_fixed() +
  scale_size(range=c(2,6)) +
  xlab(expression(paste(R[S]))) +
  ylab(expression(paste(R[Q])))

#### Make a "Target" Plot for Rel_Benefit and Rel_Benefit_per_Qday ####
df.prcc.output.subset_RelBen <- df.prcc.output.subset[is.element(df.prcc.output.subset$output, c("Relative Difference", "Relative Difference\nper Quarantine Day")), ]
df.prcc.output.subset_RelBen.melt <- melt(df.prcc.output.subset_RelBen, id.vars = c("parameter", "disease", "output"))

df.prcc.output.subset_RelBen.target <- data.frame(cbind(disease = as.character(rep(unique(df.prcc.output.subset_RelBen.melt$disease), times = length(unique(df.prcc.output.subset_RelBen.melt$parameter)))), parameter = as.character(rep(unique(df.prcc.output.subset_RelBen.melt$parameter), each = length(unique(df.prcc.output.subset_RelBen.melt$disease))))))
df.prcc.output.subset_RelBen.target$Rel_Benefit <- NA
df.prcc.output.subset_RelBen.target$Rel_Benefit_per_Qday <- NA

for (i in 1:nrow(df.prcc.output.subset_RelBen.target)){
  df.prcc.output.subset_RelBen.target[i, "Rel_Benefit"] <- df.prcc.output.subset_RelBen[as.character(df.prcc.output.subset_RelBen$disease) == as.character(df.prcc.output.subset_RelBen.target[i,"disease"]) & as.character(df.prcc.output.subset_RelBen$parameter) == as.character(df.prcc.output.subset_RelBen.target[i,"parameter"]) & as.character(df.prcc.output.subset_RelBen$output) == "Relative Difference", "coef"]
  
  df.prcc.output.subset_RelBen.target[i, "Rel_Benefit_per_Qday"] <- df.prcc.output.subset_RelBen[as.character(df.prcc.output.subset_RelBen$disease) == as.character(df.prcc.output.subset_RelBen.target[i,"disease"]) & as.character(df.prcc.output.subset_RelBen$parameter) == as.character(df.prcc.output.subset_RelBen.target[i,"parameter"]) & as.character(df.prcc.output.subset_RelBen$output) == "Relative Difference\nper Quarantine Day", "coef"]
}

ggplot() +
  geom_hline(yintercept=0, color = "darkgrey") +
  geom_vline(xintercept=0, color = "darkgrey") +
  geom_path(data=circle_0.25, aes(x,y), color = "lightgrey") +
  geom_path(data=circle_0.5, aes(x,y), color = "lightgrey") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_text(data=df.prcc.output.subset_RelBen.target[df.prcc.output.subset_RelBen.target$disease == "All Diseases",], aes(x = Rel_Benefit, y = Rel_Benefit_per_Qday, label=parameter, color = parameter, size = ((abs(Rel_Benefit)^.5)+abs(Rel_Benefit_per_Qday)^.5)), angle=0) +
  geom_point(data=df.prcc.output.subset_RelBen.target[df.prcc.output.subset_RelBen.target$disease == "All Diseases",], aes(x = Rel_Benefit, y = Rel_Benefit_per_Qday, color = parameter), shape = 5) +
  xlim(-1,1) +
  ylim(-1,1) +
  guides(size = FALSE, color = FALSE) +
  coord_fixed() +
  scale_size(range=c(2,6)) +
  xlab(expression(paste("Relative Effectiveness of Quarantine ", frac(R[S]-R[Q],R[S]), sep=""))) +
  ylab(expression(paste("Relative Effeciency of Quarantine ", (frac(R[S]-R[Q],R[S]))/d[Q], sep="")))
  

#### Make a horizontal bar plot grouped by parameter for absolute and relative comparative effectveness ####
# Pull df.prcc.output.subset from "Plot each disease, horizontal bar chart, subset of outputs" section above
df.prcc.output.subset.plot1 <- df.prcc.output.subset[df.prcc.output.subset$disease == "All Diseases" & df.prcc.output.subset$output %in% c("Relative Difference", "Absolute Difference"),]
df.prcc.output.subset.plot1 <- df.prcc.output.subset.plot1[order(-df.prcc.output.subset.plot1$coef),]

# df.prcc.output.subset.plot1 <- df.prcc.output.subset.plot1[(df.prcc.output.subset.plot1$parameter %in% c("R_0 Spread"))==0,]

df.prcc.output.subset.plot1$parameter <- factor(df.prcc.output.subset.plot1$parameter, levels = (unique(as.character(df.prcc.output.subset.plot1$parameter))), ordered = TRUE)
levels(df.prcc.output.subset.plot1$parameter )
df.prcc.output.subset.plot1$output <- factor(df.prcc.output.subset.plot1$output, levels = rev(c("Absolute Difference", "Relative Difference"))) # Reverse it so absolute difference is on top after coord flip
levels(df.prcc.output.subset.plot1$output)

plot_1 <- ggplot(data = df.prcc.output.subset.plot1, aes(y = coef, fill = output)) +
  # theme_bw() +
  theme_classic() + 
  geom_bar(aes(x = parameter), stat = "identity", position = "dodge") +
  geom_errorbar(aes(x = parameter, ymin = CImin, ymax = CImax, color = output), position = position_dodge(width=1), width=0) +
  scale_fill_brewer(type = "qual", palette = 7, name = "Comparative Metric", 
                     breaks = c("Absolute Difference", "Relative Difference"),
                     labels = c(expression(paste("Absolute Difference  ", (R[S]-R[Q]), sep="")),
                                expression(paste("Relative Difference  ", (frac(R[S]-R[Q],R[S])), sep="")))) +
  scale_color_brewer(type = "qual", palette = 7, name = "Comparative Metric", 
                    breaks = c("Absolute Difference", "Relative Difference"),
                    labels = c(expression(paste("Absolute Difference  ", (R[S]-R[Q]), sep="")),
                               expression(paste("Relative Difference  ", (frac(R[S]-R[Q],R[S])), sep="")))) +
  theme(axis.text.y = element_text(hjust = 0)) +
  ylab("Partial Rank Correlation Coefficient") +
  theme(axis.ticks.y = element_blank()) +
  scale_x_discrete(name = element_blank(),
                   labels=c(expression(paste("Fraction of\nContacts Traced ", (P[CT]), sep="")),
                            expression(paste("Isolation\nEffectiveness ", (gamma), sep="")),
                            expression(paste("Delay from Symptom\nOnset to Isolation ", (D[SM]), sep="")),
                            expression(paste("Basic Reproductive\nNumber ", (R[0]), sep="")),
                            expression(paste("Incubation\nPeriod ", (T[INC]), sep="")),
                            expression(paste("Fraction of Traced Contacts\nwho are Infected ", (P[INF]), sep="")),
                            expression(paste("Dispersion\nFactor ", (kappa), sep="")),
                            expression(paste("Delay in Tracing\na Contact ", (D[CT]), sep="")),
                            expression(paste("Time of Peak\nInfectiousness ", (tau[beta]), sep="")),
                            expression(paste("Latent Period\nOffset ", (T[OFFSET]), sep="")),
                            expression(paste("Duration of\nInfectiousness ", (d[INF]), sep="")))) +
  theme(legend.position = c(0.75, 0.75)) +
  theme(text = element_text(size=18)) +
  coord_flip()
plot_1

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PRCC_Plot1_10by7.pdf", sep=""), width = 10, height = 7)
plot(plot_1)
dev.off()

# Comparative cost-effectiveness of Q and SM
df.prcc.output.subset.plot2 <- df.prcc.output.subset[df.prcc.output.subset$disease == "All Diseases" & df.prcc.output.subset$output %in% c("Relative Difference\nper Quarantine Day", "Absolute Difference\nper Quarantine Day"),]
df.prcc.output.subset.plot2 <- df.prcc.output.subset.plot2[order(-df.prcc.output.subset.plot2$coef),]

# df.prcc.output.subset.plot2 <- df.prcc.output.subset.plot2[(df.prcc.output.subset.plot2$parameter %in% c("R_0 Spread"))==0,]

df.prcc.output.subset.plot2$parameter <- factor(df.prcc.output.subset.plot2$parameter, levels = (unique(as.character(df.prcc.output.subset.plot1$parameter))), ordered = TRUE)
levels(df.prcc.output.subset.plot2$parameter) == levels(df.prcc.output.subset.plot1$parameter)
df.prcc.output.subset.plot2$output <- factor(df.prcc.output.subset.plot2$output, levels = rev(c("Absolute Difference\nper Quarantine Day", "Relative Difference\nper Quarantine Day"))) # Reverse it so absolute difference is on top after coord flip
levels(df.prcc.output.subset.plot2$output)

ggplot(data = df.prcc.output.subset.plot2, aes(y = coef, fill = output)) +
  theme_classic() + 
  geom_bar(aes(x = parameter), stat = "identity", position = "dodge") +
  geom_errorbar(aes(x = parameter, ymin = CImin, ymax = CImax, color = output), position = position_dodge(width=1), width=0) +
  scale_fill_brewer(type = "qual", palette = 7, name = "Comparative Metric", 
                    breaks = c("Absolute Difference\nper Quarantine Day", "Relative Difference\nper Quarantine Day"),
                    labels = c(expression(paste("Absolute Difference\nper Quarantine Day  ", (R[S]-R[Q]), sep="")),
                               expression(paste("Relative Difference\nper Quarantine Day  ", (frac(R[S]-R[Q],d[Q])), sep="")))) +
  scale_color_brewer(type = "qual", palette = 7, name = "Comparative Metric", 
                     breaks = c("Absolute Difference\nper Quarantine Day", "Relative Difference\nper Quarantine Day"),
                     labels = c(expression(paste("Absolute Difference\nper Quarantine Day  ", (R[S]-R[Q]), sep="")),
                                expression(paste("Relative Difference\nper Quarantine Day  ", (frac(R[S]-R[Q],d[Q])), sep="")))) +
  theme(axis.text.y = element_text(hjust = 0)) +
  ylab("Partial Rank Correlation Coefficient") +
  theme(axis.ticks.y = element_blank()) +
  scale_x_discrete(name = element_blank(),
                   labels=c(expression(paste("Fraction of\nContacts Traced ", (P[CT]), sep="")),
                            expression(paste("Isolation\nEffectiveness ", (gamma), sep="")),
                            expression(paste("Delay from Symptom\nOnset to Isolation ", (D[SM]), sep="")),
                            expression(paste("Basic Reproductive\nNumber ", (R[0]), sep="")),
                            expression(paste("Incubation\nPeriod ", (T[INC]), sep="")),
                            expression(paste("Fraction of Traced Contacts\nwho are Infected ", (P[INF]), sep="")),
                            expression(paste("Dispersion\nFactor ", (kappa), sep="")),
                            expression(paste("Delay in Tracing\na Contact ", (D[CT]), sep="")),
                            expression(paste("Time of Peak\nInfectiousness ", (tau[beta]), sep="")),
                            expression(paste("Latent Period\nOffset ", (T[OFFSET]), sep="")),
                            expression(paste("Duration of\nInfectiousness ", (d[INF]), sep="")))) +
  theme(legend.position = c(0.25, 0.15)) +
  coord_flip()

# Impact of Q and SM
df.prcc.output.subset.plot3 <- df.prcc.output.subset[df.prcc.output.subset$disease == "All Diseases" & df.prcc.output.subset$output %in% c("Symptom Monitoring R", "Quarantine R"),]
df.prcc.output.subset.plot3 <- df.prcc.output.subset.plot3[order(-df.prcc.output.subset.plot3$coef),]

# df.prcc.output.subset.plot3 <- df.prcc.output.subset.plot3[(df.prcc.output.subset.plot3$parameter %in% c("R_0 Spread"))==0,]

df.prcc.output.subset.plot3$parameter <- factor(df.prcc.output.subset.plot3$parameter, levels = (unique(as.character(df.prcc.output.subset.plot1$parameter))), ordered = TRUE)
levels(df.prcc.output.subset.plot3$parameter) == levels(df.prcc.output.subset.plot1$parameter)
df.prcc.output.subset.plot3$output <- factor(df.prcc.output.subset.plot3$output, levels = rev(c("Symptom Monitoring R", "Quarantine R"))) # Reverse it so absolute difference is on top after coord flip
levels(df.prcc.output.subset.plot3$output)

ggplot(data = df.prcc.output.subset.plot3, aes(y = coef, fill = output)) +
  theme_classic() + 
  geom_bar(aes(x = parameter), stat = "identity", position = "dodge") +
  geom_errorbar(aes(x = parameter, ymin = CImin, ymax = CImax, color = output), position = position_dodge(width=1), width=0) +
  scale_fill_brewer(type = "qual", palette = 7, name = "Metric", 
                    breaks = c("Symptom Monitoring R", "Quarantine R"),
                    labels = c(expression(paste("Symptom Monitoring  ", (R[S]), sep="")),
                               expression(paste("Quarantine  ", (R[Q]), sep="")))) +
  scale_color_brewer(type = "qual", palette = 7, name = "Metric", 
                     breaks = c("Symptom Monitoring R", "Quarantine R"),
                     labels = c(expression(paste("Symptom Monitoring  ", (R[S]), sep="")),
                                expression(paste("Quarantine  ", (R[Q]), sep="")))) +
  theme(axis.text.y = element_text(hjust = 0)) +
  ylab("Partial Rank Correlation Coefficient") +
  theme(axis.ticks.y = element_blank()) +
  scale_x_discrete(name = element_blank(),
                   labels=c(expression(paste("Fraction of\nContacts Traced ", (P[CT]), sep="")),
                            expression(paste("Isolation\nEffectiveness ", (gamma), sep="")),
                            expression(paste("Delay from Symptom\nOnset to Isolation ", (D[SM]), sep="")),
                            expression(paste("Basic Reproductive\nNumber ", (R[0]), sep="")),
                            expression(paste("Incubation\nPeriod ", (T[INC]), sep="")),
                            expression(paste("Dispersion\nFactor ", (kappa), sep="")),
                            expression(paste("Fraction of Traced Contacts\nwho are Infected ", (P[INF]), sep="")),
                            expression(paste("Delay in Tracing\na Contact ", (D[CT]), sep="")),
                            expression(paste("Time of Peak\nInfectiousness ", (tau[beta]), sep="")),
                            expression(paste("Latent Period\nOffset ", (T[OFFSET]), sep="")),
                            expression(paste("Duration of\nInfectiousness ", (d[INF]), sep="")))) +
  theme(legend.position = c(0.75, 0.15)) +
  coord_flip()

#### Make a vertical bar plot grouped by parameter for impact of Q and SM ####
# Pull df.prcc.output.subset from "Plot each disease, horizontal bar chart, subset of outputs" section above
df.prcc.output.subset.plot2 <- df.prcc.output.subset[df.prcc.output.subset$disease == "All Diseases" & df.prcc.output.subset$output %in% c("Quarantine R", "Symptom Monitoring R"),]
df.prcc.output.subset.plot2 <- df.prcc.output.subset.plot2[order(-df.prcc.output.subset.plot2$coef),]
df.prcc.output.subset.plot2$parameter <- factor(df.prcc.output.subset.plot2$parameter, levels = unique(as.character(df.prcc.output.subset.plot2$parameter))[c(1,2,3,4,5,8,9,6,10,7,11,12)], ordered = TRUE)

ggplot(data = df.prcc.output.subset.plot2, aes(y = coef, fill = output)) +
  theme_bw() + 
  geom_bar(aes(x = parameter), stat = "identity", position = "dodge") +
  geom_errorbar(aes(x = parameter, ymin = CImin, ymax = CImax, color = output), position = position_dodge(width=1), width=0) +
  scale_fill_brewer(type = "qual", palette = 3) +
  scale_color_brewer(type = "qual", palette = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#### Save Workspace ####
date <- format(Sys.time(), "%Y%m%d")
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", date, "_PRCC.RData", sep=""))

