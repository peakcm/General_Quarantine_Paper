#### Header ####
# Confirm behavior of SMC algorithm

#### Load Libraries ####
library(MASS)
library(lhs)
library(reshape2)
library(ggplot2)

#### Source Functions ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

#### Disease: SimvirusA ####
# Essentially Ebola in an ideal case
# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "SimvirusA"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5931, 0.1697)   # This will actually be generated from the other inputs
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, .182, 999, "independent", "independent", 1)
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target", "T_inc_stretch")

parms_R_0 = list("uniform", 1.2, 1.2, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")   # Set offset
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")    # link to d_inf. Doesn't matter
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

#### Disease: SimvirusB ####
# Same as SimvirusA, but now the d_inf follows a triangular distribution and trying to fit it using a uniform distribution

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "SimvirusB"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5931, 0.1697)   # This will actually be generated from the other inputs
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, .182, 999, "independent", "independent", 1)
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target", "T_inc_stretch")

parms_R_0 = list("uniform", 1.2, 1.2, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.20)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")   # Set offset
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("triangle", 1, 15, 5, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")    # link to d_inf. Doesn't matter
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

#### Disease: SimvirusC ####
# Same as SimvirusA, infectiousness follows symptoms by 3 days

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "SimvirusC"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5931, 0.1697)   # This will actually be generated from the other inputs
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, .182, 999, "independent", "independent", 1)
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target", "T_inc_stretch")

parms_R_0 = list("uniform", 1.2, 1.2, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 3, "T_inc")   # Set offset
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")    # link to d_inf. Doesn't matter
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

#### Disease: SimvirusD ####
# Same as SimvirusA, infectiousness precedes symptoms by 3 days

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "SimvirusD"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5931, 0.1697)   # This will actually be generated from the other inputs
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, .182, 999, "independent", "independent", 1)
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target", "T_inc_stretch")

parms_R_0 = list("uniform", 1.2, 1.2, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, -3, "T_inc")   # Set offset
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")    # link to d_inf. Doesn't matter
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

#### Disease: SimvirusE ####
# Like Pertussis in an ideal case

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "SimvirusE"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.45585778, .11071164) # Approximation from te Beest reported quantiles
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("normal", 7, 1.53, 999, "independent", "independent", 1) # (10-7)/1.96 = 1.53
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target", "T_inc_stretch")

parms_R_0 = list("uniform", 1.2, 1.2, 999, "independent", "independent") 
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.45)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, -2, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 70, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

#### Initialize intervention parameters ####
parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

prob_CT <- 1

gamma <- 0.9

background_intervention <- "u"

#### Simulate chains of infections with the simulated disease ####
dispersion = 1
n_pop = 5000
Pop_alpha <- Create_Pop(n_pop, 
                        parms_T_inc, 
                        parms_T_lat, 
                        parms_d_inf, 
                        parms_d_symp, 
                        parms_R_0, 
                        parms_epsilon, 
                        generation = 1,
                        background_intervention,
                        parms_CT_delay,
                        gamma)
Pop_alpha <- observe_and_isolate_fcn(Pop_alpha, intervention = background_intervention)
children_list <- children_list_fcn(Pop_alpha, parms_pi_t$distribution, parms_pi_t$triangle_center, gamma, intervention = background_intervention, background_intervention, dispersion = dispersion)
Num_Infected <- unlist(lapply(children_list, sum))
cat('Generation 1 : n=', nrow(Pop_alpha), '. Effective Reproductive Number:', mean(Num_Infected), '. Number of infections:', sum(Num_Infected))

intervention = "u"
Pop_beta <- next_generation_fcn(Pop = Pop_alpha,
                                children_list = children_list,
                                parms_T_inc = parms_T_inc,
                                parms_T_lat = parms_T_lat,
                                parms_d_inf = parms_d_inf,
                                parms_d_symp = parms_d_symp,
                                parms_R_0 = parms_R_0,
                                parms_epsilon = parms_epsilon,
                                generation = 2,
                                parms_CT_delay = parms_CT_delay,
                                prob_CT = prob_CT,
                                gamma=gamma,
                                n_pop = n_pop)
Pop_beta <- observe_and_isolate_fcn(Pop_beta, intervention = intervention)
children_list <- children_list_fcn(Pop_beta, pi_t_distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma, intervention = intervention, background_intervention, dispersion = dispersion)
Num_Infected <- unlist(lapply(children_list, sum))
cat('Generation 2 : n=', nrow(Pop_beta), '. Effective Reproductive Number:', mean(Num_Infected), '. Number of infections:', sum(Num_Infected))

Pop_gamma <- next_generation_fcn(Pop_beta,
                                 children_list,
                                 parms_T_inc,
                                 parms_T_lat,
                                 parms_d_inf,
                                 parms_d_symp,
                                 parms_R_0,
                                 parms_epsilon,
                                 generation = 3,
                                 prob_CT = prob_CT,
                                 parms_CT_delay,
                                 gamma = gamma,
                                 n_pop)
Pop_gamma <- observe_and_isolate_fcn(Pop_gamma, intervention = intervention)
children_list <- children_list_fcn(Pop_gamma, parms_pi_t$distribution, parms_pi_t$triangle_center, gamma, intervention = intervention, background_intervention, dispersion = dispersion)
Num_Infected <- unlist(lapply(children_list, sum))
cat('Generation 3 : n=', nrow(Pop_gamma), '. Effective Reproductive Number:', mean(Num_Infected), '. Number of infections:', sum(Num_Infected))

#### Generate a serial interval distribution and fit for the simulated disease ####
sim_SI <- serial_interval_fcn(Pop_alpha, Pop_beta, parms_serial_interval, plot="False", values_returned = "SI")
(sim_SI_fit <- serial_interval_fcn(Pop_beta, Pop_gamma, parms_serial_interval, plot="True", values_returned = "fit"))

parms_serial_interval[c(2,3)] <- as.numeric(sim_SI_fit$estimate[c(1,2)])

#### Use SMC algorithm to fit pi_t, d_inf, and T_lat_offset ####
dir = "~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/"

# Ranges for particle filter
if ((disease %in% c("SimvirusE")) == 0){
  T_lat_offset.min <- -5
  T_lat_offset.max <- 5
  d_inf.min <- 1
  d_inf.max <- 40
  pi_t_triangle_center.min <- 0
  pi_t_triangle_center.max <- 1
} else {
  T_lat_offset.min <- -10
  T_lat_offset.max <- 4
  d_inf.min <- 5
  d_inf.max <- 120
  pi_t_triangle_center.min <- 0
  pi_t_triangle_center.max <- 1
}

# # Restricted
# T_lat_offset.min <- 0
# T_lat_offset.max <- 0
# d_inf.min <- 10
# d_inf.max <- 10
# pi_t_triangle_center.min <- 0.5
# pi_t_triangle_center.max <- 0.5

n_pop <- 500
times <- 500

outputs <- particle_filter_fcn(T_lat_offset.max = T_lat_offset.max,
                               T_lat_offset.min = T_lat_offset.min, 
                               d_inf.max = d_inf.max, 
                               d_inf.min = d_inf.min,
                               pi_t_triangle_center.max = pi_t_triangle_center.max,
                               pi_t_triangle_center.min = pi_t_triangle_center.min, 
                               parms_serial_interval = parms_serial_interval,
                               dir = dir,
                               ks_conv_criteria = 0.05, 
                               disease_name = disease, 
                               n_pop = n_pop, 
                               times = times,
                               num_generations = 2,
                               SMC_times = 15,
                               perturb_initial = 1/25,
                               perturb_final = 1/75)

head(outputs$data)
head(outputs$ks_conv_stat)

#### Explore bias and precision ####
df_bias_precision <- data.frame(matrix(NA, nrow = 3, ncol = 5))
names(df_bias_precision) <- c("variable", "true", "mean", "lower_ci", "upper_ci")
df_bias_precision$variable <- c("T_lat_offset", "d_inf", "pi_t_triangle_center")

if (parms_d_inf$dist == "uniform"){
  df_bias_precision[df_bias_precision$variable == "d_inf", "true"] <- (parms_d_inf$parm2 + 1)/2
}
if (parms_d_inf$dist == "triangle"){
  df_bias_precision[df_bias_precision$variable == "d_inf", "true"] <- (parms_d_inf$parm1 + parms_d_inf$parm2 + parms_d_inf$parm3)/3 # Mean of triangle distribution
}

df_bias_precision[df_bias_precision$variable == "T_lat_offset", "true"] <- parms_T_lat$anchor_value


if (parms_pi_t$distribution == "triangle"){
  df_bias_precision[df_bias_precision$variable == "pi_t_triangle_center", "true"] <- parms_pi_t$triangle_center
}

df_bias_precision$mean <- c(mean(outputs$data$T_lat_offset), mean((outputs$data$d_inf+1)/2), mean(outputs$data$pi_t_triangle_center))

df_gg <- melt(outputs$data)
tail(df_gg)

#### Plot parameter fits ####
ggplot(df_gg, aes(x = factor(variable), y = value)) +
  geom_violin(data = df_gg[df_gg$variable == "T_lat_offset",], aes(y = value)) +
  geom_boxplot(data = df_gg[df_gg$variable == "T_lat_offset",], aes(y = value), width = 0.1) +
  geom_violin(data = df_gg[df_gg$variable == "d_inf",], aes(y = (value+1)/2)) +
  geom_boxplot(data = df_gg[df_gg$variable == "d_inf",], aes(y = (value+1)/2), width = 0.1) +
  geom_violin(data = df_gg[df_gg$variable == "pi_t_triangle_center",], aes(y = value)) +
  geom_boxplot(data = df_gg[df_gg$variable == "pi_t_triangle_center",], aes(y = value), width = 0.1) +
  geom_point(data = df_bias_precision, aes(x = variable, y = true), color = "darkred", shape = 8, size = 3) +
  ggtitle("Comparison of Known and SMC-Fit parameters") +
  xlab("Parameter") +
  scale_x_discrete(labels = c("Mean Duration of Infectiousness (days)", "Peak Infectiousness\n(range 0 to 1)", "Mean offset between Latent Period\nand Incubation Period (days)"))

#### Calculate intervention performance for SMC-fit parameters ####

# Settings
n_pop = 500
num_generations <- 5
times <- 500
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.mr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.mr) <- names

background_intervention <- "u"

prob_CT <- 0.75

gamma <- 0.75

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

  # Set the parms_d_inf back to a uniform distribution from 1 to XXX
parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

data_fit <- intervention_effect_fcn(background_intervention = "u", resource_level = "medium", n_pop = n_pop, num_generations = num_generations, times = times, input_data = outputs$data, parms_serial_interval = parms_serial_interval, parms_R_0 = parms_R_0, parms_T_lat = parms_T_lat, parms_d_inf = parms_d_inf, parms_pi_t = parms_pi_t)

data_fit_melt <- melt(data_fit)

#### Calculate intervention performance for known input parameters ####
if (disease == "SimvirusA"){
  parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
  names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
}
if (disease == "SimvirusB"){
  parms_d_inf = list("triangle", 1, 15, 5, "independent", "independent")
  names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
}
if (disease == "SimvirusC"){
  parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
  names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
}
if (disease == "SimvirusD"){
  parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
  names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
}

n_pop = 500
times = 500

data_known <- intervention_effect_fcn(background_intervention = "u", resource_level = "medium", n_pop = n_pop, num_generations = num_generations, times = times, input_data = NA, parms_serial_interval = parms_serial_interval, parms_R_0 = parms_R_0, parms_T_lat = parms_T_lat, parms_d_inf = parms_d_inf, parms_pi_t = parms_pi_t)

data_known_melt <- melt(data_known)

#### Plots to compare impact of Q and SM ####
ggplot() +
  geom_violin(data = data_fit_melt[data_fit_melt$variable == "R_s",], aes(x = variable, y = value), alpha = 0.7) +
  geom_boxplot(data = data_fit_melt[data_fit_melt$variable == "R_s",], aes(x = variable, y = value), width = 0.1, alpha = 0.7) +
  geom_violin(data = data_known_melt[data_known_melt$variable == "R_s",], aes(x = variable, y = value), color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
  geom_boxplot(data = data_known_melt[data_known_melt$variable == "R_s",], aes(x = variable, y = value), width = 0.1, color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
  geom_violin(data = data_fit_melt[data_fit_melt$variable == "R_q",], aes(x = variable, y = value), alpha = 0.7) +
  geom_boxplot(data = data_fit_melt[data_fit_melt$variable == "R_q",], aes(x = variable, y = value), width = 0.1, alpha = 0.7) +
  geom_violin(data = data_known_melt[data_known_melt$variable == "R_q",], aes(x = variable, y = value), color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
  geom_boxplot(data = data_known_melt[data_known_melt$variable == "R_q",], aes(x = variable, y = value), width = 0.1, color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
    ggtitle("Intervention Impact using Known and SMC-Fit parameters") +
  xlab("Parameter") +
  scale_x_discrete(labels = c(expression(paste("Quarantine (", R[Q], ")", sep="")), expression(paste("Symptom Monitoring (", R[S], ")", sep=""))))

#### Plots to compare Abs_Benefit and Rel_Benefit ####
ggplot() +
  geom_violin(data = data_fit_melt[data_fit_melt$variable == "Abs_Benefit",], aes(x = variable, y = value), alpha = 0.7) +
  geom_boxplot(data = data_fit_melt[data_fit_melt$variable == "Abs_Benefit",], aes(x = variable, y = value), width = 0.1, alpha = 0.7) +
  geom_violin(data = data_known_melt[data_known_melt$variable == "Abs_Benefit",], aes(x = variable, y = value), color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
  geom_boxplot(data = data_known_melt[data_known_melt$variable == "Abs_Benefit",], aes(x = variable, y = value), width = 0.1, color = "forestgreen", fill = "lightgreen", alpha = 0.4) +  
  geom_violin(data = data_fit_melt[data_fit_melt$variable == "Rel_Benefit",], aes(x = variable, y = value), alpha = 0.7) +
  geom_boxplot(data = data_fit_melt[data_fit_melt$variable == "Rel_Benefit",], aes(x = variable, y = value), width = 0.1, alpha = 0.7) +
  geom_violin(data = data_known_melt[data_known_melt$variable == "Rel_Benefit",], aes(x = variable, y = value), color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
  geom_boxplot(data = data_known_melt[data_known_melt$variable == "Rel_Benefit",], aes(x = variable, y = value), width = 0.1, color = "forestgreen", fill = "lightgreen", alpha = 0.4) +
  ylim(-0.5, 1) +
  ggtitle("Comparative Efficacy using Known and SMC-Fit parameters") +
  xlab("Parameter") +
  scale_x_discrete(labels = c(expression(paste("Absolute Benefit (", R[S]-R[Q], ")", sep="")), expression(paste("Relative Benefit (", frac(R[S]-R[Q],R[S]), ")", sep=""))))

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_SMC_verification.RData", sep=""))

#### Test particle_filter_fcn ####
# Ranges for particle filter
T_lat_offset.min <- -5
T_lat_offset.max <- 5
d_inf.min <- 1
d_inf.max <- 40
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

dir = c("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/")

ks_stats <- particle_filter_fcn(T_lat_offset.max = T_lat_offset.max, T_lat_offset.min = T_lat_offset.min, d_inf.max = d_inf.max, d_inf.min = d_inf.min, pi_t_triangle_center.max = pi_t_triangle_center.max, pi_t_triangle_center.min = pi_t_triangle_center.min, parms_serial_interval = parms_serial_interval,  dir = dir, ks_conv_criteria = 0.001, disease_name = "SimvirusA", n_pop = 2000, times = 1, SMC_times = 10, perturb_initial = 1/2500, perturb_final = 1/7500)

sd(ks_stats$ks_conv_stat)/mean(ks_stats$ks_conv_stat)

#### What n_pop do we need to get a stable estimate for each set of parameters? ####
# Using the same input parameters each time (very little pertubation and range of inputs), see how the standard error of KS changes with each iteration.
# This will be helpful to determine how stable the KS measurement is and how large n_pop from each parameter set must be

n_pop_vector <- c(200, 400, 600, 800, 1000, 1200, 1400, 2000)
times <- 20
num_gen <- 2
df <- data.frame(cbind(n_pop_vector, rep(NA, length(n_pop_vector))))
names(df) <- c("n_pop", "stderr")
df$mean <- NA

T_lat_offset.min <- 2
T_lat_offset.max <- 2
d_inf.min <- 20
d_inf.max <- 20
pi_t_triangle_center.min <- .2
pi_t_triangle_center.max <- .2

for (i in 1:length(n_pop_vector)){
  outputs <- particle_filter_fcn(T_lat_offset.max = T_lat_offset.max, T_lat_offset.min = T_lat_offset.min, d_inf.max = d_inf.max, d_inf.min = d_inf.min, pi_t_triangle_center.max = pi_t_triangle_center.max, pi_t_triangle_center.min = pi_t_triangle_center.min, parms_serial_interval = parms_serial_interval,  dir = dir, ks_conv_criteria = 0.001, disease_name = "SimvirusA", n_pop = n_pop_vector[i], times = 1, num_generations = num_gen, SMC_times = times, perturb_initial = 1/2500, perturb_final = 1/7500)
  df[i, "stderr"] <- sd(outputs$ks_conv_stat)/mean(outputs$ks_conv_stat)
  df[i, "mean"] <- mean(outputs$ks_conv_stat)
}
layout(c(1))
plot(x = df$n_pop, y = df$stderr, xlab = "Number of individuals in initial chain", ylab = "Standard Error of KS")
plot(x = df$n_pop, y = df$mean, xlab = "Number of individuals in initial chain", ylab = "Mean of median KS", ylim = c(0, 0.5))

# Conclusion: n_pop needs to be 1000 for a stable estimate of KS
# This gets very slow though
# so consider increasing n_pop as you run more simulations and want to get a more precise KS value?

