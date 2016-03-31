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
n_pop = 2000
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
sd(ks_stats)/mean(ks_stats)

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

#### Use SMC algorithm to fit pi_t, d_inf, and T_lat_offset ####
dir = "~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/"

# Ranges for particle filter
T_lat_offset.min <- -5
T_lat_offset.max <- 5
d_inf.min <- 1
d_inf.max <- 40
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

# Restricted
T_lat_offset.min <- 0
T_lat_offset.max <- 0
d_inf.min <- 10
d_inf.max <- 10
pi_t_triangle_center.min <- 0.5
pi_t_triangle_center.max <- 0.5

outputs <- particle_filter_fcn(T_lat_offset.max = T_lat_offset.max, T_lat_offset.min = T_lat_offset.min, d_inf.max = d_inf.max, d_inf.min = d_inf.min, pi_t_triangle_center.max = pi_t_triangle_center.max, pi_t_triangle_center.min = pi_t_triangle_center.min, parms_serial_interval = parms_serial_interval,  dir = dir, ks_conv_criteria = 0.05, disease_name = "SimvirusA", n_pop = 1000, times = 1000, num_generations = 2, SMC_times = 15, perturb_initial = 1/25, perturb_final = 1/75)

head(outputs$data)
head(outputs$ks_conv_stat)

# Explore bias and precision
df_bias_precision <- data.frame(matrix(NA, nrow = 3, ncol = 5))
names(df_bias_precision) <- c("variable", "true", "mean", "lower_ci", "upper_ci")
df_bias_precision$variable <- c("T_lat_offset", "d_inf", "pi_t_triangle_center")
df_bias_precision$true <- c(parms_T_lat$anchor_value, parms_d_inf$parm2, parms_pi_t$triangle_center)
df_bias_precision$mean <- c(mean(outputs$data$T_lat_offset), mean(outputs$data$d_inf), mean(outputs$data$pi_t_triangle_center))
df_bias_precision$lower_ci <- c(mean(outputs$data$T_lat_offset) - 1.96*sd(outputs$data$T_lat_offset),
                                mean(outputs$data$d_inf) - 1.96*sd(outputs$data$d_inf),
                                mean(outputs$data$pi_t_triangle_center) - 1.96*sd(outputs$data$pi_t_triangle_center))
df_bias_precision$upper_ci <- c(mean(outputs$data$T_lat_offset) + 1.96*sd(outputs$data$T_lat_offset),
                                mean(outputs$data$d_inf) + 1.96*sd(outputs$data$d_inf),
                                mean(outputs$data$pi_t_triangle_center) + 1.96*sd(outputs$data$pi_t_triangle_center))

df_gg <- melt(outputs$data)
head(df_gg)

ggplot(df_gg, aes(x = factor(variable), y = value)) +
  geom_violin(data = df_gg[df_gg$variable == "T_lat_offset",], aes(y = value) ) +
  geom_violin(data = df_gg[df_gg$variable == "d_inf",], aes(y = value) ) +
  geom_violin(data = df_gg[df_gg$variable == "pi_t_triangle_center",], aes(y = value) ) +
  geom_point(data = df_bias_precision, aes(x = variable, y = true), color = "darkred", shape = 8, size = 3)

#### Compare intervention performance for input and fit parameters ####

