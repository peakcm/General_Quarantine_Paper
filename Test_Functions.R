#### Header ####
# Test Functions for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 19, 2015

#### Load Libraries ####
library(MASS)
library(parallel)
library(magrittr)

#### Define a sample of initial parameters ####
n_pop = 500

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_R_0 = list("uniform", 1.72, 1.94, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("confit", 7.5, 6.8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

prob_CT <- 0.3 #27.4 to 31.1% of case patients were in the contact registry before identification. around 30% of new cases reported having any contacts they could have infected 

parms_d_inf = list("triangle", 3, 8, 6.5, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 999, 999, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

background_intervention <- "u"

#### Test Create_Pop ####
Pop <- Create_Pop(n_pop=500, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma)
Pop[1,]
hist(Pop$T_inc)
hist(Pop$T_lat)
hist(Pop$R_0)
cat(length(which(Pop[,"R_0"]<0)), "individuals have an R_0 < 0")

#### Test observe_and_isolate_fcn ####
Pop <- observe_and_isolate_fcn(Pop, intervention = "u")

#### Test pi_t_fcn ####
pi_t <- pi_t_fcn(Pop[1,"T_lat"], Pop[1,"d_inf"], Pop[1,"t_iso"], Pop[1,"t_obs"], Pop[1, "R_0"], Pop[1, "R_0_hsb_adjusted"], 
                 gamma=gamma, 
                 distribution = parms_pi_t$distribution,
                 triangle_center = parms_pi_t$triangle_center,
                 intervention = "u",
                 background_intervention = "u")
plot(pi_t)
cat("The AUC for person 1 is", sum(pi_t), "\nThe R_0 for person 1 is", Pop[1, "R_0"])

#### Test infection_times_fcn ####
dispersion = 4
children <- infection_times_fcn(Pop[1,"T_lat"], Pop[1,"d_inf"], Pop[1,"t_iso"], Pop[1,"t_obs"], Pop[1, "R_0"], Pop[1, "R_0_hsb_adjusted"], gamma=gamma, distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, intervention = "s", background_intervention = "u", dispersion = dispersion)
plot(children, main = "Children of person 1", xlab = "days since onset of infectiousness")

#### Test children_list_fcn ####
dispersion = 5
children_list <- children_list_fcn(Pop, pi_t_distribution=parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma=0.9, intervention = "u", background_intervention = "u", dispersion = dispersion)
cat('The hour(s) on which the first person transmitted is(are)',as.integer(which(children_list[[1]]==1)) )
plot(children_list[[1]], xlab="hour")
Num_Infected <- unlist(lapply(children_list, sum))
hist(Num_Infected, breaks = seq(0, max(Num_Infected)))
cat('The Effective Reproductive Number is', mean(Num_Infected))
sum(Num_Infected)

#### Test next_generation_fcn ####
Pop_2 <- next_generation_fcn(Pop = Pop,
                             children_list = children_list,
                             parms_T_inc = parms_T_inc,
                             parms_T_lat = parms_T_lat,
                             parms_d_inf = parms_d_inf,
                             parms_d_symp = parms_d_symp,
                             parms_R_0 =  parms_R_0,
                             parms_epsilon = parms_epsilon,
                             generation = 2,
                             prob_CT = prob_CT,
                             parms_CT_delay = parms_CT_delay,
                             gamma = gamma,
                             n_pop = n_pop,
                             cap_pop = FALSE)
dim(Pop_2)

#### Test three generations of infection ####
dispersion = 2
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

intervention = "hsb"
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

#### Plot Results from alpha, beta, gamma runs ####
plot((Pop_alpha$t_iso - Pop_alpha$T_inc)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), ylab = "Days from Symptoms to Isolation", xlab = "Person")
plot((Pop_beta$t_iso - Pop_beta$T_inc)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), ylab = "Days from Symptoms to Isolation", xlab = "Person")
plot((Pop_gamma$t_iso - Pop_gamma$T_inc)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), ylab = "Days from Symptoms to Isolation", xlab = "Person")

hist((Pop_alpha$t_iso - Pop_alpha$T_inc)/24, xlim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), xlab ="Days from Symptoms to Isolation")
hist((Pop_beta$t_iso - Pop_beta$T_inc)/24, xlim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), xlab ="Days from Symptoms to Isolation")
hist((Pop_gamma$t_iso - Pop_gamma$T_inc)/24, xlim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), xlab ="Days from Symptoms to Isolation")

plot(apply(Pop_alpha, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))), Pop_alpha$t_iso, xlim=c(0, max(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))))), ylim=c(0, max(Pop_gamma$t_iso)))
plot(apply(Pop_beta, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))), Pop_beta$t_iso, xlim=c(0, max(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))))), ylim=c(0, max(Pop_gamma$t_iso)))
plot(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))), Pop_gamma$t_iso, xlim=c(0, max(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))))), ylim=c(0, max(Pop_gamma$t_iso)))

#### Test repeat_call_fcn ####
dispersion = 2
In_Out <- repeat_call_fcn(n_pop = n_pop, 
                          parms_T_inc = parms_T_inc, 
                          parms_T_lat = parms_T_lat, 
                          parms_d_inf = parms_d_inf, 
                          parms_d_symp = parms_d_symp, 
                          parms_R_0 = parms_R_0, 
                          parms_epsilon = parms_epsilon, 
                          parms_pi_t = parms_pi_t,
                          num_generations = 5,
                          background_intervention="u",
                          subseq_interventions="hsb",
                          gamma=gamma,
                          prob_CT = prob_CT,
                          parms_CT_delay = parms_CT_delay,
                          parms_serial_interval = parms_serial_interval,
                          dispersion = dispersion, 
                          cap_pop = FALSE)
In_Out$output
# plot(In_Out$output$R)

# input value for R_0
# In_Out$input$parms_R_0

#### Test serial_interval_fcn ####
layout(c(1))
# First, run the code block "Test three generations of infection" to create Pop_alpha and Pop_beta
serial_interval_fcn(Pop_alpha, Pop_beta, parms_serial_interval, plot="True")

#### Summarize distribution ####
layout(c(1))
summarize_dist_fcn(parms_serial_interval, lower = 0.025, upper = 0.975)
summarize_dist_fcn(parms_T_inc, lower = 0.025, upper = 0.975)

summarize_dist_fcn(parms_serial_interval, lower = 0.05, upper = 0.95)
summarize_dist_fcn(parms_T_inc, lower = 0.05, upper = 0.95)

#### Test overdispersion super spreading feature ####

parms_R_0 = list("uniform", 2, 2, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

par(mfrow = c(3,3))

Pop <- Create_Pop(n_pop=5000, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma) %>%
  observe_and_isolate_fcn(intervention = "u")

for (dispersion in 1:9){
  children_list <- children_list_fcn(Pop, pi_t_distribution=parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma=0.9, intervention = "u", background_intervention = "u", dispersion = dispersion)
  
  Num_Infected <- unlist(lapply(children_list, sum))
  hist(Num_Infected, breaks = seq(0, max(Num_Infected)), xlim = c(0, 20), ylim = c(0, 4000), 
       main = paste("Dispersion =", dispersion),
       xlab = "Number of Infections")
  text(x= 10, y=3000, labels = paste("Mean R =", mean(Num_Infected), "\n 95th percentile =", sort(Num_Infected)[length(Num_Infected)*0.95]))

}

#### Test magrittr pipeline ####
dog <- c(5,5) %>% sum %>% class

Pop <- Create_Pop(n_pop=500, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma)

Pop2 <- Pop %>% observe_and_isolate_fcn(intervention = "u")
