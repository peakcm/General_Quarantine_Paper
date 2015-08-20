#### Header ####
# Test Functions for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 19, 2015

#### Load Libraries ####
library(MASS)

#### Define a sample of initial parameters ####
# parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
# names(parms_T_inc) <- c("dist","parm1","parm2",  "parm3","parm3", "anchor_value", "anchor_target")

parms_T_inc = list("triangle", 1, 21, 10, "independent", "independent")
names(parms_T_inc) <- c("dist","parm1","parm2", "parm3", "anchor_value", "anchor_target")

parms_T_lat = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 5, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 5, 10, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("confit", 1.7, 1.2, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

gamma = 0.9

prob_CT = 1

background_intervention = "u"
intervention = "u"

parms_pi_t <- list("triangle", 0.20)
names(parms_pi_t) <- c("distribution","triangle_center")

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
children <- infection_times_fcn(Pop[1,"T_lat"], Pop[1,"d_inf"], Pop[1,"t_iso"], Pop[1,"t_obs"], Pop[1, "R_0"], Pop[1, "R_0_hsb_adjusted"], gamma=gamma, distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, intervention = "s", background_intervention = "u")
plot(children, main = "Children of person 1", xlab = "days since onset of infectiousness")

#### Test children_list_fcn ####
children_list <- children_list_fcn(Pop, pi_t_distribution=parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma=0.9, intervention = "u", background_intervention = "u")
cat('The hour(s) on which the first person transmitted is(are)',as.integer(which(children_list[[1]]==1)) )
plot(children_list[[1]], xlab="hour")
Num_Infected <- unlist(lapply(children_list, sum))
hist(Num_Infected, breaks = seq(0, max(Num_Infected)))
cat('The Effective Reproductive Number is', mean(Num_Infected))

#### Test next_generation_fcn ####
Pop_2 <- next_generation_fcn(Pop,
                             children_list,
                             parms_T_inc,
                             parms_T_lat,
                             parms_d_inf,
                             parms_d_symp,
                             parms_R_0,
                             parms_epsilon,
                             generation = 2,
                             prob_CT,
                             parms_CT_delay,
                             gamma)

#### Test three generations of infection ####
n_pop = 1000
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
if (nrow(Pop_alpha) > n_pop){Pop_alpha <- Pop_alpha[1:n_pop,]}
Pop_alpha <- observe_and_isolate_fcn(Pop_alpha, intervention = background_intervention)
children_list <- children_list_fcn(Pop_alpha, parms_pi_t$distribution, parms_pi_t$triangle_center, gamma, intervention = background_intervention, background_intervention)
Num_Infected <- unlist(lapply(children_list, sum))
cat('Generation 1 : n=', nrow(Pop_alpha), '. Effective Reproductive Number:', mean(Num_Infected), '. Number of infections:', sum(Num_Infected))

Pop_beta <- next_generation_fcn(Pop_alpha,
                                children_list,
                                parms_T_inc,
                                parms_T_lat,
                                parms_d_inf,
                                parms_d_symp,
                                parms_R_0,
                                parms_epsilon,
                                generation = 2,
                                prob_CT = prob_CT,
                                parms_CT_delay,
                                gamma=gamma,
                                n_pop)
if (nrow(Pop_beta) > n_pop){Pop_beta <- Pop_beta[1:n_pop,]}
Pop_beta <- observe_and_isolate_fcn(Pop_beta, intervention = intervention)
children_list <- children_list_fcn(Pop_beta, parms_pi_t$distribution, parms_pi_t$triangle_center, gamma, intervention = intervention, background_intervention)
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
if (nrow(Pop_gamma) > n_pop){Pop_gamma <- Pop_gamma[1:n_pop,]}
Pop_gamma <- observe_and_isolate_fcn(Pop_gamma, intervention = intervention)
children_list <- children_list_fcn(Pop_gamma, parms_pi_t$distribution, parms_pi_t$triangle_center, gamma, intervention = intervention, background_intervention)
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
In_Out <- repeat_call_fcn(n_pop=1000, 
                          parms_T_inc = parms_T_inc, 
                          parms_T_lat = parms_T_lat, 
                          parms_d_inf = parms_d_inf, 
                          parms_d_symp = parms_d_symp, 
                          parms_R_0 = parms_R_0, 
                          parms_epsilon = parms_epsilon, 
                          parms_pi_t = parms_pi_t,
                          num_generations = 5,
                          background_intervention="u",
                          subseq_interventions="u",
                          gamma=gamma,
                          prob_CT = prob_CT,
                          parms_CT_delay,
                          parms_serial_interval)
In_Out$output
# plot(In_Out$output$R)

# input value for R_0
# In_Out$input$parms_R_0

#### Test serial_interval_fcn ####
layout(c(1))
# First, run the code block "Test three generations of infection" to create Pop_alpha and Pop_beta
serial_interval_fcn(Pop_alpha, Pop_beta, parms_serial_interval, plot="True")

# SARS serial interval approximation
parms_serial_interval <- list("weibull", 2, 10)
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(dweibull(x, shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2),
      from=0, to=30, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")

# Ebola serial interval approximation from WHO
parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(dgamma(x, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2),
      from=0, to=40, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")

# Ebola serial interval approximation from Eichner
parms_serial_interval <- list("lognormal", 12.7, 4.31)
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(dlnorm(x, meanlog=log(parms_serial_interval$parm1), sdlog=log(sqrt(parms_serial_interval$parm2))),
      from=0, to=40, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")

# Pertussis serial interval approximation
parms_serial_interval <- list("weibull", 2, 22)
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(dweibull(x, shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2),
      from=0, to=84, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")
hist(rweibull(10000, parms_serial_interval$parm1, parms_serial_interval$parm2),breaks = seq(from=0, to=84, by=7))
summary(rweibull(10000, parms_serial_interval$parm1, parms_serial_interval$parm2))
