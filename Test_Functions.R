#### Header ####
# Test Functions for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 19, 2015

#### Load Libraries ####
library(MASS)
library(parallel)
library(magrittr)
library(ggplot2)

#### Source Functions ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

#### Define a sample of initial parameters ####
n_pop = 500

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent", 1)
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target", "T_inc_stretch")

parms_T_lat = list("triangle", 999, 999, 999, -1, "T_inc")
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

gamma <- 0.9

#### Test Create_Pop ####
Pop <- Create_Pop(n_pop=10000, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma)
Pop[1,]
hist(Pop$T_inc)
hist(Pop$T_lat)
hist(Pop$R_0)
cat(length(which(Pop[,"R_0"]<0)), "individuals have an R_0 < 0")

parms_T_inc$T_inc_stretch <- 1

#### Test observe_and_isolate_fcn ####
Pop <- observe_and_isolate_fcn(Pop, intervention = "u")
sapply(Pop, summary)
Pop[1,]

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
dispersion = 20
children <- infection_times_fcn(Pop[1,"T_lat"], Pop[1,"d_inf"], Pop[1,"t_iso"], Pop[1,"t_obs"], Pop[1, "R_0"], Pop[1, "R_0_hsb_adjusted"], gamma=gamma, distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, intervention = "s", background_intervention = "u", dispersion = dispersion)
plot(children, main = "Children of person 1", xlab = "days since onset of infectiousness")

#### Test infection_times_fcn_2 ####
dispersion = 20
children <- infection_times_fcn_2(Pop[1,"T_lat"], Pop[1,"d_inf"], Pop[1,"t_iso"], Pop[1,"t_obs"], Pop[1, "R_0"], Pop[1, "R_0_hsb_adjusted"], gamma=gamma, distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, intervention = "s", background_intervention = "u", dispersion = dispersion)
plot(children, main = "Children of person 1", xlab = "days since onset of infectiousness")

#### Test children_list_fcn ####
dispersion = 1
start <- proc.time()
children_list <- children_list_fcn(Pop, pi_t_distribution=parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma=0.9, intervention = "u", background_intervention = "u", dispersion = dispersion)
proc.time() - start
cat('The hour(s) on which the first person transmitted is(are)',as.integer(which(children_list[[1]]==1)) )
# plot(children_list[[1]], xlab="hour")
Num_Infected <- unlist(lapply(children_list, sum))
hist(Num_Infected, breaks = seq(0, max(Num_Infected), max(1,Num_Infected/10)), xlab = "Number of Offspring Generated", ylab = "Number of Infectors", main = paste("Dispersion Parameter ( k =",1/(dispersion), ")"))
cat('The Effective Reproductive Number is', mean(Num_Infected))
sum(Num_Infected)

#### Test dispersion parameter within children_list_fcn ####
dispersion_options = c(1, 2, 10, 100)
df_disp <- data.frame(disp = rep(dispersion_options, each = 1000), offspring = rep(0:999, times = length(dispersion_options)), count = NA )
for (disp in dispersion_options){
  children_list <- children_list_fcn(Pop, pi_t_distribution=parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma=0.9, intervention = "u", background_intervention = "u", dispersion = disp)
  Num_Infected <- unlist(lapply(children_list, sum))
  df_disp[df_disp$disp == disp,"count"] <- hist(Num_Infected, breaks = 0:1000)$count
}

ggplot(df_disp, aes(x = offspring, y = count, color = factor(1/disp), group = disp)) + geom_line(size = 1, alpha = .8) + scale_x_continuous(limits = c(0,10), breaks = 0:10, name = "Offspring Count") + theme_bw() + scale_color_brewer(type = "div", palette = 4, name = "Dispersion Factor (k)") + theme(panel.grid.minor = element_blank()) + ylab("Number of Agents\n(N = 10,000)") + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/dispersion.tiff", width = 3, height = 3, units = "in")

# The following should all be equal
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1000, lambda = 2, d = 1))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(100, lambda = 2, d = 10))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(10, lambda = 2, d = 100))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 1000))))

# If we want each person to have an average R of 2, then having 1000 hours of infectiousness or just 1 doesn't matter for the dispersion constant
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1000, lambda = 2/1000, d = 100)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(1000, lambda = 2/1000, d = 100))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(100, lambda = 2/100, d = 100)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(100, lambda = 2/100, d = 100))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(10, lambda = 2/10, d = 100)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(10, lambda = 2/10, d = 100))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 100)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 100))))

# Should have higher dispersion when d is higher
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 100)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 100))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 10)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 10))))
sd(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 1)))); mean(sapply(rep(1, 100000), function(x) sum(rpois.od(1, lambda = 2, d = 1))))

# 6-panel figure "Fig S4.tiff"
parms_R_0 = list("uniform", 2, 2, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

par(mfrow = c(2,3))

Pop <- Create_Pop(n_pop=100000, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma) %>%
  observe_and_isolate_fcn(intervention = "u")

for (dispersion in c(1, 2, 5, 10, 20, 100)){
  children_list <- children_list_fcn(Pop, pi_t_distribution=parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma=0.9, intervention = "u", background_intervention = "u", dispersion = dispersion)
  
  Num_Infected <- unlist(lapply(children_list, sum))
  hist(Num_Infected, breaks = seq(0, max(Num_Infected)), xlim = c(0, 20), ylim = c(0, 100000), 
       main = paste("k =", 1/dispersion),
       xlab = "Number of Infections")
  text(x= 10, y=50000, labels = paste("Mean R =", round(mean(Num_Infected), 2), "\nVariance R =", round(sd(Num_Infected)^2, 2), "\n95th percentile =", sort(Num_Infected)[length(Num_Infected)*0.95]))
  
}

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
head(Pop_2)

#### Test to see if overdispersion leads to faster extinction ####

# Set parameters
parms_R_0$parm1 <- parms_R_0$parm2 <- 0.9
reps = 1000
dispersion_values <- c(1,2,5, 10, 100)
R0_values <- c(0.75)

# Create a data frame to store data
extinct_df <- data.frame(disp = rep(dispersion_values, each = reps*length(R0_values)), 
                         R0 = rep(rep(R0_values, each = reps), times = length(dispersion_values)),
                         rep = rep(1:reps, times = length(dispersion_values)*length(R0_values)),
                         generations = NA)

# Create population
Pop_alpha <- Create_Pop(n_pop=1, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma)


for (i in 1:nrow(extinct_df)){
    
  # Initialize
  Num_Infected = 999
  counter = 0
  Pop = Pop_alpha
  parms_R_0$parm1 <- parms_R_0$parm2 <- extinct_df[i, "RO"]
  
  while (Num_Infected > 0 & counter < 20){
    Pop <- observe_and_isolate_fcn(Pop, intervention = intervention)
    children_list <- children_list_fcn(Pop, parms_pi_t$distribution, parms_pi_t$triangle_center, gamma, intervention = "u", background_intervention = "u" , dispersion = extinct_df[i,"disp"])
    Num_Infected <- unlist(lapply(children_list, sum))
    
    if (Num_Infected > 0){
      counter = counter  + 1
    }
  }
  # cat(".")
  
  extinct_df[i, "generations"] = counter
}

# Helper functions
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="black", geom=geom, ...)
}
stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 1, shape = 15,...)
}


ggplot(extinct_df, aes(x = factor(disp), y = generations, color = generations)) +
  # facet_grid(.~R0) +
  geom_jitter(size = 0.5) + 
  scale_x_discrete(name = "Dispersion Factor (k)", breaks = dispersion_values, labels = 1/dispersion_values) + 
  stat_sum_df("mean_cl_normal", geom = "errorbar") +
  stat_sum_single(mean) +
  theme_bw() +
  scale_color_gradient(low = "grey", high = "blue") +
  guides(color = FALSE) + 
  scale_y_continuous(name = "Generations until Extinction") +
  theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "Extinction_R0.9.tiff", width = 4, height = 2, units = "in")

# Summarize 
extinct_df_summary <- data.frame(disp = rep(dispersion_values, each = length(R0_values)), 
                                 R0 = rep(R0_values, times = length(dispersion_values)),
                                 median = NA,
                                 lower_quartile = NA,
                                 upper_quartile = NA,
                                 mean = NA)
for (i in 1:nrow(extinct_df_summary)){
  extinct_df_summary[i, c("median", "lower_quartile", "upper_quartile")] <- quantile(extinct_df[extinct_df$disp == extinct_df_summary[i, "disp"] &
                                                                    extinct_df$R0 == extinct_df_summary[i, "R0"], 
                                                                  "generations"])[c("25%", "50%", "75%")]
  extinct_df_summary[i, c("mean")] <- mean(extinct_df[extinct_df$disp == extinct_df_summary[i, "disp"] &
                                                                                                  extinct_df$R0 == extinct_df_summary[i, "R0"], 
                                                                                                "generations"])
}

ggplot(extinct_df_summary, aes(x = 1/disp, y = mean, group = R0, color = R0)) + geom_line()


#### Test three generations of infection ####
dispersion = 1
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

intervention = "q"
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

start <- proc.time()
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
proc.time() - start

#### Plot Results from alpha, beta, gamma runs ####
plot((Pop_alpha$t_iso - Pop_alpha$T_inc)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), ylab = "Days from Symptoms to Isolation", xlab = "Person")
plot((Pop_beta$t_iso - Pop_beta$T_inc)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), ylab = "Days from Symptoms to Isolation", xlab = "Person")
plot((Pop_gamma$t_iso - Pop_gamma$T_inc)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), ylab = "Days from Symptoms to Isolation", xlab = "Person")

plot((Pop_alpha$t_iso - Pop_alpha$t_obs)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$t_obs)/24), ylab = "Days from Observation to Isolation", xlab = "Person")
plot((Pop_beta$t_iso - Pop_beta$t_obs)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$t_obs)/24), ylab = "Days from Observation to Isolation", xlab = "Person")
plot((Pop_gamma$t_iso - Pop_gamma$t_obs)/24, ylim = c(0, max(Pop_alpha$t_iso - Pop_alpha$t_obs)/24), ylab = "Days from Observation to Isolation", xlab = "Person")

hist((Pop_alpha$t_iso - Pop_alpha$T_inc)/24, xlim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), xlab ="Days from Symptoms to Isolation")
hist((Pop_beta$t_iso - Pop_beta$T_inc)/24, xlim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), xlab ="Days from Symptoms to Isolation")
hist((Pop_gamma$t_iso - Pop_gamma$T_inc)/24, xlim = c(0, max(Pop_alpha$t_iso - Pop_alpha$T_inc)/24), xlab ="Days from Symptoms to Isolation")

plot(apply(Pop_alpha, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))), Pop_alpha$t_iso, xlim=c(0, max(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))))), ylim=c(0, max(Pop_gamma$t_iso)))
plot(apply(Pop_beta, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))), Pop_beta$t_iso, xlim=c(0, max(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))))), ylim=c(0, max(Pop_gamma$t_iso)))
plot(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))), Pop_gamma$t_iso, xlim=c(0, max(apply(Pop_gamma, 1, function(x) min(as.numeric(x['T_inc']), as.numeric(x['T_lat']))))), ylim=c(0, max(Pop_gamma$t_iso)))

#### Test repeat_call_fcn ####
start <- proc.time()
dispersion = 2
In_Out <- repeat_call_fcn(n_pop = 200, 
                          parms_T_inc = parms_T_inc, 
                          parms_T_lat = parms_T_lat, 
                          parms_d_inf = parms_d_inf, 
                          parms_d_symp = parms_d_symp, 
                          parms_R_0 = parms_R_0, 
                          parms_epsilon = parms_epsilon, 
                          parms_pi_t = parms_pi_t,
                          num_generations = 5,
                          background_intervention="u",
                          subseq_interventions="q",
                          gamma=gamma,
                          prob_CT = prob_CT,
                          parms_CT_delay = parms_CT_delay,
                          parms_serial_interval = parms_serial_interval,
                          dispersion = dispersion, 
                          cap_pop = FALSE)
proc.time() - start
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

#### Test magrittr pipeline ####
dog <- c(5,5) %>% sum %>% class

Pop <- Create_Pop(n_pop=500, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma)

Pop2 <- Pop %>% observe_and_isolate_fcn(intervention = "u")
