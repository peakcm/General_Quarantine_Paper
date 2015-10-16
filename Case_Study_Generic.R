#### Header #### 
# Case Study code for a generic SEIR acute directly transmitted disease
# Corey Peak
# Version 1.0
# October 16, 2015

#### Load Libraries ####
library(ppcor)
library(MASS)
library(lhs)
library(Hmisc)
library(sensitivity)
library(ggplot2)
library(reshape)
library(psych)

#### Disease: SARS ####

# Fixed Disease Parameters
parms_serial_interval <- list("weibull", 2, 10)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("lognormal", 4, 1.81, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 2, 15, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Ranges for particle filter
T_lat_offset.min <- -7
T_lat_offset.max <- 2
d_symp_lower.min <- 1
d_symp_lower.max <- 15
d_symp_width.min <- 0
d_symp_width.max <- 40
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

# Save and load workspaces

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/20150829_SARS_ParticleFilter.RData")
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150901_SARS_ParticleFilter_2000.RData')
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/20150829_SARS_ParticleFilter.RData')

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150904_SARS_HRSetting.RData")
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150904_SARS_HRSetting.RData')

# save.image('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150904_SARS_LRSetting.RData')
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150904_SARS_LRSetting.RData')

# save.image('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150904_plot.RData')
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/SARS/20150904_plot.RData')

#### Disease: Ebola ####

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("triangle", 1, 8, 2, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Ranges for particle filter
T_lat_offset.min <- -0.5
T_lat_offset.max <- 2
d_symp_peak.min <- 0
d_symp_peak.max <- 1
d_symp_max.min <- 15
d_symp_max.max <- 30
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

# Save and load workspaces

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/20150829_Ebola_ParticleFilter.RData")
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/20150829_Ebola_ParticleFilter.RData')

# save.image('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Ebola/20151001_plot.RData')
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Ebola/20151001_plot.RData')

#### Particle Filter with weights and threshold ####

# Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Initialize
n_pop = 1000
num_generations <- 4
times <- 1000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_lower","d_symp_upper","pi_t_triangle_center")
init_threshold <- 0.5
reduce <- 0.80
SMC_times <- 3
ks_conv_criteria <- 0.10
ks_conv_stat <- rep(NA, SMC_times+1)
subseq_interventions <- "u"
printing = FALSE

lhs <- maximinLHS(times, length(dimensions))
params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_symp_lower = lhs[,2]*(d_symp_lower.max - d_symp_lower.min) + d_symp_lower.min,
  d_symp_width = lhs[,3]*(d_symp_width.max - d_symp_width.min) + d_symp_width.min,
  pi_t_triangle_center = lhs[,4]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

start.time <- proc.time()

for (i in 1:times){
  if (printing == TRUE){cat('\nIteration',i, '\n')}
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_symp$parm1 <- params.set[i,"d_symp_lower"]
  parms_d_symp$parm2 <- parms_d_symp$parm1 + max(0, params.set[i,"d_symp_width"])
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  
  In_Out <- repeat_call_fcn(n_pop=n_pop, 
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
                            printing = printing)
  data[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
  data[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
}
print(proc.time() - start.time)

data$T_lat_offset <- params.set[,"T_lat_offset"]
data$d_symp_lower <- params.set[,"d_symp_lower"]
data$d_symp_width <- params.set[,"d_symp_width"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]

ks_conv_stat[1] <- sort(data$ks)[floor(times*0.975)] 

# Save scatterplots
pdf(file = paste("Iteration_1.pdf"))
layout(rbind(c(1,2,3),c(4,5,6),c(7,8,9)))
plot(x=data$pi_t_triangle_center, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
     ylim = c(T_lat_offset.min, T_lat_offset.max),
     xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max),
     main = paste("Iteration 1"))
plot(x=data$d_symp_lower, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
     ylim = c(T_lat_offset.min, T_lat_offset.max),
     xlim = c(d_symp_lower.min, d_symp_lower.max),
     main = paste("KS = ", round(ks_conv_stat[1],2)))
plot(x=data$d_symp_width, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
     ylim = c(T_lat_offset.min, T_lat_offset.max),
     xlim = c(d_symp_width.min, d_symp_width.max))  
plot(x=data$pi_t_triangle_center, y=data$d_symp_width, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
     ylim = c(d_symp_width.min, d_symp_width.max),
     xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max)) 
plot(x=data$d_symp_lower, y=data$d_symp_width, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
     ylim = c(d_symp_width.min, d_symp_width.max),
     xlim = c(d_symp_lower.min, d_symp_lower.max))
frame()
plot(x=data$pi_t_triangle_center, y=data$d_symp_lower, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
     ylim = c(d_symp_lower.min, d_symp_lower.max),
     xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max))
frame()
frame()
dev.off()

# Plot scatterplot with cumulative distributions
pairs.panels(data[3:6],pch=21,
             bg = rainbow(1000)[floor(data$ks*1000)+1], ellipses = FALSE, smooth=FALSE, rug=FALSE) 

T_lat_offset.perturb <- (max(data[,"T_lat_offset"]) - min(data[,"T_lat_offset"])) / 50
d_symp_lower.perturb <- (max(data[,"d_symp_lower"]) - min(data[,"d_symp_lower"])) / 50
d_symp_width.perturb <- (max(data[,"d_symp_width"]) - min(data[,"d_symp_width"])) / 50
pi_t_triangle_center.perturb <- (max(data[,"pi_t_triangle_center"]) - min(data[,"pi_t_triangle_center"])) / 50

threshold <- init_threshold

theta_pre_can <- sample(row.names(data[data$ks <= threshold,]), times, prob= (1/data[data$ks <= threshold,]$ks), replace=TRUE) #sample pre-candidate theta parameter sets from previous generation
T_lat_offset.theta <- data[theta_pre_can,"T_lat_offset"] + runif(times, min=-1*T_lat_offset.perturb, max=T_lat_offset.perturb) #perturb and propose
d_symp_lower.theta <- data[theta_pre_can,"d_symp_lower"] + runif(times, min=-1*d_symp_lower.perturb, max=d_symp_lower.perturb) #perturb and propose
d_symp_width.theta <- data[theta_pre_can,"d_symp_width"] + runif(times, min=-1*d_symp_width.perturb, max=d_symp_width.perturb) #perturb and propose
pi_t_triangle_center.theta <- data[theta_pre_can,"pi_t_triangle_center"] + runif(times, min=-1*pi_t_triangle_center.perturb, max=pi_t_triangle_center.perturb) #perturb and propose

d_symp_lower.theta[d_symp_lower.theta < 1] <- 1
d_symp_width.theta[d_symp_width.theta < 0] <- 0
pi_t_triangle_center.theta[pi_t_triangle_center.theta > 1] <- 1
pi_t_triangle_center.theta[pi_t_triangle_center.theta < 0] <- 0

SMC_break <- FALSE
SMC_counter <- 2
while (SMC_break == FALSE){
  for (i in 1:times){
    if (printing == TRUE){cat('\nIteration',i, '\n')}
    
    parms_T_lat$anchor_value <- T_lat_offset.theta[i]
    parms_d_symp$parm1 <- d_symp_lower.theta[i]
    parms_d_symp$parm2 <- parms_d_symp$parm1 + max(0, d_symp_width.theta[i])
    parms_pi_t$triangle_center <- pi_t_triangle_center.theta[i]
    
    In_Out <- repeat_call_fcn(n_pop=n_pop, 
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
                              printing = printing)
    if (subseq_interventions == background_intervention){
      data[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
  
  data$T_lat_offset <- T_lat_offset.theta
  data$d_symp_lower <- d_symp_lower.theta
  data$d_symp_width <- d_symp_width.theta
  data$pi_t_triangle_center <- pi_t_triangle_center.theta
  
  ks_conv_stat[SMC_counter] <- sort(data$ks)[floor(times*0.975)]  # Calculate the upper end of the inner 95% credible interval
  
  # Save scatterplots
  pdf(paste("Iteration_",SMC_counter,".pdf", sep = ""))
  layout(rbind(c(1,2,3),c(4,5,6),c(7,8,9)))
  plot(x=data$pi_t_triangle_center, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(T_lat_offset.min, T_lat_offset.max),
       xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max),
       main = paste("Iteration ", SMC_counter))
  plot(x=data$d_symp_lower, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(T_lat_offset.min, T_lat_offset.max),
       xlim = c(d_symp_lower.min, d_symp_lower.max),
       main = paste("KS = ", round(ks_conv_stat[SMC_counter],2)))
  plot(x=data$d_symp_width, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(T_lat_offset.min, T_lat_offset.max),
       xlim = c(d_symp_width.min, d_symp_width.max))  
  plot(x=data$pi_t_triangle_center, y=data$d_symp_width, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(d_symp_width.min, d_symp_width.max),
       xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max)) 
  plot(x=data$d_symp_lower, y=data$d_symp_width, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(d_symp_width.min, d_symp_width.max),
       xlim = c(d_symp_lower.min, d_symp_lower.max))
  frame()
  plot(x=data$pi_t_triangle_center, y=data$d_symp_lower, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(d_symp_lower.min, d_symp_lower.max),
       xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max))
  frame()
  frame()
  dev.off()
  
  # Plot scatterplot with cumulative distributions
  pairs.panels(data[3:6],pch=21,
               bg = rainbow(1000)[floor(data$ks*1000)+1], ellipses = FALSE, smooth=FALSE, rug=FALSE) 
  
  T_lat_offset.perturb <- (max(data[,"T_lat_offset"]) - min(data[,"T_lat_offset"])) / 25
  d_symp_lower.perturb <- (max(data[,"d_symp_lower"]) - min(data[,"d_symp_lower"])) / 25
  d_symp_width.perturb <- (max(data[,"d_symp_width"]) - min(data[,"d_symp_width"])) / 25
  pi_t_triangle_center.perturb <- (max(data[,"pi_t_triangle_center"]) - min(data[,"pi_t_triangle_center"])) / 25
  
  threshold <- max( ks_conv_criteria, threshold * reduce )
  
  theta_pre_can <- sample(row.names(data[data$ks <= threshold,]), times, prob= (1/data[data$ks <= threshold,]$ks), replace=TRUE) #sample pre-candidate theta parameter sets from previous generation
  T_lat_offset.theta <- data[theta_pre_can,"T_lat_offset"] + runif(times, min=-1*T_lat_offset.perturb, max=T_lat_offset.perturb) #perturb and propose
  d_symp_lower.theta <- data[theta_pre_can,"d_symp_lower"] + runif(times, min=-1*d_symp_lower.perturb, max=d_symp_lower.perturb) #perturb and propose
  d_symp_width.theta <- data[theta_pre_can,"d_symp_width"] + runif(times, min=-1*d_symp_width.perturb, max=d_symp_width.perturb) #perturb and propose
  pi_t_triangle_center.theta <- data[theta_pre_can,"pi_t_triangle_center"] + runif(times, min=-1*pi_t_triangle_center.perturb, max=pi_t_triangle_center.perturb) #perturb and propose
  
  d_symp_lower.theta[d_symp_lower.theta < 1] <- 1
  d_symp_width.theta[d_symp_width.theta < 0] <- 0
  pi_t_triangle_center.theta[pi_t_triangle_center.theta > 1] <- 1
  pi_t_triangle_center.theta[pi_t_triangle_center.theta < 0] <- 0
  
  if (ks_conv_stat[SMC_counter] <= ks_conv_criteria){
    SMC_break <- TRUE
    cat("Convergence acheived in", i, "iterations")
  }
  
  SMC_counter <- SMC_counter + 1
  if (SMC_counter >= SMC_times){
    SMC_break <- TRUE
    cat("Unable to converge by", SMC_counter, "SMC iterations")
  }
}

# write.table(ks_conv_stat, "20150829_SARS_ParticleFilter.csv")

#### Ranking of Intervention Sensitivities ####

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
times <- 2000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names

# Sample from posterior distributions for each parameter independently
params.set <- cbind(
  T_lat_offset = sample(data$T_lat_offset, size = times, replace = TRUE),
  d_symp_lower = sample(data$d_symp_lower, size = times, replace = TRUE),
  d_symp_width = sample(data$d_symp_width, size = times, replace = TRUE),
  pi_t_triangle_center = sample(data$pi_t_triangle_center, size = times, replace = TRUE) )

# Set range for other parameters to vary
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0", "dispersion")
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0.2
gamma.max <- 0.9
prob_CT.min <- 0.25
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 5
epsilon.min <- 0
epsilon.max <- 7
R_0.min <- 1
R_0.max <- 3
dispersion.min <- 1
dispersion.max <- 5

params.set <- cbind(params.set,
                    gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
                    prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
                    CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
                    epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
                    R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min,
                    dispersion = lhs[,6]*(dispersion.max - dispersion.min) + dispersion.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]))
  dispersion <- as.numeric(params.set[i, "dispersion"])
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 100
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 100
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 100
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
                              dispersion = dispersion)
    if (subseq_interventions == background_intervention){
      data[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data[,"Abs_Benefit"] <- data[,"R_s"] - data[,"R_q"]
data[,"Rel_Benefit"] <- data[,"Abs_Benefit"] / data[,"R_s"]
data[,"NNQ"] <- 1 / data[,"Abs_Benefit"]
data[data$NNQ < 1,"NNQ"] <- 1
data[data$NNQ > 9999,"NNQ"] <- 9999
data[data$NNQ == Inf,"NNQ"] <- 9999
data[,"Abs_Benefit_per_Qday"] <- data[,"Abs_Benefit"] / data[,"obs_to_iso_q"]
data$gamma_vector <- params.set[,"gamma"]
data$prob_CT_vector <- params.set[,"prob_CT"]
data$CT_delay_vector <- params.set[,"CT_delay"]
data$epsilon_vector <- params.set[,"epsilon"]
data$R_0_vector <- params.set[,"R_0"]
data$dispersion <- params.set[,"dispersion"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data$T_lat_offset_vector <- params.set[,"T_lat_offset"]

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$gamma_vector, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$gamma_vector, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$gamma_vector, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$gamma_vector, log10(data$NNQ))
plot(data$gamma_vector, data$Abs_Benefit)
plot(data$gamma_vector, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$prob_CT_vector, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$prob_CT_vector, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$prob_CT_vector, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$prob_CT_vector, log10(data$NNQ))
plot(data$prob_CT_vector, data$Abs_Benefit)
plot(data$prob_CT_vector, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$CT_delay_vector, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$CT_delay_vector, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$CT_delay_vector, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$CT_delay_vector, log10(data$NNQ))
plot(data$CT_delay_vector, data$Abs_Benefit)
plot(data$CT_delay_vector, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$epsilon_vector, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$epsilon_vector, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$epsilon_vector, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$epsilon_vector, log10(data$NNQ))
plot(data$epsilon_vector, data$Abs_Benefit)
plot(data$epsilon_vector, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$R_0_vector, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$R_0_vector, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$R_0_vector, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$R_0_vector, log10(data$NNQ))
plot(data$R_0_vector, data$Abs_Benefit)
plot(data$R_0_vector, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$pi_t_triangle_center, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$pi_t_triangle_center, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$pi_t_triangle_center, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$pi_t_triangle_center, log10(data$NNQ))
plot(data$pi_t_triangle_center, data$Abs_Benefit)
plot(data$pi_t_triangle_center, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$T_lat_offset_vector, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$T_lat_offset_vector, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$T_lat_offset_vector, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$T_lat_offset_vector, log10(data$NNQ))
plot(data$T_lat_offset_vector, data$Abs_Benefit)
plot(data$T_lat_offset_vector, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$dispersion, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$dispersion, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$dispersion, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$dispersion, log10(data$NNQ))
plot(data$dispersion, data$Abs_Benefit)
plot(data$dispersion, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(t(c(1,2)))
plot(data$pi_t_triangle_center, data$ks)
plot(data$T_lat_offset_vector, data$ks)

decile_plot_fcn(data, params.set)

data_store <- data

# Partial Rank Correlation using "ppcor" package
require(ppcor)
pcor(data[,c("R_s",names(data)[11:length(names(data))])], method=c("spearman"))$estimate[1,]
pcor(data[,c("R_s",names(data)[11:length(names(data))])], method=c("spearman"))$p.value[1,]
pcor(data[,c("Abs_Benefit",names(data)[11:length(names(data))])], method=c("spearman"))$estimate[1,]
pcor(data[,c("Abs_Benefit",names(data)[11:length(names(data))])], method=c("spearman"))$p.value[1,]
pcor(data[,c("Abs_Benefit_per_Qday", names(data)[11:length(names(data))])], method=c("spearman"))$estimate[1,]
pcor(data[,c("Abs_Benefit_per_Qday", names(data)[11:length(names(data))])], method=c("spearman"))$p.value[1,]

# Partial Rank Correlation using "sensitivity" package
require(sensitivity)
bonferroni.alpha <- 0.05/length(dimensions)
prcc_data <- data.frame(matrix(rep(NA, 7*length(names)*length(dimensions)), ncol=7)) 
names(prcc_data) <- c("output","parameter","coef","bias","stderr","CImin","CImax")
prcc_data$output <- rep(names, each = length(dimensions))
prcc_data$parameter <- rep(dimensions, times = length(names))
for (output in names){
  prcc <- pcc(data[,(length(names)+1):length(names(data))], data[,output], nboot = 100, rank=TRUE, conf=1-bonferroni.alpha)
  summary <- print(prcc)
  prcc_data[prcc_data$output == output,3:7] <- summary
}

require(ggplot2)
ggplot(prcc_data, aes(x = parameter, y= coef)) +
  facet_grid(output ~ .) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  geom_errorbar(data = prcc_data, aes(ymin = CImin, ymax = CImax), width = 0.1)

#### Case Study in High Resource Setting ####

# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0.9, 0.9, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 2000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.hr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_symp_lower = data[sample, "d_symp_lower"],
  d_symp_width = data[sample, "d_symp_width"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1, max = 10))

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_symp$parm1 <- params.set[i,"d_symp_lower"]
  parms_d_symp$parm2 <- parms_d_symp$parm1 + max(0, params.set[i,"d_symp_width"])
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
                              dispersion = dispersion)
    if (subseq_interventions == background_intervention){
      data.hr[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.hr[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.hr[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
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
data.hr$d_symp_lower_vector <- params.set[,"d_symp_lower"]
data.hr$d_symp_width_vector <- params.set[,"d_symp_width"]
data.hr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.hr$R_0_vector <- params.set[,"R_0"]
data.hr$T_lat_offset_vector <- params.set[,"T_lat_offset"]
data.hr$dispersion_vector <- params.set[,"dispersion"]

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$R_0_vector, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$R_0_vector, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$R_0_vector, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$R_0_vector, log10(data.hr$NNQ))
plot(data.hr$R_0_vector, data.hr$Abs_Benefit)
plot(data.hr$R_0_vector, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$pi_t_triangle_center, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$pi_t_triangle_center, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$pi_t_triangle_center, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$pi_t_triangle_center, log10(data.hr$NNQ))
plot(data.hr$pi_t_triangle_center, data.hr$Abs_Benefit)
plot(data.hr$pi_t_triangle_center, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$T_lat_offset_vector, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$T_lat_offset_vector, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$T_lat_offset_vector, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$T_lat_offset_vector, log10(data.hr$NNQ))
plot(data.hr$T_lat_offset_vector, data.hr$Abs_Benefit)
plot(data.hr$T_lat_offset_vector, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$d_symp_lower_vector, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_symp_lower_vector, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_symp_lower_vector, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_symp_lower_vector, log10(data.hr$NNQ))
plot(data.hr$d_symp_lower_vector, data.hr$Abs_Benefit)
plot(data.hr$d_symp_lower_vector, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$d_symp_width_vector, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_symp_width_vector, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_symp_width_vector, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_symp_width_vector, log10(data.hr$NNQ))
plot(data.hr$d_symp_width_vector, data.hr$Abs_Benefit)
plot(data.hr$d_symp_width_vector, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$dispersion_vector, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$dispersion_vector, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$dispersion_vector, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$dispersion_vector, log10(data.hr$NNQ))
plot(data.hr$dispersion_vector, data.hr$Abs_Benefit)
plot(data.hr$dispersion_vector, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

summary(data.hr$R_0)
summary(data.hr$R_hsb)
summary(data.hr$R_s)
summary(data.hr$R_q)
summary(data.hr$Abs_Benefit)
summary(data.hr$NNQ)

summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 , "R_0"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "R_hsb"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "R_s"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "R_q"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "Abs_Benefit"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "NNQ"])

quantile(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "R_s"], c(0.025, 0.50, 0.975))
quantile(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6, "R_q"], c(0.025, 0.50, 0.975))

summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset_vector < 0, "R_0"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset_vector < 0, "R_hsb"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset_vector < 0, "R_s"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset_vector < 0, "R_q"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset_vector < 0, "Abs_Benefit"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset_vector < 0, "NNQ"])

#### Case Study in Low Resource Setting ####

# Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

parms_epsilon = list("uniform", 2, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 2000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.lr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.lr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_symp_lower = data[sample, "d_symp_lower"],
  d_symp_width = data[sample, "d_symp_width"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1, max = 10))

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- params.set[i,"T_lat_offset"]
  parms_d_symp$parm1 <- params.set[i,"d_symp_lower"]
  parms_d_symp$parm2 <- parms_d_symp$parm1 + max(0, params.set[i,"d_symp_width"])
  parms_pi_t$triangle_center <- params.set[i,"pi_t_triangle_center"]
  parms_R_0$parm1 <- params.set[i,"R_0"]
  parms_R_0$parm2 <- params.set[i,"R_0"]
  dispersion <- params.set[i, "dispersion"]
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
                              dispersion = dispersion)
    if (subseq_interventions == background_intervention){
      data.lr[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.lr[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.lr[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
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
data.lr$d_symp_lower_vector <- params.set[,"d_symp_lower"]
data.lr$d_symp_width_vector <- params.set[,"d_symp_width"]
data.lr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.lr$R_0_vector <- params.set[,"R_0"]
data.lr$T_lat_offset_vector <- params.set[,"T_lat_offset"]
data.lr$dispersion_vector <- params.set[,"dispersion"]

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$R_0_vector, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$R_0_vector, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$R_0_vector, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$R_0_vector, log10(data.lr$NNQ))
plot(data.lr$R_0_vector, data.lr$Abs_Benefit)
plot(data.lr$R_0_vector, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$pi_t_triangle_center, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$pi_t_triangle_center, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$pi_t_triangle_center, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$pi_t_triangle_center, log10(data.lr$NNQ))
plot(data.lr$pi_t_triangle_center, data.lr$Abs_Benefit)
plot(data.lr$pi_t_triangle_center, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$T_lat_offset_vector, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$T_lat_offset_vector, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$T_lat_offset_vector, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$T_lat_offset_vector, log10(data.lr$NNQ))
plot(data.lr$T_lat_offset_vector, data.lr$Abs_Benefit)
plot(data.lr$T_lat_offset_vector, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$d_symp_lower_vector, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_symp_lower_vector, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_symp_lower_vector, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_symp_lower_vector, log10(data.lr$NNQ))
plot(data.lr$d_symp_lower_vector, data.lr$Abs_Benefit)
plot(data.lr$d_symp_lower_vector, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$d_symp_width_vector, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_symp_width_vector, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_symp_width_vector, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_symp_width_vector, log10(data.lr$NNQ))
plot(data.lr$d_symp_width_vector, data.lr$Abs_Benefit)
plot(data.lr$d_symp_width_vector, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$dispersion_vector, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$dispersion_vector, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$dispersion_vector, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$dispersion_vector, log10(data.lr$NNQ))
plot(data.lr$dispersion_vector, data.lr$Abs_Benefit)
plot(data.lr$dispersion_vector, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

summary(data.lr$R_0)
summary(data.lr$R_hsb)
summary(data.lr$R_s)
summary(data.lr$R_q)
summary(data.lr$Abs_Benefit)
summary(data.lr$NNQ)

summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 , "R_0"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "R_hsb"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "R_s"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "R_q"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "Abs_Benefit"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "NNQ"])

quantile(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "R_s"], c(0.025, 0.50, 0.975))
quantile(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6, "R_q"], c(0.025, 0.50, 0.975))

summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset_vector < 0, "R_0"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset_vector < 0, "R_hsb"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset_vector < 0, "R_s"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset_vector < 0, "R_q"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset_vector < 0, "Abs_Benefit"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset_vector < 0, "NNQ"])

#### Plot R_q and R_s ####

# Set range for relevant R_0 values
R_0_relevant.min <- 2.2
R_0_relevant.max <- 3.6

layout(c(1))
plot(data.hr[data.hr$R_0 > R_0_relevant.min & data.hr$R_0 < R_0_relevant.max, "R_s"], data.hr[data.hr$R_0 > R_0_relevant.min & data.hr$R_0 < R_0_relevant.max, "R_q"], col="lightblue", pch = 16, xlim = c(0, 4), ylim = c(0, 4), xlab = "R_s", ylab = "R_q")
points(data.lr$R_s, data.lr$R_q, xlim = c(0, 4), ylim = c(0, 4), col="blue", pch = 16)
points(data.lr$R_0, data.lr$R_0, xlim = c(0, 4), ylim = c(0, 4), col="white", pch = 1)
points(data.lr$R_0, data.lr$R_0, xlim = c(0, 4), ylim = c(0, 4), col="darkblue", pch = 20)

# color is R0, shape is HR or LR
# color bar along x and y for R0
# add lines X and Y

data.hr$Setting <- "HR"
data.lr$Setting <- "LR"
data.hr.lr <- rbind(data.hr, data.lr)

ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min & data.hr.lr$R_0 < R_0_relevant.max,]) +
  geom_vline(x=1, col="grey") + geom_hline(y=1, col="grey") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4, alpha = .1, fill = "yellow") + 
  annotate("rect", xmin = 1, xmax = 4, ymin = 0, ymax = 1, alpha = .1, fill = "blue") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = .1, fill = "green") +
  annotate("text", x = 2.5, y = 0.9, label = "Control with Quarantine", col = "blue") +
  annotate("text", x = 0.5, y = 3.5, label = "Control with\nSymptom Monitoring", col = "orange") +
  geom_point(aes(x=R_s, y=R_q, col = R_0, shape=Setting) ) +
  xlim(0,4) + ylim(0,4) +
  scale_colour_gradient(low="yellow", high="darkred") +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("High Resource", "Low Resource")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab("Effective Reproductive Number under Symptom Monitoring") +
  ylab("Effective Reproductive Number under Quarantine") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Disease: SARS")

ggplot(data = data.hr.lr[data.hr.lr$R_0 > R_0_relevant.min & data.hr.lr$R_0 < R_0_relevant.max,]) +
  annotate("rect", xmin = 2.2, xmax = 3.6, ymin = 0.5, ymax = 1, alpha = .1, fill = "green") +
  geom_hline(y=1, col = "grey") +
  geom_point(aes(x=R_0, y=R_s, shape = Setting), col = "darkgreen", alpha = 0.7) +
  geom_point(aes(x=R_0, y=R_q, shape = Setting), col = "blue", alpha = 0.7) +
  stat_smooth(aes(x=R_0, y=R_s, shape = Setting), method = "loess", color="darkgreen", size = 1.2) +
  stat_smooth(aes(x=R_0, y=R_q, shape = Setting), method = "loess", color="blue", size = 1.2) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(name = "Setting",
                     values = c(17, 1),
                     labels = c("High Resource", "Low Resource")) +
  guides(shape = guide_legend(reverse=TRUE)) +
  xlab(expression("Basic Reproductive Number R" [0])) + ylab(expression("Effective Reproductive Number R" [e])) +
  ggtitle("Disease: SARS")
