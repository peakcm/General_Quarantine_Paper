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

#### Load Workspaces ####
desired_root <- "20151021_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_PRCC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", v, "_Plots.RData", sep=""))

# If workspaces are in their own folder, named the same as the root
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_PRCC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_Plots.RData", sep=""))

#### Disease: SARS ####

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "SARS"
root <- paste(date, disease, sep = "_")

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

parms_d_inf = list("uniform", 1, 8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")


# Ranges for particle filter
T_lat_offset.min <- -7
T_lat_offset.max <- 2
d_inf.min <- 10
d_inf.max <- 40
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

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "Ebola"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent") # approximation from WHO
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent") 
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 8, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Ranges for particle filter
T_lat_offset.min <- 0
T_lat_offset.max <- 2
d_inf.min <- 10
d_inf.max <- 25
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

#### Particle Filter ####
setwd(dir = "~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/")

# Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Initialize
n_pop = 200
num_generations <- 4
times <- 1000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_inf","pi_t_triangle_center")
adaptive_thresh <- 0.80
SMC_times <- 15
perturb <- seq(from = 1/25, to = 1/50, length.out = SMC_times)
ks_conv_criteria <- 0.10 # convergence if median ks is within [ks_conv_criteria] percent of previous two 
ks_conv_stat <- rep(NA, SMC_times)
subseq_interventions <- "u"
printing = FALSE

# Run particle filter
SMC_break <- FALSE
SMC_counter <- 1
while (SMC_break == FALSE){
  cat('\nSMC iteration',SMC_counter, '\n')
  
  if (SMC_counter == 1){    
    lhs <- maximinLHS(times, length(dimensions))
    T_lat_offset.theta <- lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min
    d_inf.theta <- lhs[,2]*(d_inf.max - d_inf.min) + d_inf.min
    pi_t_triangle_center.theta <- lhs[,3]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min
  }
    
  for (i in 1:times){
    cat(".")
    if (i%%10 == 0){cat("|")}
    
    parms_T_lat$anchor_value <- T_lat_offset.theta[i]
    parms_d_inf$parm2 <- d_inf.theta[i]
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
      data[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
  }
  
  data$T_lat_offset <- T_lat_offset.theta
  data$d_inf <- d_inf.theta
  data$pi_t_triangle_center <- pi_t_triangle_center.theta
  data$weight <- (1/data$ks) / sum(1/data$ks)
  
  ks_conv_stat[SMC_counter] <- median(data$ks)
  
  # Save scatterplots
  pdf(paste(root, "_SMC_Iteration_",SMC_counter,".pdf", sep = ""))
  layout(rbind(c(1,2,3),c(4,5,6),c(7,8,9)))
  hist(data$T_lat_offset, xlim = c(T_lat_offset.min, T_lat_offset.max),
       main = "T_lat_offset", xlab = "days")
  plot(x=data$d_inf, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(T_lat_offset.min, T_lat_offset.max),
       xlim = c(d_inf.min, d_inf.max),
       main = paste("KS = ", round(ks_conv_stat[SMC_counter],3)), xlab = "days", ylab = "days")
  plot(x=data$pi_t_triangle_center, y=data$T_lat_offset, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(T_lat_offset.min, T_lat_offset.max),
       xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max),
       main = paste("Iteration ", SMC_counter), xlab = "proportion", ylab = "days")
  plot(x=data$T_lat_offset, y=data$d_inf, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(d_inf.min, d_inf.max),
       xlim = c(T_lat_offset.min, T_lat_offset.max), xlab = "days", ylab = "days")
  hist(data$d_inf, xlim = c(d_inf.min, d_inf.max),
       main = "d_inf", xlab = "days")  
  plot(x=data$pi_t_triangle_center, y=data$d_inf, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(d_inf.min, d_inf.max),
       xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max), xlab = "proportion", ylab = "days")
  plot(x=data$T_lat_offset, y=data$pi_t_triangle_center, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(pi_t_triangle_center.min, pi_t_triangle_center.max),
       xlim = c(T_lat_offset.min, T_lat_offset.max), xlab = "days", ylab = "proportion") 
  plot(x=data$d_inf, y=data$pi_t_triangle_center, col = rainbow(1000)[floor(data$ks*1000)+1], pch=16,
       ylim = c(pi_t_triangle_center.min, pi_t_triangle_center.max),
       xlim = c(d_inf.min, d_inf.max), xlab = "days", ylab = "proportion")
  hist(data$pi_t_triangle_center, xlim = c(pi_t_triangle_center.min, pi_t_triangle_center.max),
       main = "pi_t_triangle_center", xlab = "proportion")  
  dev.off()
  
  T_lat_offset.perturb <- (max(data[,"T_lat_offset"]) - min(data[,"T_lat_offset"])) * perturb[SMC_counter]
  d_inf.perturb <- (max(data[,"d_inf"]) - min(data[,"d_inf"])) * perturb[SMC_counter]
  pi_t_triangle_center.perturb <- (max(data[,"pi_t_triangle_center"]) - min(data[,"pi_t_triangle_center"])) * perturb[SMC_counter]
  
  #Adaptive threshold
  sorted.ks <- sort(data$ks)
  threshold <- sorted.ks[round(adaptive_thresh*length(sorted.ks))]
  
  #sample pre-candidate theta parameter sets from previous generation
  theta_pre_can <- sample(row.names(data[data$ks <= threshold,]), times, prob= data[data$ks <= threshold,"weight"], replace=TRUE)
  
  #perturb and propose
  T_lat_offset.theta <- sapply(data[theta_pre_can, "T_lat_offset"], function(x) rnorm(n=1, mean=x, sd=(T_lat_offset.max - T_lat_offset.min) * perturb[SMC_counter]))
  d_inf.theta <- sapply(data[theta_pre_can, "d_inf"], function(x) rnorm(n=1, mean=x, sd=(d_inf.max - d_inf.min) * perturb[SMC_counter]))
  pi_t_triangle_center.theta <- sapply(data[theta_pre_can, "pi_t_triangle_center"], function(x) rnorm(n=1, mean=x, sd=(pi_t_triangle_center.max - pi_t_triangle_center.min) * perturb[SMC_counter]))
  
  # Restrict range of candidates to the original range for the disease
  T_lat_offset.theta[T_lat_offset.theta < T_lat_offset.min] <- T_lat_offset.min
  T_lat_offset.theta[T_lat_offset.theta > T_lat_offset.max] <- T_lat_offset.max
  d_inf.theta[d_inf.theta < d_inf.min] <- d_inf.min
  d_inf.theta[d_inf.theta > d_inf.max] <- d_inf.max
  pi_t_triangle_center.theta[pi_t_triangle_center.theta < pi_t_triangle_center.min] <- pi_t_triangle_center.min
  pi_t_triangle_center.theta[pi_t_triangle_center.theta > pi_t_triangle_center.max] <- pi_t_triangle_center.max
  
  cat('\nMedian KS is ', ks_conv_stat[SMC_counter], "\n")
  
  # Check for "convergence"
  if (SMC_counter > 3){
    if (abs(ks_conv_stat[SMC_counter] - ks_conv_stat[SMC_counter-1])/ks_conv_stat[SMC_counter-1] <= ks_conv_criteria){
      if (abs(ks_conv_stat[SMC_counter] - ks_conv_stat[SMC_counter-2])/ks_conv_stat[SMC_counter-2] <= ks_conv_criteria){
        SMC_break <- TRUE
        cat("\nConvergence acheived in", SMC_counter, "iterations")
        break
      }
    }
  }
  
  if (SMC_counter >= SMC_times){
    SMC_break <- TRUE
    cat("\nUnable to converge by", SMC_counter, "SMC iterations")
    break
  }
  
  SMC_counter <- SMC_counter + 1
  
}

# Save tables and output files
write.table(ks_conv_stat, paste(root,"ks_conv_stat.csv", sep="_"))
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_SMC.RData", sep=""))

#### PRCC Ranking of Intervention Sensitivities ####

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
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.prcc <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.prcc) <- names

# Sample from posterior distributions for each parameter independently
params.set <- cbind(
  T_lat_offset = sample(data$T_lat_offset, size = times, replace = TRUE),
  d_inf = sample(data$d_inf, size = times, replace = TRUE),
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
                    dispersion = lhs[,6]*(dispersion.max - dispersion.min) + dispersion.min,
                    R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min)
params.set <- data.frame(params.set)

i=1
while (i <= times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_inf$parm2 <- as.numeric(params.set[i,"d_inf"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]))
  dispersion <- as.numeric(params.set[i, "dispersion"])
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
data.prcc[,"Rel_Benefit"] <- data.prcc[,"Abs_Benefit"] / data.prcc[,"R_s"]
data.prcc[,"NNQ"] <- 1 / data.prcc[,"Abs_Benefit"]
data.prcc[data.prcc$NNQ < 1,"NNQ"] <- 1
data.prcc[data.prcc$NNQ > 9999,"NNQ"] <- 9999
data.prcc[data.prcc$NNQ == Inf,"NNQ"] <- 9999
data.prcc[,"Abs_Benefit_per_Qday"] <- data.prcc[,"Abs_Benefit"] / data.prcc[,"obs_to_iso_q"]
data.prcc$gamma <- params.set[,"gamma"]
data.prcc$prob_CT <- params.set[,"prob_CT"]
data.prcc$CT_delay <- params.set[,"CT_delay"]
data.prcc$epsilon <- params.set[,"epsilon"]
data.prcc$R_0 <- params.set[,"R_0"]
data.prcc$dispersion <- params.set[,"dispersion"]
data.prcc$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.prcc$T_lat_offset <- params.set[,"T_lat_offset"]
data.prcc$d_inf <- params.set[,"d_inf"]

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$gamma, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$gamma, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$gamma, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$gamma, log10(data.prcc$NNQ))
plot(data.prcc$gamma, data.prcc$Abs_Benefit)
plot(data.prcc$gamma, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$prob_CT, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$prob_CT, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$prob_CT, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$prob_CT, log10(data.prcc$NNQ))
plot(data.prcc$prob_CT, data.prcc$Abs_Benefit)
plot(data.prcc$prob_CT, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$CT_delay, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$CT_delay, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$CT_delay, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$CT_delay, log10(data.prcc$NNQ))
plot(data.prcc$CT_delay, data.prcc$Abs_Benefit)
plot(data.prcc$CT_delay, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$epsilon, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$epsilon, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$epsilon, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$epsilon, log10(data.prcc$NNQ))
plot(data.prcc$epsilon, data.prcc$Abs_Benefit)
plot(data.prcc$epsilon, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$R_0, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$R_0, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$R_0, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$R_0, log10(data.prcc$NNQ))
plot(data.prcc$R_0, data.prcc$Abs_Benefit)
plot(data.prcc$R_0, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$pi_t_triangle_center, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$pi_t_triangle_center, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$pi_t_triangle_center, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$pi_t_triangle_center, log10(data.prcc$NNQ))
plot(data.prcc$pi_t_triangle_center, data.prcc$Abs_Benefit)
plot(data.prcc$pi_t_triangle_center, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$T_lat_offset, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$T_lat_offset, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$T_lat_offset, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$T_lat_offset, log10(data.prcc$NNQ))
plot(data.prcc$T_lat_offset, data.prcc$Abs_Benefit)
plot(data.prcc$T_lat_offset, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.prcc$dispersion, data.prcc$R_0, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$dispersion, data.prcc$R_s, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$dispersion, data.prcc$R_q, ylim=c(0, max(data.prcc$R_0)))
plot(data.prcc$dispersion, log10(data.prcc$NNQ))
plot(data.prcc$dispersion, data.prcc$Abs_Benefit)
plot(data.prcc$dispersion, data.prcc$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(t(c(1,2)))
plot(data.prcc$pi_t_triangle_center, data.prcc$ks)
plot(data.prcc$T_lat_offset, data.prcc$ks)

decile_plot_fcn(data.prcc, params.set)

data.prcc_store <- data.prcc

# Partial Rank Correlation using "ppcor" package
require(ppcor)
pcor(data.prcc[,c("R_s",names(data.prcc)[11:length(names(data.prcc))])], method=c("spearman"))$estimate[1,]
pcor(data.prcc[,c("R_s",names(data.prcc)[11:length(names(data.prcc))])], method=c("spearman"))$p.value[1,]
pcor(data.prcc[,c("Abs_Benefit",names(data.prcc)[11:length(names(data.prcc))])], method=c("spearman"))$estimate[1,]
pcor(data.prcc[,c("Abs_Benefit",names(data.prcc)[11:length(names(data.prcc))])], method=c("spearman"))$p.value[1,]
pcor(data.prcc[,c("Abs_Benefit_per_Qday", names(data.prcc)[11:length(names(data.prcc))])], method=c("spearman"))$estimate[1,]
pcor(data.prcc[,c("Abs_Benefit_per_Qday", names(data.prcc)[11:length(names(data.prcc))])], method=c("spearman"))$p.value[1,]

# Partial Rank Correlation using "sensitivity" package
require(sensitivity)
bonferroni.alpha <- 0.05/length(dimensions)
data.prcc_sensitivity <- data.frame(matrix(rep(NA, 7*length(names)*(ncol(params.set)-1)), ncol=7)) 
names(data.prcc_sensitivity) <- c("output","parameter","coef","bias","stderr","CImin","CImax")
data.prcc_sensitivity$output <- rep(names, each = (ncol(params.set)-1))
data.prcc_sensitivity$parameter <- rep(names(params.set)[-ncol(params.set)], times = length(names))
for (output in names){
  prcc <- pcc(data.prcc[,(length(names)+1):length(names(data.prcc))], data.prcc[,output], nboot = 100, rank=TRUE, conf=1-bonferroni.alpha)
  summary <- print(prcc)
  data.prcc_sensitivity[data.prcc_sensitivity$output == output,3:7] <- summary
}

require(ggplot2)
ggplot(prcc_data.prcc, aes(x = parameter, y= coef)) +
  facet_grid(output ~ .) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() +
  geom_errorbar(data.prcc = prcc_data.prcc, aes(ymin = CImin, ymax = CImax), width = 0.1)

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/", root, "_PRCC.RData", sep=""))

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
times <- 50
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.hr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.hr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1, max = 10))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  
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
                              dispersion = dispersion,
                              printing = printing)
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
data.hr$d_inf <- params.set[,"d_inf"]
data.hr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.hr$R_0 <- params.set[,"R_0"]
data.hr$T_lat_offset <- params.set[,"T_lat_offset"]
data.hr$dispersion <- params.set[,"dispersion"]

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$R_0, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$R_0, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$R_0, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$R_0, log10(data.hr$NNQ))
plot(data.hr$R_0, data.hr$Abs_Benefit)
plot(data.hr$R_0, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$pi_t_triangle_center, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$pi_t_triangle_center, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$pi_t_triangle_center, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$pi_t_triangle_center, log10(data.hr$NNQ))
plot(data.hr$pi_t_triangle_center, data.hr$Abs_Benefit)
plot(data.hr$pi_t_triangle_center, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$T_lat_offset, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$T_lat_offset, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$T_lat_offset, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$T_lat_offset, log10(data.hr$NNQ))
plot(data.hr$T_lat_offset, data.hr$Abs_Benefit)
plot(data.hr$T_lat_offset, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.hr$d_inf, data.hr$R_0, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_inf, data.hr$R_s, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_inf, data.hr$R_q, ylim=c(0, max(data.hr$R_0)))
plot(data.hr$d_inf, log10(data.hr$NNQ))
plot(data.hr$d_inf, data.hr$Abs_Benefit)
plot(data.hr$d_inf, data.hr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

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

summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset > 0, "R_0"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset > 0, "R_hsb"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset > 0, "R_s"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset > 0, "R_q"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset > 0, "Abs_Benefit"])
summary(data.hr[data.hr$R_0 > 2.2 & data.hr$R_0 < 3.6 & data.hr$T_lat_offset > 0, "NNQ"])

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/", root, "_HR.RData", sep=""))

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
times <- 50
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.lr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.lr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1, max = 10))

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  
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
                              dispersion = dispersion,
                              printing = printing)
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
data.lr$d_inf <- params.set[,"d_inf"]
data.lr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.lr$R_0 <- params.set[,"R_0"]
data.lr$T_lat_offset <- params.set[,"T_lat_offset"]
data.lr$dispersion <- params.set[,"dispersion"]

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$R_0, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$R_0, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$R_0, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$R_0, log10(data.lr$NNQ))
plot(data.lr$R_0, data.lr$Abs_Benefit)
plot(data.lr$R_0, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$pi_t_triangle_center, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$pi_t_triangle_center, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$pi_t_triangle_center, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$pi_t_triangle_center, log10(data.lr$NNQ))
plot(data.lr$pi_t_triangle_center, data.lr$Abs_Benefit)
plot(data.lr$pi_t_triangle_center, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$T_lat_offset, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$T_lat_offset, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$T_lat_offset, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$T_lat_offset, log10(data.lr$NNQ))
plot(data.lr$T_lat_offset, data.lr$Abs_Benefit)
plot(data.lr$T_lat_offset, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$d_inf, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_inf, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_inf, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$d_inf, log10(data.lr$NNQ))
plot(data.lr$v, data.lr$Abs_Benefit)
plot(data.lr$d_inf, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data.lr$dispersion, data.lr$R_0, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$dispersion, data.lr$R_s, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$dispersion, data.lr$R_q, ylim=c(0, max(data.lr$R_0)))
plot(data.lr$dispersion, log10(data.lr$NNQ))
plot(data.lr$dispersion, data.lr$Abs_Benefit)
plot(data.lr$dispersion, data.lr$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

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

summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset > 0, "R_0"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset > 0, "R_hsb"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset > 0, "R_s"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset > 0, "R_q"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset > 0, "Abs_Benefit"])
summary(data.lr[data.lr$R_0 > 2.2 & data.lr$R_0 < 3.6 & data.lr$T_lat_offset > 0, "NNQ"])

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/", root, "_LR.RData", sep=""))

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
  ggtitle(paste("Disease: ", disease))

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
  ggtitle(paste("Disease: ", disease))

save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/", root, "_Plots.RData", sep=""))

