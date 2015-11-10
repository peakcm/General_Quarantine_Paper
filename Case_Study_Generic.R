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

#### Source Functions ####
source("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/Functions.R")

#### Load Workspaces ####
desired_root <- "20151024_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"

# If workspaces are in main folder
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_SMC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_PRCC.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_HR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_LR.RData", sep=""))
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "_Plots.RData", sep=""))

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
T_lat_offset.min <- 0
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
parms_serial_interval <- list("gamma", 2.5931, 0.1697) # Althaus 2015 Lancet ID
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

#### Disease: MERS ####

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "MERS"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("lognormal", 7.6, 1.77)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("lognormal", 5.2, 1.7, 999, "independent", "independent")
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
T_lat_offset.min <- -4
T_lat_offset.max <- 4
d_inf.min <- 5
d_inf.max <- 25
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

#### Disease: Pertussis ####

# Note: probably substantially slower because d_inf is so much longer, for each individual the rpois is drawn that many more times.

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "Pertussis"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.45585778, .11071164) # Approximation from te Beest reported quantiles
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("normal", 7, 1.53, 999, "independent", "independent") # (10-7)/1.96 = 1.53
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
T_lat_offset.min <- -10
T_lat_offset.max <- 4
d_inf.min <- 5
d_inf.max <- 14*7
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

#### Disease: Hepatitis A ####

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "HepatitisA"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 65.478798, 2.438040) # Approximation from Vink cf Brodribb
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 143.325515, 4.911873, 999, "independent", "independent") #Using data from Pickles 1930
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
T_lat_offset.min <- -14
T_lat_offset.max <- 4
d_inf.min <- 2
d_inf.max <- 20
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

#### Disease: Smallpox ####

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "Smallpox"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("lognormal", 15.53418, 1.253332) 
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("lognormal", 11.82245, 1.1853, 999, "independent", "independent") 
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
T_lat_offset.min <- -5
T_lat_offset.max <- 5
d_inf.min <- 2
d_inf.max <- 30
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

#### Disease: InfluenzaA ####

# Name the trial
date <- format(Sys.time(), "%Y%m%d")
disease <- "InfluenzaA"
root <- paste(date, disease, sep = "_")

# Fixed Disease Parameters
parms_serial_interval <- list("normal", 2.2, 0.8) 
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("lognormal", 1.4, 1.5, 999, "independent", "independent") # Using Vink 2014, cf Fine 2003
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
T_lat_offset.min <- -2
T_lat_offset.max <- 2
d_inf.min <- 2
d_inf.max <- 14
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
    if (i%%100 == 0){cat("\n")}
    
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
        cat("\nConvergence achieved in", SMC_counter, "iterations")
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

#### Case Study in High Resource Setting ####

# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 500
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
  R_0 = runif(n = times, min = 1.72, max = 1.94)) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
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
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
data.hr$R_0_input <- params.set[,"R_0"]
data.hr$T_lat_offset <- params.set[,"T_lat_offset"]
data.hr$dispersion <- params.set[,"dispersion"]

# Plot each of the covariate - outcome scatterplots
for (covariate in names(data.hr)[11:15]){
  panel_plot_fcn(data = data.hr, covariate = covariate)
  # cat ("Press [enter] to continue")
  # line <- readline()
}

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

# save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_HR.RData", sep=""))

#### Case Study in Middle Resource Setting ####

# Interventions
background_intervention <- "u"

prob_CT <- 0.75

gamma <- 0.75

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 500
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.mr <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.mr) <- names

# sample from joint posterior distribution
sample <- sample(x = row.names(data), size = times, replace = FALSE)
params.set <- cbind(
  T_lat_offset = data[sample, "T_lat_offset"],
  d_inf = data[sample, "d_inf"],
  pi_t_triangle_center = data[sample, "pi_t_triangle_center"],
  dispersion = runif(n=times, min = 1, max = 1),
  R_0 = runif(n = times, min = 1.72, max = 1.94)) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
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
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
      data.mr[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.mr[i,"ks"]  <- weighted.mean(x=In_Out$output[2:nrow(In_Out$output),"ks"], w=In_Out$output[2:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "hsb"){
      data.mr[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.mr[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.mr[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.mr[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
}

data.mr[,"Abs_Benefit"] <- data.mr[,"R_s"] - data.mr[,"R_q"]
data.mr[,"Rel_Benefit"] <- data.mr[,"Abs_Benefit"] / data.mr[,"R_s"]
data.mr[,"NNQ"] <- 1 / data.mr[,"Abs_Benefit"]
data.mr[data.mr$NNQ < 1,"NNQ"] <- 1
data.mr[data.mr$NNQ > 9999,"NNQ"] <- 9999
data.mr[data.mr$NNQ == Inf,"NNQ"] <- 9999
data.mr[,"Abs_Benefit_per_Qday"] <- data.mr[,"Abs_Benefit"] / data.mr[,"obs_to_iso_q"]
data.mr$d_inf <- params.set[,"d_inf"]
data.mr$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data.mr$R_0_input <- params.set[,"R_0"]
data.mr$T_lat_offset <- params.set[,"T_lat_offset"]
data.mr$dispersion <- params.set[,"dispersion"]

# Plot each of the covariate - outcome scatterplots
for (covariate in names(data.mr)[11:15]){
  panel_plot_fcn(data = data.mr, covariate = covariate)
  # cat ("Press [enter] to continue")
  # line <- readline()
}

summary(data.mr$R_0)
summary(data.mr$R_hsb)
summary(data.mr$R_s)
summary(data.mr$R_q)
summary(data.mr$Abs_Benefit)
summary(data.mr$NNQ)

summary(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6 , "R_0"])
summary(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "R_hsb"])
summary(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "R_s"])
summary(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "R_q"])
summary(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "Abs_Benefit"])
summary(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "NNQ"])

quantile(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "R_s"], c(0.025, 0.50, 0.975))
quantile(data.mr[data.mr$R_0 > 2.2 & data.mr$R_0 < 3.6, "R_q"], c(0.025, 0.50, 0.975))

# save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_MR.RData", sep=""))

#### Case Study in Low Resource Setting ####

# Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

parms_epsilon = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 500
num_generations <- 5
times <- 500
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
  R_0 = runif(n = times, min = 1.72, max = 1.94)) # note this is changed

for (i in 1:times){
  cat(".")
  if (i%%10 == 0){cat("|")}
  if (i%%100 == 0){cat("\n")}
  
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
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
data.lr$R_0_input <- params.set[,"R_0"]
data.lr$T_lat_offset <- params.set[,"T_lat_offset"]
data.lr$dispersion <- params.set[,"dispersion"]

# Plot each of the covariate - outcome scatterplots
for (covariate in names(data.lr)[11:15]){
  panel_plot_fcn(data = data.lr, covariate = covariate)
  # cat ("Press [enter] to continue")
  # line <- readline()
}

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

# save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_LR.RData", sep=""))

#### Save Workspace Image ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Plots.RData", sep=""))

