#### Header #### 
# SARS Case Study for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 25, 2015

#### Load Libraries ####
library(ppcor)
library(MASS)
library(lhs)
library(sensitivity)
library(ggplot2)
library(Hmisc)
library(reshape)

#### Step 1a: Choose range of parameters ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_12_20150824_5000reps.RData')
set.seed(1)

# Fixed Disease Parameters
parms_serial_interval <- list("weibull", 2, 10)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("weibull", 1.7, 6, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 5, 15, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Initialize
n_pop = 100
num_generations <- 5
times <- 5000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_lower","d_symp_upper","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -3
T_lat_offset.max <- 3
d_symp_lower.min <- 1
d_symp_lower.max <- 5
d_symp_upper.min <- 6
d_symp_upper.max <- 16
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_symp_lower = lhs[,2]*(d_symp_lower.max - d_symp_lower.min) + d_symp_lower.min,
  d_symp_upper = lhs[,3]*(d_symp_upper.max - d_symp_upper.min) + d_symp_upper.min,
  pi_t_triangle_center = lhs[,4]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_symp$parm1 <- as.numeric(params.set[i,"d_symp_lower"])
  parms_d_symp$parm2 <- as.numeric(params.set[i,"d_symp_upper"])
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  
  for (subseq_interventions in c(background_intervention)){      
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
                              parms_serial_interval)
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

data$T_lat_offset <- params.set[,"T_lat_offset"]
data$d_symp_lower <- params.set[,"d_symp_lower"]
data$d_symp_upper <- params.set[,"d_symp_upper"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval partial rank correlation
pcor(data[,c("ks","T_lat_offset","d_symp_lower","d_symp_upper","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_offset","d_symp_lower","d_symp_upper","pi_t_triangle_center")])$p.value[1,]

# Serial Interval Inspection Scatter
layout(cbind(c(1),c(2),c(3),c(4)))
plot(data$T_lat_offset, data$ks )
plot(data$d_symp_lower, data$ks)
plot(data$d_symp_upper, data$ks)
plot(data$pi_t_triangle_center, data$ks)

# Serial Interval Inspection Deciles
decile_plot_fcn(data, params.set)

# Narrow range in 10% increments
data <- data[data$T_lat_offset < sort(data$T_lat_offset)[round(nrow(data)*.9)],]
data <- data[data$T_lat_offset > sort(data$T_lat_offset)[round(nrow(data)*.1)],]
data <- data[data$d_symp_lower < sort(data$d_symp_lower)[round(nrow(data)*.9)],]
data <- data[data$d_symp_lower > sort(data$d_symp_lower)[round(nrow(data)*.1)],]
data <- data[data$d_symp_upper < sort(data$d_symp_upper)[round(nrow(data)*.9)],]
data <- data[data$d_symp_upper > sort(data$d_symp_upper)[round(nrow(data)*.1)],]
data <- data[data$pi_t_triangle_center < sort(data$pi_t_triangle_center)[round(nrow(data)*.9)],]
data <- data[data$pi_t_triangle_center > sort(data$pi_t_triangle_center)[round(nrow(data)*.1)],]

# Try to narrow the range of input parameters to improve serial interval fit
apply(data[,names(data.frame(params.set))], 2, min)
apply(data[,names(data.frame(params.set))], 2, max)

data_trim <- data

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_12_20150824_5000reps.RData")

#### Step 1b: Narrow down range of parameters ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_13_20150825_2500reps.RData')

set.seed(13)

# Fixed Disease Parameters
parms_serial_interval <- list("weibull", 2, 10)
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

parms_d_symp = list("uniform", 3, 8, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 100
num_generations <- 5
times <- 1000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_spread","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -1.25
T_lat_offset.max <- 0.5
d_symp_spread.min <- 2
d_symp_spread.max <- 6
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 0.5

params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_symp_spread = lhs[,2]*(d_symp_spread.max - d_symp_spread.min) + d_symp_spread.min,
  pi_t_triangle_center = lhs[,3]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_symp[c("parm1", "parm2")] <- c(as.numeric(7-params.set[i,"d_symp_spread"]), as.numeric(7+params.set[i,"d_symp_spread"]))
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  
  for (subseq_interventions in c(background_intervention)){      
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
                              parms_serial_interval)
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

data$T_lat_offset <- params.set[,"T_lat_offset"]
data$d_symp_spread <- params.set[,"d_symp_spread"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval partial rank correlation
pcor(data[,c("ks","T_lat_offset", "d_symp_spread","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_offset","d_symp_spread","pi_t_triangle_center")])$p.value[1,]

# Serial Interval Inspection Scatter
layout(cbind(c(1),c(2),c(3)))
plot(data$T_lat_offset, data$ks)
plot(data$d_symp_spread, data$ks)
plot(data$pi_t_triangle_center, data$ks)

# Serial Interval Inspection Deciles
decile_plot_fcn(data, params.set)

# Narrow range in 10% increments
data <- data[data$T_lat_offset < sort(data$T_lat_offset)[round(nrow(data)*.9)],]
data <- data[data$T_lat_offset > sort(data$T_lat_offset)[round(nrow(data)*.1)],]
data <- data[data$d_symp_spread < sort(data$d_symp_spread)[round(nrow(data)*.9)],]
data <- data[data$d_symp_spread > sort(data$d_symp_spread)[round(nrow(data)*.1)],]
data <- data[data$pi_t_triangle_center < sort(data$pi_t_triangle_center)[round(nrow(data)*.9)],]
data <- data[data$pi_t_triangle_center > sort(data$pi_t_triangle_center)[round(nrow(data)*.1)],]

# Try to narrow the range of input parameters to improve serial interval fit
apply(data[,names(data.frame(params.set))], 2, min)
apply(data[,names(data.frame(params.set))], 2, max)

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_13_20150825_2500reps.RData")

#### Ranking of Intervention Sensitivities ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_3_20150821.RData')

set.seed(3)

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 14, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 2000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0","pi_t_triangle_center", "T_lat_offset")

require(lhs)
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
pi_t_triangle_center.min <- 0.00
pi_t_triangle_center.max <- 0.50
T_lat_offset.min <- -0.5
T_lat_offset.max <- 0.5

params.set <- cbind(
  gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
  prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
  CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
  epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
  R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min,
  pi_t_triangle_center = lhs[,6]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min,
  T_lat_offset = lhs[,7]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]))
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 100
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 100
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 100
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
                              parms_serial_interval)
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

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_3_20150821.RData")

#### Case Study in High Resource Setting ####

set.seed(14)

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 3, 8, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0.9, 0.9, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 10
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("R_0","pi_t_triangle_center", "T_lat_offset")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

R_0.min <- 1
R_0.max <- 3
pi_t_triangle_center.min <- 0.00
pi_t_triangle_center.max <- 0.50
T_lat_offset.min <- -0.5
T_lat_offset.max <- 0.5

params.set <- cbind(
  R_0 = lhs[,1]*(R_0.max - R_0.min) + R_0.min,
  pi_t_triangle_center = lhs[,2]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min,
  T_lat_offset = lhs[,3]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]))
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 300
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 300
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 300
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
                              parms_serial_interval)
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
data$R_0_vector <- params.set[,"R_0"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data$T_lat_offset_vector <- params.set[,"T_lat_offset"]

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

summary(data$R_0)
summary(data$R_hsb)
summary(data$R_s)
summary(data$R_q)
summary(data$Abs_Benefit)
summary(data$NNQ)

summary(data[data$T_lat_offset_vector > 0, "R_0"])
summary(data[data$T_lat_offset_vector > 0, "R_hsb"])
summary(data[data$T_lat_offset_vector > 0, "R_s"])
summary(data[data$T_lat_offset_vector > 0, "R_q"])
summary(data[data$T_lat_offset_vector > 0, "Abs_Benefit"])
summary(data[data$T_lat_offset_vector > 0, "NNQ"])

summary(data[data$T_lat_offset_vector < 0, "R_0"])
summary(data[data$T_lat_offset_vector < 0, "R_hsb"])
summary(data[data$T_lat_offset_vector < 0, "R_s"])
summary(data[data$T_lat_offset_vector < 0, "R_q"])
summary(data[data$T_lat_offset_vector < 0, "Abs_Benefit"])
summary(data[data$T_lat_offset_vector < 0, "NNQ"])

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

data$resource_level <- NA
data[data$gamma_vector > 0.85 & data$gamma_vector < 0.95 &
       data$prob_CT_vector > 0.85 & data$prob_CT_vector < 0.95 &
       data$CT_delay_vector > 0.0 & data$CT_delay_vector < 0.1 &
       data$epsilon_vector > 0.0 & data$epsilon_vector < 0.5,
     "resource_level"] <- "high"
data[data$gamma_vector > 0.45 & data$gamma_vector < 0.55 &
       data$prob_CT_vector > 0.45 & data$prob_CT_vector < 0.55 &
       data$CT_delay_vector > 0.5 & data$CT_delay_vector < 1.5 &
       data$epsilon_vector > 1.5 & data$epsilon_vector < 2.5,
     "resource_level"] <- "low"

#### Case Study in Low Resource Setting ####

set.seed(15)

# Fixed Disease Parameters
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Variable Disease Parameters
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 3, 8, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

parms_epsilon = list("uniform", 2, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 10
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("R_0","pi_t_triangle_center", "T_lat_offset")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

R_0.min <- 1
R_0.max <- 3
pi_t_triangle_center.min <- 0.00
pi_t_triangle_center.max <- 0.50
T_lat_offset.min <- -0.5
T_lat_offset.max <- 0.5

params.set <- cbind(
  R_0 = lhs[,1]*(R_0.max - R_0.min) + R_0.min,
  pi_t_triangle_center = lhs[,2]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min,
  T_lat_offset = lhs[,3]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]))
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 300
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 300
    } else if (subseq_interventions == "s" | subseq_interventions == "q" & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
      n_pop_input <- 300
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
                              parms_serial_interval)
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
data$R_0_vector <- params.set[,"R_0"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data$T_lat_offset_vector <- params.set[,"T_lat_offset"]

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

summary(data$R_0)
summary(data$R_hsb)
summary(data$R_s)
summary(data$R_q)
summary(data$Abs_Benefit)
summary(data$NNQ)

summary(data[data$T_lat_offset_vector > 0, "R_0"])
summary(data[data$T_lat_offset_vector > 0, "R_hsb"])
summary(data[data$T_lat_offset_vector > 0, "R_s"])
summary(data[data$T_lat_offset_vector > 0, "R_q"])
summary(data[data$T_lat_offset_vector > 0, "Abs_Benefit"])
summary(data[data$T_lat_offset_vector > 0, "NNQ"])

summary(data[data$T_lat_offset_vector < 0, "R_0"])
summary(data[data$T_lat_offset_vector < 0, "R_hsb"])
summary(data[data$T_lat_offset_vector < 0, "R_s"])
summary(data[data$T_lat_offset_vector < 0, "R_q"])
summary(data[data$T_lat_offset_vector < 0, "Abs_Benefit"])
summary(data[data$T_lat_offset_vector < 0, "NNQ"])

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
