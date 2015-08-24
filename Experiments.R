#### Header #### 
# Experiments for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 19, 2015

#### Load Libraries ####
library(ppcor)
library(MASS)
library(lhs)
library(sensitivity)
library(ggplot2)
library(Hmisc)

#### Experiment 3: Ebola Case Study ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_3_20150821.RData')
# Need to redo

set.seed(3)

background_intervention <- "u"

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# Ebola serial interval approximation
parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# WHO NEJM
# Liberia R0 1.83 (1.72, 1.94)
parms_R_0 = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# All countries: Incubation period
parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
# Set T_lat = T_inc
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dixon 2015 EID
prob_CT <- 0.3 #27.4 to 31.1% of case patients were in the contact registry before identification. around 30% of new cases reported having any contacts they could have infected 

# Less sure
parms_d_symp = list("triangle", 1, 14, 7, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 11, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
# Set d_symp = d_inf

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 50
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
pi_t_triangle_center.max <- 0.30
T_lat_offset.min <- -3
T_lat_offset.max <- 3

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

data_store <- data
data <- data[data$pi_t_triangle_center > 0.2 & data$pi_t_triangle_center < 0.5,]

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

save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_3_20150821.RData")

#### Experiment 4: Does variance in parms_R_0 matter? #### 
# Compare output abs_benefit, R_s, etc. when the R_0 is set to the mean for all people
#   versus when it has a range of up to 0 to 2*R_0
# If we find that the output doesn't depend on the R_0 variance, then we can set all people in a population to the mean R_0.
#   There is some variability in the R_0 within a population run due to the poisson distribution of drawing each persons number of cases
#   You can add variability in the desired R_0 via lhs 
#   Therefore, it doesn't matter if your child has a different R_0 draw from you
# This will help with inferring the R_0_hsb that is desired under different background interventions

# Conclusions: The variance of the R_0 input for a population does not matter. Only the mean R_0 matters.   

# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_4_28072015.RData')

set.seed(2346342)

background_intervention <- "u"

parms_pi_t <- list("uniform", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_R_0 = list("uniform", 1, 5, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_T_inc = list("uniform", 1, 20, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_lat = list("uniform", 1, 20, 999, "independent", "independent")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

prob_CT <- 0.5
gamma <- 0.8

parms_d_symp = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 1000
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("R_0.mean","R_0.var")

lhs <- maximinLHS(times, length(dimensions))

R_0.mean.min <- 0.5
R_0.mean.max <- 5
R_0.var.min  <- 0
R_0.var.max  <- 1

params.set <- cbind(
  R_0.mean = lhs[,1]*(R_0.mean.max - R_0.mean.min) + R_0.mean.min,
  R_0.var  = lhs[,2]*(R_0.var.max - R_0.var.min) + R_0.var.min)
params.set[,"R_0.var"] <- params.set[,"R_0.mean"] * params.set[,"R_0.var"]
params.set[params.set[,"R_0.var"]<0] <- 0
params.set <- params.set[order(params.set[,"R_0.mean"]),]
params.set[seq(1,times) %% 2 == 0,"R_0.var"] <- 0  #seq variance to zero for the even numbered rows

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  # Spread of R_0 is always 0.5. The mean is changing
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0.mean"]) - as.numeric(params.set[i,"R_0.var"]) , as.numeric(params.set[i,"R_0.mean"]) + as.numeric(params.set[i,"R_0.var"]))
  
  if (params.set[i,"R_0.mean"] - params.set[i,"R_0.var"] < 0) { parms_R_0[c("parm1")] <- c(0)}
  
  for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
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

data[,"Abs_Benefit"] <- data[,"R_s"] - data[,"R_q"]
data[,"Rel_Benefit"] <- data[,"Abs_Benefit"] / data[,"R_s"]
data[,"NNQ"] <- 1 / data[,"Abs_Benefit"]
data[data$NNQ < 1,"NNQ"] <- 1
data[data$NNQ == Inf,"NNQ"] <- 9999
data[data$NNQ > 9999,"NNQ"] <- 9999
data[,"Abs_Benefit_per_Qday"] <- data[,"Abs_Benefit"] / data[,"obs_to_iso_q"]
data$R_0_mean_vector <- params.set[,"R_0.mean"]
data$R_0_var_vector <- params.set[,"R_0.var"]

data$group <- 1
data[data$R_0_var_vector == 0,"group"] <- 0
data$group <- as.factor(data$group)

ggplot(data, aes(x=R_0_mean_vector, group=group)) +
  geom_point(aes(y=R_0, shape=group, color=group), size=2) + theme_bw()
ggplot(data, aes(x=R_0_mean_vector, group=group)) +
  geom_point(aes(y=R_s, shape=group, color=group), size=2) + theme_bw()
ggplot(data, aes(x=R_0_mean_vector, group=group)) +
  geom_point(aes(y=R_q, shape=group, color=group), size=2) + theme_bw()
ggplot(data, aes(x=R_0_mean_vector, group=group)) +
  geom_point(aes(y=Abs_Benefit, shape=group, color=group), size=2) + theme_bw()

summary(data[data$group == 1,"R_0"])
summary(data[data$group == 0,"R_0"])
t.test(data$R_0 ~data$group)

summary(data[data$group == 1,"R_s"])
summary(data[data$group == 0,"R_s"])
t.test(data$R_s ~data$group)

summary(data[data$group == 1,"Abs_Benefit"])
summary(data[data$group == 0,"Abs_Benefit"])
t.test(data$Abs_Benefit ~data$group)

# Partial Rank Correlation using "ppcor" package
pcor(data[,c("Abs_Benefit",names(data)[10:(length(names(data))-1)])], method=c("spearman"))$estimate[1,]
pcor(data[,c("Abs_Benefit",names(data)[10:(length(names(data))-1)])], method=c("spearman"))$p.value[1,]
pcor(data[,c("Abs_Benefit_per_Qday", names(data)[10:(length(names(data))-1)])], method=c("spearman"))$estimate[1,]
pcor(data[,c("Abs_Benefit_per_Qday", names(data)[10:(length(names(data))-1)])], method=c("spearman"))$p.value[1,]
#after controlling for the mean, does the variance matter? Check p-value

# save.image('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_4_28072015.RData')

#### Experiment 5: SARS Case Study #### 
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_5_20150807_sars1000reps.RData')

# Need to redo

set.seed(5)

background_intervention = "u"

# SARS serial interval approximation
parms_serial_interval <- list("weibull", 2, 10)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# [Lipsitch 2003 Science]
# Extensive contact tracing in Hong Kong has failed to identify a known symptomatic SARS contact for 8.6% of reported cases (14)
Prob_CT = 0.914
# quarantine is assumed to be 100% effective for those contacts who are found before they become infectious.
gamma = 1
# Best estimates for R_0 are from 2.2 to 3.6
parms_R_0 = list("uniform", 2.2, 3.6, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
# Infectiousness increases until day 10 then decreases
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# [Donnelly 2003 Lancet]
# mean incubation period of the disease is estimated to be 6·4 days (95% CI 5·2–7·7).
# eye-balling figure 2a
parms_T_inc = list("weibull", 1.7, 6, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Unsure
parms_T_lat = list("gamma", 1.75, 0.182, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 5, 20, 10, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("confit", 25, 10, 999, -3, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 200
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0","pi_t_triangle_center", "T_lat_offset", "d_inf_offset")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

# Refer to Experiment 6b to tune these ranges
gamma.min <- 0.5
gamma.max <- 1
prob_CT.min <- 0.25
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 5
epsilon.min <- 0
epsilon.max <- 7
R_0.min <- 1
R_0.max <- 5
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 0.4
T_lat_offset.min <- -1.5
T_lat_offset.max <- 0
d_inf_offset.min <- -3
d_inf_offset.max <- 1

params.set <- cbind(
  gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
  prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
  CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
  epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
  R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min,
  pi_t_triangle_center = lhs[,6]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min,
  T_lat_offset = lhs[,7]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_inf_offset = lhs[,8]*(d_inf_offset.max - d_inf_offset.min) + d_inf_offset.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_inf$anchor_value <- as.numeric(params.set[i,"d_inf_offset"])
  
  #R_0 has no variance within trials
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set[i,"R_0"]), as.numeric(params.set[i,"R_0"]))
  
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

data$gamma_vector <- params.set[,"gamma"]
data$prob_CT_vector <- params.set[,"prob_CT"]
data$CT_delay_vector <- params.set[,"CT_delay"]
data$epsilon_vector <- params.set[,"epsilon"]
data$R_0_vector <- params.set[,"R_0"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data$T_lat_offset <- params.set[,"T_lat_offset"]
data$d_inf_offset <- params.set[,"d_inf_offset"]

data_store <- data
data <- data[is.na(data$R_s) == 0,] # There are 11 observations with NA for R_s oddly. Look into why this happened for R_s but not R_q

data[,"Abs_Benefit"] <- data[,"R_s"] - data[,"R_q"]
data[,"Rel_Benefit"] <- data[,"Abs_Benefit"] / data[,"R_s"]
data[,"NNQ"] <- 1 / data[,"Abs_Benefit"]
data[data$NNQ < 1,"NNQ"] <- 1
data[data$NNQ > 9999,"NNQ"] <- 9999
data[data$NNQ == Inf,"NNQ"] <- 9999
data[,"Abs_Benefit_per_Qday"] <- data[,"Abs_Benefit"] / data[,"obs_to_iso_q"]

data <- data[data$pi_t_triangle_center > 0.3 & data$pi_t_triangle_center < 0.5,]

# Serial Interval Inspection
layout(cbind(c(1), c(2), c(3)))
plot(data$pi_t_triangle_center, data$ks)
plot(data$T_lat_offset, data$ks)
plot(data$d_inf_offset, data$ks)
pcor(data[,c("ks","pi_t_triangle_center","T_lat_offset","d_inf_offset")], method=c("spearman"))$estimate[1,]
pcor(data[,c("ks","pi_t_triangle_center","T_lat_offset","d_inf_offset")], method=c("spearman"))$p.value[1,]

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
plot(data$T_lat_offset, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$T_lat_offset, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$T_lat_offset, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$T_lat_offset, log10(data$NNQ))
plot(data$T_lat_offset, data$Abs_Benefit)
plot(data$T_lat_offset, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

layout(cbind(c(1,2,3),c(4,5,6)))
plot(data$d_inf_offset, data$R_0, ylim=c(0, max(data$R_0)))
plot(data$d_inf_offset, data$R_s, ylim=c(0, max(data$R_0)))
plot(data$d_inf_offset, data$R_q, ylim=c(0, max(data$R_0)))
plot(data$d_inf_offset, log10(data$NNQ))
plot(data$d_inf_offset, data$Abs_Benefit)
plot(data$d_inf_offset, data$Abs_Benefit_per_Qday, ylim=c(-0.4, 0.4))

# Partial Rank Correlation using "ppcor" package
require(ppcor)
pcor(data[,c("R_s",names(data)[10:length(names(data))])], method=c("spearman"))$estimate[1,]
pcor(data[,c("R_s",names(data)[10:length(names(data))])], method=c("spearman"))$p.value[1,]
pcor(data[,c("Abs_Benefit",names(data)[10:length(names(data))])], method=c("spearman"))$estimate[1,]
pcor(data[,c("Abs_Benefit",names(data)[10:length(names(data))])], method=c("spearman"))$p.value[1,]
pcor(data[,c("Abs_Benefit_per_Qday", names(data)[10:length(names(data))])], method=c("spearman"))$estimate[1,]
pcor(data[,c("Abs_Benefit_per_Qday", names(data)[10:length(names(data))])], method=c("spearman"))$p.value[1,]

# Partial Rank Correlation using "sensitivity" package
require(sensitivity)
bonferroni.alpha <- 0.05/length(dimensions)
prcc_data <- data.frame(matrix(rep(NA, 7*length(names)*length(dimensions)), ncol=7)) 
names(prcc_data) <- c("output","parameter","coef","bias","stderr","CImin","CImax")
prcc_data$output <- rep(names, each = length(dimensions))
prcc_data$parameter <- rep(dimensions, times = length(names))
for (output in names){
  prcc <- pcc(data[,(length(names)+1):length(names(data))], data[,output], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
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

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_5_20150807_sars1000reps.RData")

#### Experiment 6: Optimize serial interval for SARS Case study #### 
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_6_20150807_4000reps.RData')

set.seed(3456)

background_intervention = "u"

# SARS serial interval approximation
parms_serial_interval <- list("weibull", 2, 10)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# [Lipsitch 2003 Science]
# Extensive contact tracing in Hong Kong has failed to identify a known symptomatic SARS contact for 8.6% of reported cases (14)
prob_CT = 0.914
# quarantine is assumed to be 100% effective for those contacts who are found before they become infectious.
gamma = 1
# Best estimates for R_0 are from 2.2 to 3.6
parms_R_0 = list("uniform", 2.2, 3.6, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
# Infectiousness increases until day 10 then decreases
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# [Donnelly 2003 Lancet]
# mean incubation period of the disease is estimated to be 6·4 days (95% CI 5·2–7·7).
# eye-balling figure 2a
parms_T_inc = list("uniform", 4, 8, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Unsure
parms_T_lat = list("uniform", 4, 8, 999, "independent", "independent")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 5, 15, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 5, 15, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 200
num_generations <- 5
times <- 4000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_lower","T_lat_upper", "T_inc_lower","T_inc_upper","d_inf_lower","d_inf_upper","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_lower.min <- 1
T_lat_lower.max <- 4
T_lat_upper.min <- 5
T_lat_upper.max <- 15
T_inc_lower.min <- 1
T_inc_lower.max <- 6
T_inc_upper.min <- 6.5
T_inc_upper.max <- 15
d_inf_lower.min <- 2
d_inf_lower.max <- 9.5
d_inf_upper.min <- 10
d_inf_upper.max <- 20
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 0.25


params.set <- cbind(
  T_lat_lower = lhs[,1]*(T_lat_lower.max - T_lat_lower.min) + T_lat_lower.min,
  T_lat_upper = lhs[,2]*(T_lat_upper.max - T_lat_upper.min) + T_lat_upper.min,
  T_inc_lower = lhs[,3]*(T_inc_lower.max - T_inc_lower.min) + T_inc_lower.min,
  T_inc_upper = lhs[,4]*(T_inc_upper.max - T_inc_upper.min) + T_inc_upper.min,
  d_inf_lower = lhs[,5]*(d_inf_lower.max - d_inf_lower.min) + d_inf_lower.min,
  d_inf_upper = lhs[,6]*(d_inf_upper.max - d_inf_upper.min) + d_inf_upper.min,
  pi_t_triangle_center = lhs[,7]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$parm1 <- as.numeric(params.set[i,"T_lat_lower"])
  parms_T_lat$parm2 <- as.numeric(params.set[i,"T_lat_upper"])
  parms_T_inc$parm1 <- as.numeric(params.set[i,"T_inc_lower"])
  parms_T_inc$parm2 <- as.numeric(params.set[i,"T_inc_upper"])
  parms_d_inf$parm1 <- as.numeric(params.set[i,"d_inf_lower"])
  parms_d_inf$parm2 <- as.numeric(params.set[i,"d_inf_upper"])
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

data$T_lat_lower <- params.set[,"T_lat_lower"]
data$T_lat_upper <- params.set[,"T_lat_upper"]
data$T_inc_lower <- params.set[,"T_inc_lower"]
data$T_inc_upper <- params.set[,"T_inc_upper"]
data$d_inf_lower <- params.set[,"d_inf_lower"]
data$d_inf_upper <- params.set[,"d_inf_upper"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval Inspection
layout(cbind(c(1),c(2),c(3),c(4),c(5), c(6), c(7)))
plot(data$T_lat_lower, data$ks)
plot(data$T_lat_upper, data$ks)
plot(data$T_inc_lower, data$ks)
plot(data$T_inc_upper, data$ks)
plot(data$d_inf_lower, data$ks)
plot(data$d_inf_upper, data$ks)
plot(data$pi_t_triangle_center, data$ks)
pcor(data[,c("ks","T_lat_lower","T_lat_upper","T_inc_lower","T_inc_upper","d_inf_lower","d_inf_upper","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_lower","T_lat_upper","T_inc_lower","T_inc_upper","d_inf_lower","d_inf_upper","pi_t_triangle_center")])$p.value[1,]

# Try to narrow the range of input parameters to improve serial interval fit
data <- data[data$T_lat_upper < 10,] 
data <- data[data$T_lat_upper > 6,] 
data <- data[data$T_lat_upper < 8.5,]
data <- data[data$T_inc_upper > 8,]
data <- data[data$T_inc_lower < 5.5,]
data <- data[data$T_lat_lower > 1.25,]
data <- data[data$pi_t_triangle_center < 0.2,]

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_6_20150807_4000reps.RData")

#### Experiment 6b: Optimize serial interval for SARS Case study, but anchor T_lat and d_inf to T_inc and d_symp ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_6b_20150807_4000reps.RData')

set.seed(45464)

background_intervention = "u"

# SARS serial interval approximation
parms_serial_interval <- list("weibull", 2, 10)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# [Lipsitch 2003 Science]
# Extensive contact tracing in Hong Kong has failed to identify a known symptomatic SARS contact for 8.6% of reported cases (14)
Prob_CT = 0.914
# quarantine is assumed to be 100% effective for those contacts who are found before they become infectious.
gamma = 1
# Best estimates for R_0 are from 2.2 to 3.6
parms_R_0 = list("uniform", 2.2, 3.6, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")
# Infectiousness increases until day 10 then decreases
parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# [Donnelly 2003 Lancet]
# mean incubation period of the disease is estimated to be 6·4 days (95% CI 5·2–7·7).
# eye-balling figure 2a
parms_T_inc = list("weibull", 1.7, 6, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Anchor T_lat to T_inc
parms_T_lat = list("gamma", 1.75, 0.182, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 5, 16, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Anchor d_inf to d_symp
parms_d_inf = list("confit", 25, 10, 999, -3, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 100
num_generations <- 5
times <- 4000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_lower","d_symp_upper","d_inf_offset","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -3
T_lat_offset.max <- 2
d_symp_lower.min <- 3.1
d_symp_lower.max <- 9.5
d_symp_upper.min <- 10
d_symp_upper.max <- 20
d_inf_offset.min <- -3
d_inf_offset.max <- 1
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 0.5

params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_symp_lower = lhs[,2]*(d_symp_lower.max - d_symp_lower.min) + d_symp_lower.min,
  d_symp_upper = lhs[,3]*(d_symp_upper.max - d_symp_upper.min) + d_symp_upper.min,
  d_inf_offset = lhs[,4]*(d_inf_offset.max - d_inf_offset.min) + d_inf_offset.min,
  pi_t_triangle_center = lhs[,5]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_symp$parm1 <- as.numeric(params.set[i,"d_symp_lower"])
  parms_d_symp$parm2 <- as.numeric(params.set[i,"d_symp_upper"])
  parms_d_inf$anchor_value <- as.numeric(params.set[i,"d_inf_offset"])
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
data$d_inf_offset <- params.set[,"d_inf_offset"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval Inspection
layout(cbind(c(1),c(2),c(3),c(4),c(5)))
plot(data$T_lat_offset, data$ks)
plot(data$d_symp_lower, data$ks)
plot(data$d_symp_upper, data$ks)
plot(data$d_inf_offset, data$ks)
plot(data$pi_t_triangle_center, data$ks)
pcor(data[,c("ks","T_lat_offset","d_symp_lower","d_symp_upper","d_inf_offset","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_offset","d_symp_lower","d_symp_upper","d_inf_offset","pi_t_triangle_center")])$p.value[1,]

data <- data[data$T_lat_offset > -2 & data$T_lat_offset < 0,]
data <- data[data$d_symp_upper > 12,]
data <- data[data$T_lat_offset > -1.5,]
data <- data[data$pi_t_triangle_center < 0.4,]
data <- data[data$d_symp_lower > 4,]

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_6b_20150807_4000reps.RData")

#### Experiment 7: Does the prcc rank of a parameter change if the range of that parameter changes? #### 
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_7_20150806.RData')

set.seed(5647)

background_intervention <- "u"

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_R_0 = list("uniform", 0.6, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_inc = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_lat = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

prob_CT <- 0.3 #27.4 to 31.1% of case patients were in the contact registry before identification. around 30% of new cases reported having any contacts they could have infected 

parms_d_inf = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 10, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 2000
names <- c("R_0", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0_input","pi_t_triangle_center", "T_inc", "T_lat", "d_inf", "d_symp")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0.1
gamma.max <- 0.9
prob_CT.min <- 0.1
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 15
epsilon.min <- 0
epsilon.max <- 15
R_0.min <- 1
R_0.max <- 5
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1
T_inc.min <- 2
T_inc.max <- 20
T_lat.min <- 2
T_lat.max <- 20
d_inf.min <- 2
d_inf.max <- 20
d_symp.min <- 2
d_symp.max <- 20

params.set <- cbind(
  gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
  prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
  CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
  epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
  R_0_input = lhs[,5]*(R_0.max - R_0.min) + R_0.min,
  pi_t_triangle_center = lhs[,6]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min,
  T_inc = lhs[,7]*(T_inc.max - T_inc.min) + T_inc.min,
  T_lat = lhs[,8]*(T_lat.max - T_lat.min) + T_lat.min,
  d_inf = lhs[,9]*(d_inf.max - d_inf.min) + d_inf.min,
  d_symp = lhs[,10]*(d_symp.max - d_symp.min) + d_symp.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  parms_R_0$parm2 <- as.numeric(params.set[i,"R_0_input"])
  parms_pi_t$triangle_center <- as.numeric(params.set[i,"pi_t_triangle_center"])
  parms_T_inc$parm2 <- as.numeric(params.set[i,"T_inc"])
  parms_T_lat$parm2 <- as.numeric(params.set[i,"T_lat"])
  parms_d_inf$parm2 <- as.numeric(params.set[i,"d_inf"])
  parms_d_symp$parm2 <- as.numeric(params.set[i,"d_symp"])
  
  for (subseq_interventions in c(background_intervention, "s","q")){      
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

data[,"Abs_Benefit"] <- data[,"R_s"] - data[,"R_q"]
data[,"Rel_Benefit"] <- data[,"Abs_Benefit"] / data[,"R_s"]
data[,"NNQ"] <- 1 / data[,"Abs_Benefit"]
data[data$NNQ < 1,"NNQ"] <- 1
data[data$NNQ > 9999,"NNQ"] <- 9999
data[data$NNQ == Inf,"NNQ"] <- 9999
data[,"Abs_Benefit_per_Qday"] <- data[,"Abs_Benefit"] / data[,"obs_to_iso_q"]
data$gamma <- params.set[,"gamma"]
data$prob_CT <- params.set[,"prob_CT"]
data$CT_delay <- params.set[,"CT_delay"]
data$epsilon <- params.set[,"epsilon"]
data$R_0_input <- params.set[,"R_0_input"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data$T_inc <- params.set[,"T_inc"]
data$T_lat <- params.set[,"T_lat"]
data$d_inf <- params.set[,"d_inf"]
data$d_symp <- params.set[,"d_symp"]

bonferroni.alpha <- 0.05/length(dimensions)
range_levels = c(0.25, 0.50, 0.75, 1)
prcc_data <- data.frame(matrix(rep(NA, 4*length(names)*length(dimensions)*length(range_levels)), ncol=4)) 
names(prcc_data) <- c("output","parameter","range","coef")
prcc_data$output <- rep(names, each = length(dimensions)*length(range_levels))
prcc_data$parameter <- rep(dimensions, times = length(names)*length(range_levels))
prcc_data$range <- rep(rep(range_levels, each=10), times=10) 
for (output in names){
  for (range in range_levels){
    for (parm in dimensions){
      data.short <- data
      data.short <- data.short[data.short[,as.character(parm)] < max(params.set[,parm])*range, ]
      prcc <- pcc(data.short[,(length(names)+1):length(names(data.short))], data.short[,output], nboot = 0, rank=TRUE, conf=1-bonferroni.alpha)
      summary <- print(prcc)$original
      prcc_data[prcc_data$output == output & prcc_data$range == range,"coef"] <- summary
    }
  }
}

require(ggplot2)
ggplot(prcc_data, aes(x = range, y= coef, )) +
  facet_grid(output ~ parameter) +
  geom_point() +
  geom_hline(yintercept=0, color="red", size=0.25) +
  theme_bw() 

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_7_20150806.RData")

#### Experiment 8: Pertussis Case Study #### 
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_8_20150806_pertussis1000reps.RData')

set.seed(2344)

background_intervention <- "u"

# Also consider uniform triangle Try to match serial interval when I know it
parms_pi_t <- list("uniform", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# Eyeball fit. Very rough.
parms_serial_interval <- list("weibull", 2, 22)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# Dunno
parms_R_0 = list("uniform", 1.1, 10, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# 7-10 days normally. Range is 4-21 though
parms_T_inc = list("triangle", 4, 21, 8.5, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Starts with onset of prodrome. Consider also T_lat preceding T_inc by up to 2 weeks (insidious onset)
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# They say around 21 days. I made up a +/- 5 day buffer
parms_d_inf = list("uniform", 15, 26, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# usually 4-11 weeks. peak 7 weeks. up to 16 weeks
parms_d_symp = list("triangle", 4*7, 7*7, 16*7, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
prob_CT <- 0.5

# Dunno
parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 500
num_generations <- 5
times <- 10
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0.5
gamma.max <- 1.0
prob_CT.min <- 0.25
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 5
epsilon.min <- 0
epsilon.max <- 7
R_0.min <- 1.5
R_0.max <- 6

params.set <- cbind(
  gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
  prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
  CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
  epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
  R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  
  # Set R_0 to be constant within trials
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
  prcc <- pcc(data[,(length(names)+1):length(names(data))], data[,output], nboot = 10, rank=TRUE, conf=1-bonferroni.alpha)
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

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_8_20150806_pertussis1000reps.RData")

#### Experiment 9: Measles Case Study ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_9_20150806_measles1000reps.RData')

set.seed(2344)

background_intervention <- "u"

# Also consider uniform triangle Try to match serial interval when I know it
parms_pi_t <- list("triangle_center", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# Dunno
parms_serial_interval <- list("unknown", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# Dunno
parms_R_0 = list("uniform", 1.1, 10, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# 7-10 days normally. Range is 4-21 though
parms_T_inc = list("triangle", 4, 21, 8.5, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Starts with onset of prodrome. Consider also T_lat preceding T_inc by up to 2 weeks (insidious onset)
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# They say around 21 days. I made up a +/- 5 day buffer
parms_d_inf = list("uniform", 15, 26, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# usually 4-11 weeks. peak 7 weeks. up to 16 weeks
parms_d_symp = list("triangle", 4*7, 7*7, 16*7, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
prob_CT <- 0.5

# Dunno
parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 200
num_generations <- 5
times <- 20
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("gamma","prob_CT","CT_delay","epsilon","R_0")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

gamma.min <- 0.5
gamma.max <- 1.0
prob_CT.min <- 0.25
prob_CT.max <- 1
CT_delay.min <- 0
CT_delay.max <- 5
epsilon.min <- 0
epsilon.max <- 7
R_0.min <- 1.5
R_0.max <- 6


params.set <- cbind(
  gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
  prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
  CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
  epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
  R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  gamma <- as.numeric(params.set[i,"gamma"])
  prob_CT <- as.numeric(params.set[i,"prob_CT"])
  parms_CT_delay$parm2 <- as.numeric(params.set[i,"CT_delay"])
  parms_epsilon$parm2 <- as.numeric(params.set[i,"epsilon"])
  
  # Set R_0 to be constant within trials
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
  prcc <- pcc(data[,(length(names)+1):length(names(data))], data[,output], nboot = 10, rank=TRUE, conf=1-bonferroni.alpha)
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

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_9_20150806_measles1000reps.RData")

#### Experiment 10: Optimize serial interval for Pertussis Case study ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_10_20150808_4000reps.RData')

# Treat the prodrome as the onset of symptoms. Because, for a contact under SM, they will be suspecting pertussis

set.seed(45464)

background_intervention = "u"

# Also consider uniform triangle Try to match serial interval when I know it
parms_pi_t <- list("uniform", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# Eyeball fit. Very rough.
parms_serial_interval <- list("weibull", 2, 22)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# Dunno
parms_R_0 = list("uniform", 1.1, 10, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# 7-10 days normally. Range is 4-21 though
parms_T_inc = list("triangle", 4, 21, 8.5, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Starts with onset of prodrome. Consider also T_lat preceding T_inc by up to 2 weeks (insidious onset)
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# They say around 21 days. I made up a +/- 5 day buffer
parms_d_inf = list("uniform", 15, 26, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# usually 4-11 weeks. peak 7 weeks. up to 16 weeks
parms_d_symp = list("triangle", 4*7, 7*7, 16*7, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
prob_CT <- 0.5

# Dunno
parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 50
num_generations <- 5
times <- 4000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_lower","d_symp_upper","d_inf_offset","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -4
T_lat_offset.max <- 4
d_inf_lower.min <- 5
d_inf_lower.max <- 20
d_inf_upper.min <- 21
d_inf_upper.max <- 30
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_inf_lower = lhs[,2]*(d_inf_lower.max - d_inf_lower.min) + d_inf_lower.min,
  d_inf_upper = lhs[,3]*(d_inf_upper.max - d_inf_upper.min) + d_inf_upper.min,
  pi_t_triangle_center = lhs[,4]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_inf$parm1 <- as.numeric(params.set[i,"d_inf_lower"])
  parms_d_inf$parm2 <- as.numeric(params.set[i,"d_inf_upper"])
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
data$d_inf_lower <- params.set[,"d_inf_lower"]
data$d_inf_upper <- params.set[,"d_inf_upper"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval Inspection
layout(cbind(c(1),c(2),c(3),c(4)))
plot(data$T_lat_offset, data$ks)
plot(data$d_inf_lower, data$ks)
plot(data$d_inf_upper, data$ks)
plot(data$pi_t_triangle_center, data$ks)
pcor(data[,c("ks","T_lat_offset","d_inf_lower","d_inf_upper","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_offset","d_inf_lower","d_inf_upper","pi_t_triangle_center")])$p.value[1,]

data <- data[data$T_lat_offset > -2 & data$T_lat_offset < 1,]
data <- data[data$T_lat_offset > -1.5 & data$T_lat_offset < -0.5,]
data <- data[data$d_inf_lower > 6,]
data <- data[data$d_inf_upper < 28,]

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_10_20150808_4000reps.RData")

#### Experiment 11: Optimize serial interval for Measles Case study ####
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_10_20150808_4000reps.RData')

# Assume that onset of prodrome is the onset of symptoms

set.seed(11)

background_intervention = "u"

# Also consider uniform triangle Try to match serial interval when I know it
parms_pi_t <- list("uniform", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# No knowledge
parms_serial_interval <- list("unknown", 2, 22)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# Dunno
parms_R_0 = list("uniform", 1.1, 10, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Average 10-12 days until prodrome onset. Rash is 2-4 days later. Range is 7-21 days for rash
parms_T_inc = list("uniform", 4, 21, 8.5, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Starts with onset of prodrome.
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# 8 days
parms_d_inf = list("uniform", 6, 10, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# prodrome lasts 2-4 days and rash lasts 5-6 days
parms_d_symp = list("uniform", 7, 10, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
prob_CT <- 0.5

# Dunno
parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Dunno
parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 50
num_generations <- 5
times <- 4000
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_lower","d_symp_upper","d_inf_offset","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -4
T_lat_offset.max <- 4
d_inf_lower.min <- 3
d_inf_lower.max <- 7
d_inf_upper.min <- 8
d_inf_upper.max <- 14
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 1

params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  d_inf_lower = lhs[,2]*(d_inf_lower.max - d_inf_lower.min) + d_inf_lower.min,
  d_inf_upper = lhs[,3]*(d_inf_upper.max - d_inf_upper.min) + d_inf_upper.min,
  pi_t_triangle_center = lhs[,4]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
  parms_d_inf$parm1 <- as.numeric(params.set[i,"d_inf_lower"])
  parms_d_inf$parm2 <- as.numeric(params.set[i,"d_inf_upper"])
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
data$d_inf_lower <- params.set[,"d_inf_lower"]
data$d_inf_upper <- params.set[,"d_inf_upper"]
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval Inspection
layout(cbind(c(1),c(2),c(3),c(4)))
plot(data$T_lat_offset, data$ks)
plot(data$d_inf_lower, data$ks)
plot(data$d_inf_upper, data$ks)
plot(data$pi_t_triangle_center, data$ks)
pcor(data[,c("ks","T_lat_offset","d_inf_lower","d_inf_upper","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_offset","d_inf_lower","d_inf_upper","pi_t_triangle_center")])$p.value[1,]

# data <- data[data$T_lat_offset > -4 & data$T_lat_offset < 2,]
# data <- data[data$d_symp_upper > 12,]
# data <- data[data$T_lat_offset > -1.5,]
# data <- data[data$pi_t_triangle_center < 0.4,]
# data <- data[data$d_symp_lower > 4,]

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_10_20150808_4000reps.RData")

#### Experiment 12: Optimize serial interval for Ebola Case study #### 
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_12_20150809_4000reps.RData')

set.seed(12)

background_intervention = "u"

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# Ebola serial interval approximation from WHO
parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# WHO NEJM
# Liberia R0 1.83 (1.72, 1.94)
parms_R_0 = list("uniform", 1.72, 1.94, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
# All countries: Incubation period
parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")
# Set T_lat = T_inc
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# # Ebola incubation period from Eichner
# parms_T_inc <- list("lognormal", 12.7, 4.31)
# names(parms_T_inc) <- c("dist","parm1","parm2")

# Dixon 2015 EID
prob_CT <- 0.3 #27.4 to 31.1% of case patients were in the contact registry before identification. around 30% of new cases reported having any contacts they could have infected 

# Chowell 2014 and Lekone 2006 use ~6.5 days
parms_d_symp = list("uniform", 3, 8, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Set d_symp = d_inf
parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 50
num_generations <- 5
times <- 500
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset", "d_symp_lower","d_symp_upper","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -2
T_lat_offset.max <- 2
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
data <- decile_fcn(data, params.set)
decile_plot_fcn(data, params.set)

# Narrow range in 10% increments
data <- data[data$T_lat_offset < sort(data$T_lat_offset)[round(nrow(data)*.9)],]
data <- data[data$T_lat_offset > sort(data$T_lat_offset)[round(nrow(data)*.1)],]
data <- data[data$d_symp_upper > sort(data$d_symp_upper)[round(nrow(data)*.1)],]
data <- data[data$pi_t_triangle_center < sort(data$pi_t_triangle_center)[round(nrow(data)*.9)],]
data <- data[data$d_symp_upper < sort(data$d_symp_upper)[round(nrow(data)*.9)],]
data <- data[data$d_symp_lower < sort(data$d_symp_lower)[round(nrow(data)*.9)],]

# Try to narrow the range of input parameters to improve serial interval fit
apply(data[,names(data.frame(params.set))], 2, min)
apply(data[,names(data.frame(params.set))], 2, max)

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_12_20150807_4000reps.RData")

#### Experiment 13: Optimize serial interval for Ebola Case study (Step 2) #### 
# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_13_20150809_4000reps.RData')

set.seed(13)

background_intervention = "u"

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

# Ebola serial interval approximation from WHO
parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

# WHO NEJM
# Liberia R0 1.83 (1.72, 1.94)
parms_R_0 = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# All countries: Incubation period
parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Set T_lat = T_inc
parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# # Ebola incubation period from Eichner
# parms_T_inc <- list("lognormal", 12.7, 4.31)
# names(parms_T_inc) <- c("dist","parm1","parm2")

# Dixon 2015 EID
prob_CT <- 0.3 #27.4 to 31.1% of case patients were in the contact registry before identification. around 30% of new cases reported having any contacts they could have infected 

# Chowell 2014 and Lekone 2006 use ~6.5 days
parms_d_symp = list("uniform", 1, 14, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

# Set d_symp = d_inf
parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

n_pop = 100
num_generations <- 5
times <- 50
names <- c("R_0","ks")
data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data) <- names
dimensions <- c("T_lat_offset","pi_t_triangle_center")

require(lhs)
lhs <- maximinLHS(times, length(dimensions))

T_lat_offset.min <- -1.5
T_lat_offset.max <- 0.5
pi_t_triangle_center.min <- 0
pi_t_triangle_center.max <- 0.3

params.set <- cbind(
  T_lat_offset = lhs[,1]*(T_lat_offset.max - T_lat_offset.min) + T_lat_offset.min,
  pi_t_triangle_center = lhs[,2]*(pi_t_triangle_center.max - pi_t_triangle_center.min) + pi_t_triangle_center.min)

for (i in 1:times){
  cat('\nIteration',i, '\n')
  
  parms_T_lat$anchor_value <- as.numeric(params.set[i,"T_lat_offset"])
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
data$pi_t_triangle_center <- params.set[,"pi_t_triangle_center"]
data_store <- data

# Serial Interval partial rank correlation
pcor(data[,c("ks","T_lat_offset","pi_t_triangle_center")])$estimate[1,]
pcor(data[,c("ks","T_lat_offset","pi_t_triangle_center")])$p.value[1,]

# Serial Interval Inspection Scatter
layout(cbind(c(1),c(2),c(3)))
plot(data$T_lat_offset, data$ks)
plot(data$pi_t_triangle_center, data$ks)

# Serial Interval Inspection Deciles
data <- decile_fcn(data, params.set)
decile_plot_fcn(data, params.set)

# Narrow range in 10% increments
data <- data[data$T_lat_offset < sort(data$T_lat_offset)[round(nrow(data)*.9)],]
data <- data[data$T_lat_offset > sort(data$T_lat_offset)[round(nrow(data)*.1)],]
data <- data[data$d_symp_upper > sort(data$d_symp_upper)[round(nrow(data)*.1)],]
data <- data[data$pi_t_triangle_center < sort(data$pi_t_triangle_center)[round(nrow(data)*.9)],]
data <- data[data$d_symp_upper < sort(data$d_symp_upper)[round(nrow(data)*.9)],]
data <- data[data$d_symp_lower < sort(data$d_symp_lower)[round(nrow(data)*.9)],]

# Try to narrow the range of input parameters to improve serial interval fit
apply(data[,names(data.frame(params.set))], 2, min)
apply(data[,names(data.frame(params.set))], 2, max)

# save.image("~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_13_20150807_4000reps.RData")
