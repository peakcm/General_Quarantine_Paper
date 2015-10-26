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
library(reshape)

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

#### Experiment 10: Optimize serial interval and incubation period for Pertussis Case study ####

# Treat the prodrome as the onset of symptoms. Because, for a contact under SM, they will be suspecting pertussis

# serial interval approximation
# Data from te Beest 2014
n = 239
mean = 25.2
q_05 <- 3
q_25 <- 10
q_50 <- 20
q_75 <- 31
q_95 <- 69.1

sim_dist_1 <- c(rep(q_05, 5), rep(q_25, 25), rep(q_50, 25), rep(q_75, 25), rep(q_95, 5))

hist(sim_dist_1, breaks = seq(0, 70, 5), freq = FALSE)

fitdistr(sim_dist_1, densfun = "gamma")

# serial interval approximation
parms_serial_interval <- list("gamma", 2.45585778, .11071164) # Approximation from te Beest reported quantiles
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(add=TRUE, dgamma(x, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2),
      from=0, to=70, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")
summary(rgamma(1000000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))
hist(rgamma(1000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2), breaks = seq(0, 150, 5))

# Data from te Stocks 1933
# Note that this data is most likely biased towards early periods since simultaneous infection events likely account for many of the short intervals.
day <- c(0,1,2,3,4,7, 13.5, 19, 24, 31, 38, 66, 135)
counts <- c(465, 26, 12, 18, 35, 286, 231, 42, 61, 34, 18, 32, 31)
freq <- c(465, 26, 12, 18, 35, 57.2, 28.9, 14, 8.7, 4.9, 2.6, 1.6, 0.3)

stocks_data <- cbind(day, counts, freq)

# Using counts
sim_dist_2 <- c()
for (i in 1:length(day)){sim_dist_2 <- c(sim_dist_2, rep(day[i], times=counts[i]))}

hist(sim_dist_2[sim_dist_2 > 2], breaks = seq(0, 136, 2), freq = FALSE)
hist(sim_dist_2[sim_dist_2 > 2 & sim_dist_2 <=48], breaks = seq(0, 48, 1), freq = FALSE)

fitdistr(sim_dist_2[sim_dist_2 > 2], densfun = "gamma")

# serial interval approximation
parms_serial_interval <- list("gamma", 1.2959, 0.06541) # Approximation from Stocks 1933 table
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(add=TRUE, dgamma(x, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2),
      from=0, to=48, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")
summary(rgamma(1000000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))
hist(rgamma(1000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))

# Using "freq"
sim_dist_3 <- c()
for (i in 1:length(day)){sim_dist_3 <- c(sim_dist_3, rep(day[i], times=round(freq[i])))}

hist(sim_dist_3[sim_dist_3 > 2], breaks = seq(0, 136, 2), freq = FALSE)
hist(sim_dist_3[sim_dist_3 > 2 & sim_dist_3 <= 48], breaks = seq(0, 48, 1), freq = FALSE)

fitdistr(sim_dist_3[sim_dist_3 > 2], densfun = "gamma")

# serial interval approximation
parms_serial_interval <- list("gamma", 1.89698045, 0.17471522) # Approximation from Stocks 1933 table
names(parms_serial_interval) <- c("dist","parm1","parm2")

curve(add=TRUE, dgamma(x, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2),
      from=0, to=84, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")
summary(rgamma(1000000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))
hist(rgamma(1000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))

# Using Vink 2014 AJE
parms_serial_interval <- list("gamma", 2, 0.06) # Approximation from Vink's data and distribution
names(parms_serial_interval) <- c("dist","parm1","parm2")

hist(rgamma(1000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2), breaks = seq(0, 150, 5), freq=FALSE)
curve(add=TRUE, dgamma(x, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2),
      from=0, to=150, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")
summary(rgamma(1000000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))
sort(rgamma(1000000, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))[1000000*0.95]


# Incubation period
# 7-10 days normally. Range is 4-21 though
parms_T_inc = list("normal", 7, (10-7)/1.96, 999, "independent", "independent") #Gordon
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

curve(add=FALSE, dnorm(x, mean = parms_T_inc$parm1, sd =parms_T_inc$parm2),
      from=0, to=21, col="green", lwd=2,
      main = "Testing Desired Distribution", xlab = "Serial Interval (Days)", ylab = "Desired Distribution")
summary(rnorm(1000000, mean=parms_T_inc$parm1, sd=parms_T_inc$parm2))
sort(rnorm(1000000, mean=parms_T_inc$parm1, sd=parms_T_inc$parm2))[1000000*0.025]
sort(rnorm(1000000, mean=parms_T_inc$parm1, sd=parms_T_inc$parm2))[1000000*0.975]

# Also consider uniform triangle Try to match serial interval when I know it
parms_pi_t <- list("uniform", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

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

#### Experiment 16 ####
# Make a spaghetti plot with many trials starting with one individual and tracing how many are in his infection tree as a function of generations. Then apply an intervention in generation 4 or so. Color code by u, hsb, s, q

# load('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_16_20151001.RData')
# save.image('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Experiments_16_20151001.RData')

# Ebola-ish
parms_serial_interval <- list("gamma", 2.5, 0.2) # approximation from WHO
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 1, 3, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_T_lat = list("triangle", 999, 999, 999, -2, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("triangle", 1, 15, 8, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 3, 8, 999, 0, "d_symp")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 1, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Settings
n_pop = 50
num_generations_u <- 3
num_generations_intervention <- 4
times = 25
names <- c("generation", "count", "trial", "intervention")
data <- data.frame(matrix(rep(NA, length(names)*times*(num_generations_u + num_generations_intervention - 1)*4), nrow=times*(num_generations_u + num_generations_intervention - 1)*4))
names(data) <- names
data$generation <- rep(seq(1, (num_generations_u + num_generations_intervention - 1)), times = times*4)
data$trial <- rep(seq(1:(times*4)), each = 6)
data$intervention <- rep(c("u", "hsb", "s", "q"), each = times * (num_generations_u + num_generations_intervention - 1))
data$intervention <- factor(data$intervention, levels = c("u", "hsb", "s", "q"))

for (i in 1:times){
  for (subseq_interventions in c("u", "hsb", "s","q")){      
    In_Out.pre <- repeat_call_fcn(n_pop, 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations = num_generations_u,
                              background_intervention,
                              subseq_interventions = "u",
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              cap_pop = FALSE,
                              min_infections = 1)
    pre <- In_Out.pre$output[,"n"]
    
    In_Out.post <- repeat_call_fcn(n_pop = In_Out.pre$output[num_generations_u,"n"], 
                              parms_T_inc, 
                              parms_T_lat, 
                              parms_d_inf, 
                              parms_d_symp, 
                              parms_R_0, 
                              parms_epsilon, 
                              parms_pi_t,
                              num_generations = num_generations_intervention,
                              background_intervention,
                              subseq_interventions,
                              gamma,
                              prob_CT,
                              parms_CT_delay,
                              parms_serial_interval,
                              cap_pop = FALSE,
                              min_infections = 1)
    post <- In_Out.post$output[,"n"]
    
    count <- c(pre, post[-1])
    
    data[is.na(data$count) == 1 & data$intervention == subseq_interventions, "count"][1:length(count)] <- count
  }
}

ggplot(data, aes(x=generation, y=count, group=trial, color=intervention)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(x=4, col = "grey") +
  geom_line(alpha = 0.3) +
  stat_smooth(aes(group = intervention), size = 2, method = "loess", se = FALSE) +
  scale_x_continuous(breaks = seq(1:length(count))) +
  xlab("Generation") + ylab("Incident Cases") +
  scale_color_discrete(name="Intervention",
                    breaks=c("u", "hsb", "s", "q"),
                    labels=c("None", "Health-\nSeeking\nBehavior", "Symptom\nMonitoring", "Quarantine")) +
  theme(legend.direction = "vertical", 
        legend.position = "right",
        legend.key=element_rect(size=5, color="white"),
        legend.key.size = unit(2, "lines"))
  
