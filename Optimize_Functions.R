#### Header ####
# Optimize Functions for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 19, 2015

#### Load Libraries ####
library(reshape2)
library(ggplot2)
library(shiny)
library(lineprof)

#### Optimize the number of generations ####
ptm <- proc.time()

set.seed(847)

num_generations <- 10
times <- 250

background_intervention <- "u"

parms_T_inc = list("uniform", 1, 40, 99, "independent", "independent")
names(parms_T_inc) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_T_lat = list("uniform", 1, 40, 999, "independent", "independent")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 20, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 20, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 0.8, 15, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

T_inc_vector <- runif(times, min=5, max=40)
T_lat_vector <- runif(times, min=5, max=40)
d_inf_vector <- runif(times, min=5, max=20)
d_symp_vector <- runif(times, min=5, max=20)
R_0_vector <- runif(times, min=1.2, max=15)
gamma_vector <- runif(times, min=0.5, max=1)
prob_CT_vector <- runif(times, min=0.5, max=1)
CT_delay_vector <- runif(times, min=0, max=7)
epsilon_vector <- runif(times, min=0, max=2)

layout(rbind(c(1,2,3,4),c(5,6,7,8)))
for (pi_t_distribution in c("linear_increase","uniform")){
  for (subseq_interventions in c("u","hsb","s","q")){
    
    names <- c(1:num_generations, "intervention")
    data_R <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
    names(data_R) <- names
    
    names <- c(1:num_generations, "intervention")
    data_R_vs_latter <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
    names(data_R_vs_latter) <- names
    
    for (i in 1:times){
      cat('\nIteration',i, '\n')
      
      parms_T_inc$parm2 <- T_inc_vector[i]
      parms_T_lat$parm2 <- T_lat_vector[i]
      parms_d_inf$parm2 <- d_inf_vector[i]
      parms_d_symp$parm2 <- d_symp_vector[i]
      parms_R_0$parm2 <- R_0_vector[i]
      gamma <- gamma_vector[i]
      prob_CT <- prob_CT_vector[i]
      parms_CT_delay$parm2 <- CT_delay_vector[i]
      parms_epsilon$parm2 <- epsilon_vector[i]
      
      In_Out <- repeat_call_fcn(n_pop=200, 
                                parms_T_inc, 
                                parms_T_lat, 
                                parms_d_inf, 
                                parms_d_symp, 
                                parms_R_0, 
                                parms_epsilon, 
                                num_generations,
                                background_intervention,
                                subseq_interventions,
                                pi_t_distribution,
                                gamma,
                                prob_CT,
                                parms_CT_delay)
      
      data_R[i,"intervention"] <- subseq_interventions
      data_R_vs_latter[i,"intervention"] <- subseq_interventions
      data_R[i,1:length(In_Out$output$R)] <- In_Out$output$R
      for (j in 1:length(In_Out$output$R)){
        data_R_vs_latter[i,j] <- (data_R[i,j] - mean(as.numeric(data_R[i,j:length(In_Out$output$R)]))) / mean(as.numeric(data_R[i,j:length(In_Out$output$R)]))
      }
    }
    
    data_R_vs_latter_melt <- melt(data_R_vs_latter[,1:num_generations])
    plot(data_R_vs_latter_melt$variable, data_R_vs_latter_melt$value, las=3,
         ylab = "Difference from remaining (%)",
         xlab = "Generation",
         main = c(pi_t_distribution, subseq_interventions))
  }
}

proc.time() - ptm
# 10 generations, 250 iterations, 200 initial population, 5x max population
# Time elapsed = 23693.673 which is around 6.6 hours
# Output saved as Optimize_num_generations_250reps_npop200
# Conclusion: at least 3 generations are sufficient for "u" and "hsb"
# Conclusion: at least 4 generations are sufficient for "s" and "q"
# Conclusion: pi_t as uniform or linear_increase makes no difference

#### Optimize the initial population size ####
#Choose a population size that is big enough such that super effective interventions still have enough generations of data to analyze before the number drops below 20
#Choose a population size taht is small enough to allow quick analysis
#If the intervention brings the R down to 0.25, and we start with 800 people, then we will have 50 at the third generation and can still analyze it

ptm <- proc.time()

set.seed(22223)

num_generations <- 5
times <- 500

background_intervention <- "u"

parms_T_inc = list("uniform", 1, 40, 99, "independent", "independent")
names(parms_T_inc) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_T_lat = list("uniform", 1, 40, 999, "independent", "independent")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("uniform", 1, 20, 999, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 1, 20, 999, "independent", "independent")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_R_0 = list("uniform", 0.8, 15, 999, "independent", "independent")
names(parms_R_0) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

T_inc_vector <- runif(times, min=1, max=40)
T_lat_vector <- runif(times, min=1, max=40)
d_inf_vector <- runif(times, min=1, max=20)
d_symp_vector <- runif(times, min=1, max=20)
R_0_vector <- runif(times, min=1.1, max=2)
gamma_vector <- runif(times, min=0.75, max=1)
prob_CT_vector <- runif(times, min=0.75, max=1)
CT_delay_vector <- runif(times, min=0, max=1)
epsilon_vector <- runif(times, min=0, max=1)

n_pop_vector <- c(500)
pi_t_distribution_vector <- c("linear_increase")
subseq_interventions_vector <- c("q")

names <- c("R","n_pop","intervention","pi_t", "iteration")
data <- data.frame(matrix(rep(NA, length(names)*times*length(n_pop_vector)*length(pi_t_distribution_vector)*length(subseq_interventions_vector)), ncol=length(names)))
names(data) <- names

row=1
for (pi_t_distribution in pi_t_distribution_vector){
  for (subseq_interventions in subseq_interventions_vector){  
    for (i in 1:times){
      
      parms_T_inc$parm2 <- T_inc_vector[i]
      parms_T_lat$parm2 <- T_lat_vector[i]
      parms_d_inf$parm2 <- d_inf_vector[i]
      parms_d_symp$parm2 <- d_symp_vector[i]
      parms_R_0$parm2 <- R_0_vector[i]
      gamma <- gamma_vector[i]
      prob_CT <- prob_CT_vector[i]
      parms_CT_delay$parm2 <- CT_delay_vector[i]
      parms_epsilon$parm2 <- epsilon_vector[i]
      
      for (n_pop in n_pop_vector){
        cat('\nIteration',i,' n_pop=',n_pop,' intervention=',subseq_interventions, ' pi_t=', pi_t_distribution,'\n')
        In_Out <- repeat_call_fcn(n_pop=n_pop, 
                                  parms_T_inc, 
                                  parms_T_lat, 
                                  parms_d_inf, 
                                  parms_d_symp, 
                                  parms_R_0, 
                                  parms_epsilon, 
                                  num_generations,
                                  background_intervention,
                                  subseq_interventions,
                                  pi_t_distribution,
                                  gamma,
                                  prob_CT,
                                  parms_CT_delay)
        
        data[row,"n_pop"] <- n_pop
        data[row,"intervention"] <- subseq_interventions
        data[row,"pi_t"] <- pi_t_distribution
        data[row,"R"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
        data[row,"iteration"] <- as.character(i)
        
        row = row+1
      }
    }
  }
}

data_melt <- melt(data, id=c("pi_t","intervention","n_pop", "iteration"))
ggplot(data_melt, aes(x=n_pop, y=value, group=iteration, color=iteration)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(intervention ~ pi_t) +
  geom_line() +
  ylab("Equilibrium Reproductive Number") +
  xlab("Initial Population Size") +
  ggtitle("Optimize Initial Population Size")

proc.time() - ptm

data_melt[is.na(data_melt$value)==1,]

# Conclusion: People were dropped in trials when started with 100, 200, 300, 400 people. 
# Conclusion: Noone was dropped after 500 iterations of 500 people.

#### Optimize code speed ####
# Profile line seeds
source('~/Dropbox/Ebola/General_Quarantine_Paper/R_Code/Functions.R', echo=TRUE)
l <- lineprof(repeat_call_fcn(n_pop=500, 
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
                              parms_serial_interval))
shine(l)
