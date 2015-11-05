#### Header ####
# Make a spaghetti plot with many trials starting with one individual and tracing how many are in her infection tree as a function of generations. Then apply an intervention in generation 4 or so. Color code by u, hsb, s, q

#### Load PlotTrajectories Workspace ####
desired_root <- "20151028_InfluenzaA"
# desired_root <- "20151028_InfluenzaA"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTrajectories.RData", sep=""))

#### Load SMC Workspaces ####
desired_root <- "20151028_InfluenzaA"
# desired_root <- "20151028_InfluenzaA"
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))

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

parms_R_0 = list("uniform", 1.83, 1.83, 999, "independent", "independent") 
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
background_intervention = "u"

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

parms_R_0 = list("uniform", 1.46, 1.46, 999, "independent", "independent") 
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

#### Intervention ####
prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

#### Model settings ###
n_pop = 100
num_generations_u <- 2
num_generations_intervention <- 4
times = 25
names <- c("generation", "count", "trial", "intervention")
data_traj <- data.frame(matrix(rep(NA, length(names)*times*(num_generations_u + num_generations_intervention - 1)*4), nrow=times*(num_generations_u + num_generations_intervention - 1)*4))
names(data_traj) <- names
data_traj$generation <- rep(seq(1, (num_generations_u + num_generations_intervention - 1)), times = times*4)
data_traj$trial <- rep(seq(1:(times*4)), each = (num_generations_u + num_generations_intervention - 1))
data_traj$intervention <- rep(c("u", "hsb", "s", "q"), each = times * (num_generations_u + num_generations_intervention - 1))
data_traj$intervention <- factor(data_traj$intervention, levels = c("u", "hsb", "s", "q"))

#### Draw the median from SMC to represent each iteration
# This is better than each iteration having a different set of attributes because the stochasticity in the model can be shown instead of just the uncertainty in the model inputs.
parms_T_lat$anchor_value <- median(data[,"T_lat_offset"])
parms_d_inf$parm2 <- median(data[,"d_inf"])
parms_pi_t$triangle_center <- median(data[,"pi_t_triangle_center"])

#### Run model ####
for (i in 1:times){
  cat("\ni = ", i)
  
  for (subseq_interventions in c("u", "hsb", "s","q")){      
    In_Out.pre <- repeat_call_fcn(n_pop = n_pop, 
                                  parms_T_inc = parms_T_inc, 
                                  parms_T_lat = parms_T_lat, 
                                  parms_d_inf = parms_d_inf, 
                                  parms_d_symp = parms_d_symp, 
                                  parms_R_0 = parms_R_0, 
                                  parms_epsilon = parms_epsilon, 
                                  parms_pi_t = parms_pi_t,
                                  num_generations = num_generations_u,
                                  background_intervention = background_intervention,
                                  subseq_interventions = "u",
                                  gamma = gamma,
                                  prob_CT = prob_CT,
                                  parms_CT_delay = parms_CT_delay,
                                  parms_serial_interval = parms_serial_interval,
                                  cap_pop = FALSE,
                                  min_infections = 1,
                                  printing = FALSE,
                                  dispersion = 1)
    pre <- In_Out.pre$output[,"n"]
    
    In_Out.post <- repeat_call_fcn(n_pop = In_Out.pre$output[num_generations_u,"n"], 
                                   parms_T_inc = parms_T_inc, 
                                   parms_T_lat = parms_T_lat, 
                                   parms_d_inf = parms_d_inf, 
                                   parms_d_symp = parms_d_symp, 
                                   parms_R_0 = parms_R_0, 
                                   parms_epsilon = parms_epsilon, 
                                   parms_pi_t = parms_pi_t,
                                   num_generations = num_generations_intervention,
                                   background_intervention = background_intervention,
                                   subseq_interventions = subseq_interventions,
                                   gamma = gamma,
                                   prob_CT = prob_CT,
                                   parms_CT_delay = parms_CT_delay,
                                   parms_serial_interval = parms_serial_interval,
                                   cap_pop = FALSE,
                                   min_infections = 1,
                                   printing = FALSE,
                                   dispersion = 1)
    post <- In_Out.post$output[,"n"]
    
    count <- c(pre, post[-1])
    if (length(count) < (num_generations_u + num_generations_intervention - 1)){
      cat(" [Error: count is too short]")
      diff <- (num_generations_u + num_generations_intervention - 1) - length(count)
      count <- c(count, rep(0, diff))
    }
    data_traj[is.na(data_traj$count) == 1 & data_traj$intervention == subseq_interventions, "count"][1:length(count)] <- count
  }
}

#### Plot Results ####
plot_traj <- ggplot(data_traj, aes(x=generation, y=count, group=trial, color=intervention)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(x=(num_generations_u+1), col = "grey") +
  geom_line(alpha = 0.3) +
  stat_smooth(aes(group = intervention), size = 2, method = "loess", se = FALSE) +
  scale_x_continuous(breaks = seq(1:length(count))) +
  xlab("Generation") + ylab("Incident Cases") +
  scale_color_manual(name="Intervention",
                     values = c("coral3", "mediumturquoise", "darkgoldenrod", "cornflowerblue"),
                       breaks=c("u", "hsb", "s", "q"),
                       labels=c("None", "Health-\nSeeking\nBehavior", "Symptom\nMonitoring", "Quarantine")) +
  theme(legend.direction = "vertical", 
        legend.position = "right",
        legend.key=element_rect(size=5, color="white"))+
  ggtitle(paste(disease)) +
  scale_y_continuous(breaks = seq(0, 1400, by=100))
plot_traj

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTraj.pdf", sep=""))
plot(plot_traj)
dev.off()

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotTrajectories.RData", sep=""))



