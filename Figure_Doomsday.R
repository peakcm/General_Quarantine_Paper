#### Header ####
# Figure for doomsday (T_lat_offset) on Y and (R_0) on X

#### Load libraries ####
library(ggplot2)

#### Load data ####
# Use the SMC parameter space
desired_root <- "20151024_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"
root <- desired_root
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_root, "/", desired_root, "_SMC.RData", sep=""))

#### Load Workspace ####
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "/", root, "_Doomsday.RData", sep=""))
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "/", root, "_PlotDoomsday_Ideal.RData", sep=""))

#### Define Parms ####
# Set a range of T_lat_offset we're interested in
T_lat_offset.min <- -12
T_lat_offset.max <- 12

# Set a range of R_0
R_0.min <- 1
R_0.max <- 10

dispersion = 1

#### Define Interventions ####
# Optimal Resource Interventions
background_intervention = "u"

prob_CT <- 1

gamma <- 1

parms_epsilon = list("uniform", 0.5, 0.5, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 0, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# High Resource Interventions
background_intervention = "u"

prob_CT <- 0.9

gamma <- 0.9

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Mid Resource Interventions
background_intervention <- "u"

prob_CT <- 0.75

gamma <- 0.75

parms_epsilon = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 2, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

# Low Resource Interventions
background_intervention <- "u"

prob_CT <- 0.5

gamma <- 0.5

parms_epsilon = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 4, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

#### Resample ####
# Resample from _SMC except without respecting joint distribution and new T_lat_offset

# Initialize
n_pop = 500
num_generations <- 5
times <- 500
names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday", "ks")
data.doomsday <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
names(data.doomsday) <- names

# Sample from posterior distributions for each parameter independently
params.set.doomsday <- cbind(
  T_lat_offset = runif(n = times, min = T_lat_offset.min, max = T_lat_offset.max),
  d_inf = sample(data$d_inf, size = times, replace = TRUE),
  pi_t_triangle_center = sample(data$pi_t_triangle_center, size = times, replace = TRUE),
  R_0 = runif(n = times, min = R_0.min, max = R_0.max))

params.set.doomsday <- data.frame(params.set.doomsday)

i=1
x=1
while (i <= times){
  
  parms_pi_t$triangle_center <- as.numeric(params.set.doomsday[i,"pi_t_triangle_center"])
  parms_T_lat$anchor_value <- as.numeric(params.set.doomsday[i,"T_lat_offset"])
  parms_d_inf$parm2 <- as.numeric(params.set.doomsday[i,"d_inf"])
  parms_R_0[c("parm1","parm2")] <- c(as.numeric(params.set.doomsday[i,"R_0"]), as.numeric(params.set.doomsday[i,"R_0"]))

  for (subseq_interventions in c("hsb", "s","q")){      
    
    if (subseq_interventions == background_intervention & parms_R_0$parm1 > 1){
      n_pop_input <- 200
    } else if (subseq_interventions == "hsb" & parms_R_0$parm1 * (1-gamma) > 1){ 
      n_pop_input <- 200
    } else if ((subseq_interventions == "s" | subseq_interventions == "q") & parms_R_0$parm1 * (1-gamma*prob_CT) > 1.1){
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
    data.doomsday[i,"R_0"] <- parms_R_0[c("parm1")]
    if (subseq_interventions == "hsb"){
      data.doomsday[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "s"){
      data.doomsday[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
    }
    if (subseq_interventions == "q"){
      data.doomsday[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
      data.doomsday[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
    }
  }
  
  if (is.na(data.doomsday[i,"R_0"])==1 |
      is.na(data.doomsday[i,"R_hsb"])==1 |
      is.na(data.doomsday[i,"R_s"])==1 |
      is.na(data.doomsday[i,"R_q"])==1){
    i=i  #re-run that set
    cat("x")
    x = x+1
    if (x > 3){
      data.doomsday[i,c("R_0", "R_hsb","R_s","R_q")] <- 0
      i = i+1
      x=1
    }
  } else {
    i = i+1
    x = 1
    cat(".")
    if (i%%10 == 0){cat("|")}
    if (i%%100 == 0){cat("\n")}
  }
  
}

# Check for missing data
if (sum(is.na(data.doomsday[,1:4]))>0){cat("Something's missing")}

data.doomsday[,"Abs_Benefit"] <- data.doomsday[,"R_s"] - data.doomsday[,"R_q"]
data.doomsday[,"Rel_Benefit"] <- data.doomsday[,"Abs_Benefit"] / data.doomsday[,"R_s"]
data.doomsday[,"NNQ"] <- 1 / data.doomsday[,"Abs_Benefit"]
data.doomsday[data.doomsday$NNQ < 1,"NNQ"] <- 1
data.doomsday[data.doomsday$NNQ > 9999,"NNQ"] <- 9999
data.doomsday[data.doomsday$NNQ == Inf,"NNQ"] <- 9999
data.doomsday[,"Abs_Benefit_per_Qday"] <- data.doomsday[,"Abs_Benefit"] / data.doomsday[,"obs_to_iso_q"]
data.doomsday$R_0 <- params.set.doomsday[,"R_0"]
data.doomsday$pi_t_triangle_center <- params.set.doomsday[,"pi_t_triangle_center"]
data.doomsday$T_lat_offset <- params.set.doomsday[,"T_lat_offset"]
data.doomsday$d_inf <- params.set.doomsday[,"d_inf"]

#### Plot Orientation 1 ####
plot_doomsday <- ggplot(data.doomsday) +
  geom_point(data = data.doomsday[data.doomsday$R_q < 1,], aes(R_0, -T_lat_offset), color = "cornflowerblue", size = 3) +
  geom_point(data = data.doomsday[data.doomsday$R_s < 1,], aes(R_0, -T_lat_offset), color = "darkgoldenrod", size = 3) +  
  geom_point(data = data.doomsday[data.doomsday$R_hsb < 1,], aes(R_0, -T_lat_offset), color = "mediumturquoise", size = 3) +
  theme_bw() + xlim(0, R_0.max) 
plot_doomsday

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotDoomsday_Ideal_1.pdf", sep=""))
plot(plot_doomsday)
dev.off()

#### Plot Orientation 2 ####
point_size <- 0.3
plot_doomsday <- ggplot(data.doomsday) +
  geom_point(data = data.doomsday[data.doomsday$R_q < 1,], aes(-T_lat_offset, R_0), color = "cornflowerblue", size = point_size) +
  geom_point(data = data.doomsday[data.doomsday$R_s < 1,], aes(-T_lat_offset, R_0), color = "darkgoldenrod", size = point_size) +  
  geom_point(data = data.doomsday[data.doomsday$R_hsb < 1,], aes(-T_lat_offset, R_0), color = "mediumturquoise", size = point_size) +
  theme_bw() + xlim(c(min(-data.doomsday$T_lat_offset), max(-data.doomsday$T_lat_offset))) +
  ylim(c(0, 10)) + 
  theme(text = element_text(size=8)) +
  xlab(expression(paste(T[OFFSET], " (days)"))) +
  ylab(expression(R[0]))
plot_doomsday

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotDoomsday_Ideal_2.pdf", sep=""), width = 2, height = 2)
plot(plot_doomsday)
dev.off()

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotDoomsday_Ideal_2_5by5.pdf", sep=""), width = 5, height = 5)
plot(plot_doomsday)
dev.off()

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotDoomsday_Ideal.RData", sep=""))
