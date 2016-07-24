# Figure_Fraser2004

#### Load libraries ####
library(ggplot2)

#### Define root ####
desired_root <- "20151024_Ebola" # Paste the desired root here "YYYYMMDD_DISEASE"
root <- desired_root

#### Load _Figure_Doomsday.RData Workspace ####
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "/", root, "_PlotDoomsday_Ideal.RData", sep=""))

#### Load Fraser.R Workspace #### 
load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "/", root, "_Fraser.RData", sep=""))

#### Calculate Fraser 2004 Theta value ####
# Theta is the proprtion of infections that occur prior to symptom onset
data.doomsday$theta <- NA

for (i in 1:nrow(data.doomsday)){
  a <- 0
  b <- data.doomsday[i, "d_inf"]
  c <- data.doomsday[i, "pi_t_triangle_center"] * data.doomsday[i, "d_inf"]
  x <- -(data.doomsday[i, "T_lat_offset"])   # Time of symptom onset
  
  if (data.doomsday[i, "T_lat_offset"] >= 0){   # If infectious after symptom onset
    data.doomsday[i, "theta"] = 0
  } else if (b <= x){   # If end of infectiousness is before symptom onset
    data.doomsday[i, "theta"] = 1
  } else {
    if (c <= x){   # If the peak occurs before onset of symptoms
      # Calculate the right hand side triangle
      height <- 2*(b-x) / ((b-a)*(b-c))
      base <- b-x
      rhs <- 0.5*base*height
      data.doomsday[i, "theta"] <- 1-rhs
    }
    if (c > x){
      # Calculate the left hand side triangle
      height <- 2*(x-a) / ((b-a)*(c-a))
      base <- x
      lhs <- 0.5*base*height
      data.doomsday[i, "theta"] <- lhs
    }
  }
}

#### Fraser 2004 Plot ####
plot_doomsday <- ggplot(data.doomsday) +
  geom_point(data = data.doomsday[data.doomsday$R_q < 1,], aes(theta, R_0), color = "cornflowerblue", size = 3) +
  geom_point(data = data.doomsday[data.doomsday$R_s < 1,], aes(theta, R_0), color = "darkgoldenrod", size = 3) +  
  geom_point(data = data.doomsday[data.doomsday$R_hsb < 1,], aes(theta, R_0), color = "mediumturquoise", size = 3) +
  theme_bw() + xlim(c(0,1)) +
  ylim(c(0, max(data.doomsday$R_0)))
plot_doomsday

#### Zero adapted Fraser 2004 plot ####
# Need some way of showing what's happening in the theta = 0. Need to show the proportion of trials where there was control.
# Make a heatmap with the minimum sufficient control method for each grid space of theta and R_0. 
# So, the theta=1 colum may be green first, then brown, then blue, then nothing when no intervention is sufficient.
# Set the threshold to 90% of trials, there is control.

# Find max R_0 that HSB maintains 90% control
layout(c(1))
plot(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_hsb > 1, "R_0"], col="red")
points(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_hsb < 1, "R_0"], col="green", add = "TRUE")

layout(cbind(c(1,2)))
no_control <- hist(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_hsb >= 1, "R_0"], breaks = seq(1, 10, 0.2))
control <- hist(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_hsb < 1, "R_0"], breaks = seq(1, 10, 0.2))

control$density / (control$density + no_control$density)
HSB_max <- control$breaks[max(which(control$density / (control$density + no_control$density) > 0.90))]

# Find max R_0 that SM maintains 90% control
plot(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_s > 1, "R_0"], col="red")
points(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_s < 1, "R_0"], col="green", add = "TRUE")

layout(cbind(c(1,2)))
no_control <- hist(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_s >= 1, "R_0"], breaks = seq(1, 10, 0.2))
control <- hist(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_s < 1, "R_0"], breaks = seq(1, 10, 0.2))

if (sum(no_control$counts)==0){
  S_max <- max(control$breaks)
} else {
  control$density / (control$density + no_control$density)
  S_max <- control$breaks[max(which(control$density / (control$density + no_control$density) > 0.90))]
}

# Find max R_0 that Q maintains 90% control
plot(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_q > 1, "R_0"], col="red")
points(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_q < 1, "R_0"], col="green", add = "TRUE")

layout(cbind(c(1,2)))
no_control <- hist(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_q >= 1, "R_0"], breaks = seq(1, 10, 0.2))
control <- hist(data.doomsday[data.doomsday$theta == 0 & data.doomsday$R_q < 1, "R_0"], breaks = seq(1, 10, 0.2))

if (sum(no_control$counts)==0){
  Q_max <- max(control$breaks)
} else {
  control$density / (control$density + no_control$density)
  Q_max <- control$breaks[max(which(control$density / (control$density + no_control$density) > 0.90))]
}

# Plot
plot_fraser <- ggplot() +
  geom_point(data = data.doomsday[data.doomsday$R_q < 1 & data.doomsday$R_0 < Q_max & data.doomsday$theta == 0,],  aes(theta-0.02, R_0), color = "cornflowerblue", size = 1.4) +
  geom_point(data = data.doomsday[data.doomsday$R_q < 1 & data.doomsday$theta > 0,], aes(theta, R_0), color = "cornflowerblue", size = 1.4) +
  geom_point(data = data.doomsday[data.doomsday$R_s < 1 & data.doomsday$R_0 < S_max & data.doomsday$theta == 0,],  aes(theta-0.01, R_0), color = "darkgoldenrod", size = 1.4) +
  geom_point(data = data.doomsday[data.doomsday$R_s < 1 & data.doomsday$theta > 0,], aes(theta, R_0), color = "darkgoldenrod", size = 1.4) +  
  geom_point(data = data.doomsday[data.doomsday$R_hsb < 1 & data.doomsday$R_0 < HSB_max & data.doomsday$theta == 0,],  aes(theta, R_0), color = "mediumturquoise", size = 1.4) +
  geom_point(data = data.doomsday[data.doomsday$R_hsb < 1 & data.doomsday$theta > 0,], aes(theta, R_0), color = "mediumturquoise", size = 1.4) +
  theme_bw() + 
  xlim(c(-0.02,1)) +
  xlab(expression(paste("Proportion of infections that occur prior to\nsymptoms or by asymptomatic infection (",theta,")", sep = ""))) +
  theme(axis.title.x=element_text(margin=margin(20,0,0,0))) +
  ylim(c(0, 10)) +
  ylab(expression(R[0]))
plot_fraser

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotFraser.pdf", sep=""), width = 3, height = 3)
plot(plot_fraser)
dev.off()

pdf(file=paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_PlotFraser_5by5.pdf", sep=""), width = 5, height = 5)
plot(plot_fraser)
dev.off()

#### Save _Fraser.RData Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "_Fraser.RData", sep=""))

