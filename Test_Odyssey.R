# Test for odyssey

cat("hello world")

# First, try to run Functions.R
source("Functions.R")

#### Define a sample of initial parameters ####
n_pop = 500

parms_pi_t <- list("triangle", 0.50)
names(parms_pi_t) <- c("distribution","triangle_center")

parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

parms_R_0 = list("uniform", 1.72, 1.94, 999, "independent", "independent")
names(parms_R_0) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_inc = list("gamma", 1.75, 0.182, 999, "independent", "independent")
names(parms_T_inc) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_T_lat = list("triangle", 999, 999, 999, 0, "T_inc")
names(parms_T_lat) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_inf = list("confit", 7.5, 6.8, 999, "independent", "independent")
names(parms_d_inf) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

prob_CT <- 0.3 #27.4 to 31.1% of case patients were in the contact registry before identification. around 30% of new cases reported having any contacts they could have infected 

parms_d_inf = list("triangle", 3, 8, 6.5, "independent", "independent")
names(parms_d_inf) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_d_symp = list("uniform", 999, 999, 999, 0, "d_inf")
names(parms_d_symp) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_epsilon = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_epsilon) <- c("dist","parm1","parm2",  "parm3","anchor_value", "anchor_target")

parms_CT_delay = list("uniform", 0, 1, 999, "independent", "independent")
names(parms_CT_delay) <- c("dist", "parm1", "parm2",  "parm3","anchor_value", "anchor_target")

parms_serial_interval <- list("gamma", 2.5, 0.2)
names(parms_serial_interval) <- c("dist","parm1","parm2")

background_intervention <- "u"

#### Test Create_Pop ####
Pop <- Create_Pop(n_pop=500, parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation = 1, background_intervention = "u", parms_CT_delay, gamma)
Pop[1,]
hist(Pop$T_inc)
hist(Pop$T_lat)
hist(Pop$R_0)
cat(length(which(Pop[,"R_0"]<0)), "individuals have an R_0 < 0")

