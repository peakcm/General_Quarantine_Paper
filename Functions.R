#### Header #### 
# Functions for Generalized Study of Quarantine vs Symptom Monitoring
# Corey Peak
# Version 1.0
# August 19, 2015
# Development Branch

#### Notes to work on #### 
# Measure elasticity. (% change in R for % change in attribute) (% change in R for % change in attribute variance)
#
# pi_t distributions:
# exp_1.1, exp_1.2, exp_1.3, exp_e do not seem to produce the right number of people in the second generation
# ebola data-driven distribution
#
# Add a duration of quaratine and symptom monitoring. don't isolate people if symptom onset is after T_obs + d_CT
# start the timer from when they are placed under S or Q. Alternatively, start timer at day of infection (as if they could guess)
# compare abs_benefit per Q day under conditions where we modify prob_CT, d_CT, and epsilon (how frequently you check ppl)

#### Draw_Dist_fcn #### 
# Draw from Distributions of Disease Attributes
Draw_Dist_fcn <- function(Vector, distribution, parm1, parm2, parm3){
  # Unit is days
  if (distribution == "gamma"){
    Draw_Dist <- rgamma(Vector, shape=parm1, rate=parm2)
    Draw_Dist[Draw_Dist < 0] <- 0
  }
  
  if (distribution == "uniform"){
    Draw_Dist <- runif(Vector, min=parm1, max=parm2)
    Draw_Dist[Draw_Dist < 0] <- 0
  }
  
  if (distribution == "confit"){
    # Use data from a confidence interval to recreate a normal distribution
    # parm1 is the point estimate
    # parm2 is the UPPER bound of the confidence interval
    Draw_Dist <- rnorm(Vector, mean=parm1, sd=(abs(parm2 - parm1))/1.96)
    Draw_Dist[Draw_Dist < 0] <- 0
    if (parm2 > parm1){Draw_Dist[Draw_Dist > parm2] <- parm2} #truncate at the upper bound of the confidence interval
  }
  
  if (distribution == "triangle"){
    # parm1 is a, the minimum
    # parm2 is b, the maximum
    # parm3 is c, the most likely
    if (parm3==999){cat('Error Draw_Dist_fcn: Must include parm3 if calling a triangle distribution')}
    intermediate <- runif(Vector,0,1) # A vector of values that correspond to u in Marc's equation from EPI 260
    Draw_Dist <- rep(NA, length(Vector))
    Draw_Dist[intermediate < ((parm3 - parm1)/(parm2-parm1))] <- parm1 + sqrt(intermediate[intermediate < ((parm3 - parm1)/(parm2-parm1))]*(parm2-parm1)*(parm3-parm1))
    Draw_Dist[intermediate >= ((parm3 - parm1)/(parm2-parm1))] <- parm2 - sqrt((1-intermediate[intermediate >= ((parm3 - parm1)/(parm2-parm1))])*(parm2-parm3)*(parm2-parm1))
  }
  
  if (distribution == "weibull"){
    Draw_Dist <- rweibull(Vector, shape=parm1, scale=parm2)
    Draw_Dist[Draw_Dist < 0] <- 0
  }
  
  if (distribution == "lognormal"){
    Draw_Dist <- rlnorm(Vector, meanlog = parm1, sdlog = parm2)
    Draw_Dist[Draw_Dist < 0] <- 0
  }
  return(Draw_Dist)
}

#### Anchor_fcn ####
# Define some parameter by anchoring to another parameter by X days before or after
Anchor_fcn <- function(Pop, anchor_value, anchor_target){
  #input anchor value in days
  output <- round(Pop[,anchor_target] + anchor_value*24)
}

#### Create_Pop ####
# Create population characteristics
Create_Pop <- function(n_pop,
                       parms_T_inc,
                       parms_T_lat,
                       parms_d_inf,
                       parms_d_symp,
                       parms_R_0,
                       parms_epsilon,
                       generation,
                       background_intervention,
                       parms_CT_delay,
                       gamma){
  names <- c("ID", "T_inc", "T_lat", "d_inf", "d_symp", "epsilon", "R_0", "R_0_hsb_adjusted", "generation", "infector", "t_obs_infector", "background_intervention", "CT_delay")
  Pop <- data.frame(matrix(rep(NA, n_pop*length(names)), ncol=length(names)))
  names(Pop) <- names
  Pop$ID      <- 1:nrow(Pop)
  
  # For time attributes, input in terms of days, then convert to hours, then round to the nearest hour.
  # If the parameters are drawn independently of other parameters, do that first
  if (parms_T_inc$anchor_value == "independent"){ Pop$T_inc   <- round(24*Draw_Dist_fcn(rep(NA, n_pop), parms_T_inc$dist, parms_T_inc$parm1, parms_T_inc$parm2, parms_T_inc$parm3)) }
  if (parms_T_lat$anchor_value == "independent"){ Pop$T_lat   <- round(24*Draw_Dist_fcn(rep(NA, n_pop), parms_T_lat$dist, parms_T_lat$parm1, parms_T_lat$parm2, parms_T_lat$parm3)) }
  if (parms_d_inf$anchor_value == "independent"){ Pop$d_inf   <- round(24*Draw_Dist_fcn(rep(NA, n_pop), parms_d_inf$dist, parms_d_inf$parm1, parms_d_inf$parm2, parms_d_inf$parm3)) }
  if (parms_d_symp$anchor_value == "independent"){ Pop$d_symp  <- round(24*Draw_Dist_fcn(rep(NA, n_pop), parms_d_symp$dist, parms_d_symp$parm1, parms_d_symp$parm2, parms_d_symp$parm2)) }
  if (parms_epsilon$anchor_value == "independent"){ Pop$epsilon <- round(24*Draw_Dist_fcn(rep(NA, n_pop), parms_epsilon$dist, parms_epsilon$parm1, parms_epsilon$parm2, parms_epsilon$parm3)) }
  if (parms_R_0$anchor_value == "independent"){ Pop$R_0     <- Draw_Dist_fcn(rep(NA, n_pop), parms_R_0$dist, parms_R_0$parm1, parms_R_0$parm2, parms_R_0$parm3) }
  if (parms_CT_delay$anchor_value == "independent"){ Pop$CT_delay <- round(24*Draw_Dist_fcn(rep(NA, n_pop), parms_CT_delay$dist, parms_CT_delay$parm1, parms_CT_delay$parm2, parms_CT_delay$parm3)) }
  
  # If the parameters are anchored to other parameters, add or subtract anchor_value from the anchor_target
  if (parms_T_inc$anchor_value != "independent"){Pop$T_inc <- Anchor_fcn(Pop, parms_T_inc$anchor_value, parms_T_inc$anchor_target)}
  if (parms_T_lat$anchor_value != "independent"){Pop$T_lat <- Anchor_fcn(Pop, parms_T_lat$anchor_value, parms_T_lat$anchor_target)}
  if (parms_d_inf$anchor_value != "independent"){Pop$d_inf <- Anchor_fcn(Pop, parms_d_inf$anchor_value, parms_d_inf$anchor_target)}
  if (parms_d_symp$anchor_value != "independent"){Pop$d_symp <- Anchor_fcn(Pop, parms_d_symp$anchor_value, parms_d_symp$anchor_target)}
  if (parms_epsilon$anchor_value != "independent"){Pop$epsilon <- Anchor_fcn(Pop, parms_epsilon$anchor_value, parms_epsilon$anchor_target)}
  if (parms_R_0$anchor_value != "independent"){Pop$R_0 <- Anchor_fcn(Pop, parms_R_0$anchor_value, parms_R_0$anchor_target)}
  if (parms_CT_delay$anchor_value != "independent"){Pop$CT_delay <- Anchor_fcn(Pop, parms_CT_delay$anchor_value, parms_CT_delay$anchor_target)}
  
  Pop$generation <- generation
  Pop$background_intervention <- background_intervention
  
  if (generation == 1){
    Pop$t_infection <- 0
  } else {Pop$t_infection <- NA}
  return(Pop)
}

#### t_obs and t_iso ####
# Time of Observation and Isolation
# There are four routes to isolation
t_obs_u_fcn <- function(T_inc, d_symp, T_lat, d_inf){
  # "u" denotes unisolated infection
  # Their offspring are only listed for contact tracing ("observed")
  # at the end of the infector's duration of disease (symptoms OR infectiousness)
  t_obs_u <- max( (T_inc + d_symp), (T_lat + d_inf) )
  return(t_obs_u)
}

t_obs_hsb_fcn <- function(T_inc, d_symp){
  # "hsb" denotes health seeking behavior
  # We assume an individual will seek health care at a random time during symptomatic disease
  if (is.na(d_symp)==1) {cat("Error 1: t_obs_hsb_fcn ")}
  if (is.na(T_inc)==1) {cat("Error 2: t_obs_hsb_fcn ")}
  if (T_inc >= (T_inc + d_symp)) {cat("Error 3: t_obs_hsb_fcn ")}
  t_obs_hsb <- runif(1, min = T_inc, max = (T_inc + d_symp) )
  return(t_obs_hsb)
}

t_obs_fcn <- function(t_infection, t_obs_alpha, CT_delay){
  # Both symptom monitoring "s" and quarantine "q" follow the same rule to choose when the individual begins observation
  # An individual "beta" is observed whichever is later: when "alpha" infects "beta" or when "alpha" is observed
  # Infections by alpha that occur before alpha is observed are themselves observed when alpha is observed
  # Infections by alpha after alpha is observed are immediately observed because these infections occured in a healthcare setting
  if (is.na(t_infection)==1) {cat("Error 1: t_obs_fcn ")}
  if (is.na(t_obs_alpha)==1) {cat("Error 2: t_obs_fcn ")}
  t_obs <- max(t_infection, t_obs_alpha) + CT_delay
  return(t_obs)
}

t_iso_s_fcn <- function(T_inc, T_lat, t_obs, d_symp, d_inf, epsilon, background_intervention){
  # "s" denotes symptom monitoring
  # The first generation cannot be subject to symptom monitoring
  # t_iso_s is the time an individual is isolated due to symptom monitoring
  # Isolation due to symptom monitoring will occur at at time 1 or 2, whichever is earler:
    # time 1 is the latter of: epsilon after symptom onset and time of observation
    # time 2 is either the end of disease ("u") or following "hsb"
  if (is.na(d_symp)==1) {cat("Error 1: t_iso_s_fcn ")}
  if (is.na(T_inc)==1) {cat("Error 2: t_iso_s_fcn ")}
  if (is.na(epsilon)==1) {cat("Error 3: t_iso_s_fcn ")}
  if (is.na(t_obs)==1) {cat("Error 4: t_iso_s_fcn ")}
  if (T_inc >= (T_inc + d_symp)) {cat("Error 5: t_iso_s_fcn ")}
  
  if (background_intervention == "u"){
    t_iso_s <- min( max((T_inc + epsilon), t_obs), max((T_inc + d_symp), (T_lat + d_inf)) )
  } else if (background_intervention == "hsb"){
    t_iso_s <- min( max((T_inc + epsilon), t_obs), runif(1, min = T_inc, max = (T_inc + d_symp) ) )
  }
  return(t_iso_s)
}

t_iso_q_fcn <- function(T_inc, T_lat, t_obs, d_symp, d_inf, background_intervention){
  # "q" denotes quarantine
  # t_iso_1 is the time an individual is isolated due to quarantine
  # Isolation due to quarantine will occur at at time 1 or 2, whichever is earler.
  # 1 is the latter of: either symptom or infectiousness onset and time of observation
  # 2 is the a random point of the symptomatic period due to health seeking behavior
  if (is.na(d_symp)==1) {cat("Error 1: t_iso_q_fcn")}
  if (is.na(T_inc)==1) {cat("Error 2: t_iso_q_fcn")}
  if (is.na(T_lat)==1) {cat("Error 3: t_iso_q_fcn")}
  if (is.na(t_obs)==1) {cat("Error 4: t_iso_q_fcn")}
  if (T_inc >= (T_inc + d_symp)) {cat("Error5: t_iso_q_fcn")}
  
  if (background_intervention == "u"){
    t_iso_q <- min( max(T_inc, t_obs), max((T_inc + d_symp), (T_lat + d_inf)) )
  } else if (background_intervention == "hsb"){
    t_iso_q <- min( max(T_inc, t_obs), runif(1, min = T_inc, max = (T_inc + d_symp) ) )
  }
  
  return(t_iso_q)
}

#### observe_and_isolate_fcn ####
# Create observation and isolation times for Population 
observe_and_isolate_fcn <- function(Pop, intervention){
  if (intervention == "u"){
    Pop$t_obs <- apply(Pop, 1, function(x) t_obs_u_fcn(as.numeric(x['T_inc']), as.numeric(x['d_symp']), as.numeric(x['T_lat']), as.numeric(x['d_inf'])))
    Pop$t_iso <- Pop$t_obs
  }
  if (intervention == "hsb"){
    Pop$t_obs <- round(apply(Pop, 1, function(x) t_obs_hsb_fcn(as.numeric(x['T_inc']), as.numeric(x['d_symp']))))
    Pop$t_iso <- Pop$t_obs
  }
  if (intervention == "s"){
    # The first generation cannot be subject to quarantine unless we set t_obs = 0 for all
    Pop$t_obs <- apply(Pop, 1, function(x) t_obs_fcn(as.numeric(x['t_infection']), as.numeric(x['t_obs_infector']), as.numeric(x['CT_delay'])))
    Pop$t_iso <- apply(Pop, 1, function(x) t_iso_s_fcn(as.numeric(x['T_inc']), as.numeric(x['T_lat']), as.numeric(x['t_obs']), as.numeric(x['d_symp']), as.numeric(x['d_inf']), as.numeric(x['epsilon']), x['background_intervention']))
  }
  if (intervention == "q"){
    # The first generation cannot be subject to quarantine unless we set t_obs = 0 for all
    Pop$t_obs <- apply(Pop, 1, function(x) t_obs_fcn(as.numeric(x['t_infection']), as.numeric(x['t_obs_infector']), as.numeric(x['CT_delay'])))
    Pop$t_iso <- apply(Pop, 1, function(x) t_iso_q_fcn(as.numeric(x['T_inc']), as.numeric(x['T_lat']), as.numeric(x['t_obs']), as.numeric(x['d_symp']), as.numeric(x['d_inf']), x['background_intervention']))
  }
  Pop[Pop$t_iso < Pop$t_obs,"t_obs"] <- Pop[Pop$t_iso < Pop$t_obs,"t_iso"]    # For those who were isolated due to health seeking behavior before they were observed through Q or S
  return(Pop)
}

#### pi_t_fcn ####
# Hourly R_0 considering isolation times
pi_t_fcn <- function(T_lat, d_inf, t_iso, t_obs, R_0, R_0_hsb_adjusted, gamma, distribution, triangle_center, intervention, background_intervention){   
  # Distribute the individual's R_0 across their duration of infectiousness (d_inf)
  # according to the chosen distribution of relative infectiousness (pi_t)
  # Input should be in hours
  
  if (R_0 < 0){cat("Error 1: pi_t_fcn")}
  if (d_inf <= 0){cat("Error 2: pi_t_fcn")}
  if (d_inf != round(d_inf)){cat("Error 3: pi_t_fcn")}
  if (triangle_center < 0){cat("Error 4: pi_t_fcn")}
  if (triangle_center > 1){cat("Error 5: pi_t_fcn")}
  
  # Create the unisolated distribution of pi_t
  if (distribution == "uniform"){
    pi_t <- rep(R_0 / (d_inf), d_inf)
  
  } else if (distribution == "linear_increase"){
    pi_t <- seq(1, d_inf)
    pi_t <- R_0 * (2*pi_t-1) / (d_inf)^2
  
  } else if (distribution == "triangle"){
    if (floor(d_inf*triangle_center) == 0){ # linearly decreasing
      pi_t <- seq(d_inf, 1)
      pi_t <- R_0 * (2*pi_t-1) / (d_inf)^2
    } else if (ceiling(d_inf*triangle_center) == d_inf){ # linearly increasing
      pi_t <- seq(1, d_inf)
      pi_t <- R_0 * (2*pi_t-1) / (d_inf)^2
    } else if (triangle_center > 0 & triangle_center < 1){
      pi_t.early <- seq(from=1, to=floor(d_inf*triangle_center), by=1)
      pi_t.late  <- seq(from=d_inf - floor(d_inf*triangle_center), to=1)
      pi_t.late  <- pi_t.late / ( (d_inf - floor(d_inf*triangle_center)) / floor(d_inf*triangle_center) )
      pi_t       <- R_0 * ( c(pi_t.early, pi_t.late) / sum(c(pi_t.early, pi_t.late)) )
    } else if (triangle_center < 0){ # This is a trick so that we can set it to be uniform
      pi_t <- rep(R_0 / (d_inf), d_inf) # Uniform Distribution
    }
  
  } else if (distribution == "exp_1.1"){
    pi_t <- 1.1^(1:d_inf)
    pi_t <- R_0 * pi_t / sum(pi_t)
  
  } else if (distribution == "exp_1.2"){
    pi_t <- 1.2^(1:d_inf)
    pi_t <- R_0 * pi_t / sum(pi_t)
    
  } else if (distribution == "exp_1.3"){
    pi_t <- 1.3^(1:d_inf)
    pi_t <- R_0 * pi_t / sum(pi_t)
    
  } else if (distribution == "exp_e"){
    pi_t <- exp(1:d_inf)
    pi_t <- R_0 * pi_t / sum(pi_t)
    
  } else {cat("Error 4: pi_t_fcn you must define an appropriate distribution (eg. uniform, linear_increase)")}
  
  # Now scale pi_t down by (1-gamma) after isolation (or beginning of quarantine)
  if (intervention == "q"){    # If you're under quarantine
    if (T_lat < t_obs){  # If you're infectious before you're observed under quarantine 
      if ( (T_lat + d_inf) > t_obs){ # If you are still infectious by the time you are observed under quarantine
        pi_t[seq( from = (t_obs - T_lat + 1), to = d_inf )] <- (1-gamma) * pi_t[seq( from = (t_obs - T_lat + 1), to = d_inf )]
      }
    } else if (T_lat >= t_obs){ # If you're infectious after observed under quarantine
      pi_t <- (1-gamma) * pi_t
    } else {pi_t <- pi_t}
  } else{   # If you're not under quarantine
    if (T_lat < t_iso){  # If you're infectious before you're isolated 
      if ( (T_lat + d_inf) > t_iso){ # If you are still infectious by the time you are isolated
        pi_t[seq( from = (t_iso - T_lat + 1), to = d_inf )] <- (1-gamma) * pi_t[seq( from = (t_iso - T_lat + 1), to = d_inf )]
      }
    } else if (T_lat >= t_iso){ # If you're infectious after you're isolated
      pi_t <- (1-gamma) * pi_t
    } else {pi_t <- pi_t}
  }
  
  return(pi_t)
}

#### infection_times_fcn ####
# Infection times from one individual
infection_times_fcn <- function(T_lat, d_inf, t_iso, t_obs, R_0, R_0_hsb_adjusted, gamma, distribution, triangle_center, intervention, background_intervention){
  if (sum( sum(is.na(T_lat), is.na(d_inf), is.na(R_0))) > 0){cat("Error 1: infection_times_fcn")}
  pi_t <- pi_t_fcn(T_lat, d_inf, t_iso, t_obs, R_0, R_0_hsb_adjusted, gamma, distribution, triangle_center, intervention, background_intervention)
  children <- rep(NA, length(pi_t))
  if (sum(is.na(pi_t)) > 0 | sum(pi_t < 0) > 0){cat("Error 2: infection_times_fcn")}
  children <- sapply(pi_t, function(x) rpois(1, lambda = x))
  return(children)
}

#### children_list_fcn ####
# Create a list of the children of the population
children_list_fcn <- function(Pop, pi_t_distribution, triangle_center, gamma, intervention, background_intervention){
  children_list <- as.list(seq(1:nrow(Pop)))
  children_list <- apply(Pop, 1, function(x) list(infection_times_fcn(as.numeric(x['T_lat']), as.numeric(x['d_inf']), as.numeric(x['t_iso']), as.numeric(x['t_obs']), as.numeric(x['R_0']),as.numeric(x['R_0_hsb_adjusted']), gamma, pi_t_distribution, triangle_center, intervention, background_intervention)))
  children_list <- lapply(children_list, "[[", 1)  #This removes one layer of [[ ]] from the list
  return(children_list)
}

#### next_generation_fcn ####
# Create the next generation
next_generation_fcn <- function(Pop,
                                children_list,
                                parms_T_inc,
                                parms_T_lat,
                                parms_d_inf,
                                parms_d_symp,
                                parms_R_0,
                                parms_epsilon,
                                generation,
                                prob_CT,
                                parms_CT_delay,
                                gamma,
                                n_pop){
  n_pop_next <- sum(unlist(lapply(children_list, sum)))
  Pop_2 <- Create_Pop(n_pop_next , parms_T_inc, parms_T_lat, parms_d_inf, parms_d_symp, parms_R_0, parms_epsilon, generation, Pop[1,"background_intervention"], parms_CT_delay, gamma)
  index = 1
  for (i in 1:length(children_list)){     # for each list of children
    if (index <= n_pop){
      if (sum(children_list[[i]]==1)>0){    # if there is at least 1 onward infection by person i
        count <- sum(children_list[[i]]==1)   
        for (j in seq(from=0, to=count-1)){           # for each child j infected by person i
          Pop_2[index+j, "infector"] <- i
          Pop_2[index+j, "t_obs_infector"] <- Pop[i,"t_obs"]
          Pop_2[index+j, "t_infection"] <- which(children_list[[i]]==1)[j+1] + Pop[i, "T_lat"]     
        }
        index <- index + count     
      }
    }
  }
  Pop_2 <- Pop_2[is.na(Pop_2$T_inc)==0,]
  Pop_2 <- Pop_2[is.na(Pop_2$d_symp)==0,]
  Pop_2 <- Pop_2[is.na(Pop_2$infector)==0,]
  Pop_2 <- Pop_2[is.na(Pop_2$t_infection)==0,]
  Pop_2[,"T_inc"] <- ( Pop_2[,"T_inc"] + Pop_2[,"t_infection"] )
  Pop_2[,"T_lat"] <- ( Pop_2[,"T_lat"] + Pop_2[,"t_infection"] )
  if (prob_CT < 1){
    identified <- rbinom(length(Pop_2$t_obs_infector), 1, prob = prob_CT)
    Pop_2[which(identified==0), "t_obs_infector"] <- 999999999
  }
  return(Pop_2)
}

#### repeat_call_fcn ####
repeat_call_fcn <- function(n_pop, 
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
{
  input <- list(n_pop, parms_T_inc, 
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
  names(input) <- c("n_pop", 
                    "parms_T_inc", 
                    "parms_T_lat", 
                    "parms_d_inf", 
                    "parms_d_symp", 
                    "parms_R_0", 
                    "parms_epsilon",
                    "parms_pi_t",
                    "num_generations",
                    "background_intervention",
                    "subseq_interventions",
                    "gamma",
                    "prob_CT",
                    "parms_CT_delay",
                    "parms_serial_interval")
  
  Pop_1 <- Create_Pop(n_pop = n_pop, 
                      parms_T_inc = parms_T_inc, 
                      parms_T_lat = parms_T_lat, 
                      parms_d_inf = parms_d_inf, 
                      parms_d_symp = parms_d_symp, 
                      parms_R_0 = parms_R_0, 
                      parms_epsilon = parms_epsilon, 
                      generation = 1,
                      background_intervention,
                      parms_CT_delay,
                      gamma)
  Pop_1 <- observe_and_isolate_fcn(Pop = Pop_1, intervention = background_intervention)
  children_list <- children_list_fcn(Pop = Pop_1, pi_t_distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma = gamma, intervention = background_intervention, background_intervention = background_intervention)
  Num_Infected <- unlist(lapply(children_list, sum))
  cat('Generation 1 : n=', nrow(Pop_1), 'Effective Reproductive Number:', mean(Num_Infected))
  
  names <- c("n", "R", "obs_to_iso", "prop_lat_before_obs", "ks")
  output <- data.frame(matrix(rep(NA, num_generations*length(names)), nrow=num_generations))
  names(output) <- names
  output[1,"n"] <- nrow(Pop_1)
  output[1,"R"] <- mean(Num_Infected)
  output[1,"obs_to_iso"] <- mean(Pop_1$t_iso - Pop_1$t_obs)
  output[1,"prop_lat_before_obs"] <- mean((Pop_1$T_lat < Pop_1$t_obs) == 1)
  
  Pop_prev <- Pop_1
  
  if (num_generations > 1){
    for (g in seq(from=2, to=num_generations)){
      if (sum(Num_Infected) > 20){
        Pop_next <- NA
        Pop_next <- next_generation_fcn(Pop = Pop_prev,
                                        children_list,
                                        parms_T_inc = parms_T_inc,
                                        parms_T_lat = parms_T_lat,
                                        parms_d_inf = parms_d_inf,
                                        parms_d_symp = parms_d_symp,
                                        parms_R_0 = parms_R_0,
                                        parms_epsilon = parms_epsilon,
                                        generation = g,
                                        prob_CT = prob_CT,
                                        parms_CT_delay,
                                        gamma,
                                        n_pop)

        if (nrow(Pop_next) > (n_pop)){   #don't let the populations get bigger than the initial
          Pop_next <- Pop_next[1:(n_pop),]
        }
        
        Pop_next <- observe_and_isolate_fcn(Pop = Pop_next, intervention = subseq_interventions)
        
        children_list <- NA
        children_list <- children_list_fcn(Pop = Pop_next, pi_t_distribution = parms_pi_t$distribution, triangle_center = parms_pi_t$triangle_center, gamma = gamma, intervention = subseq_interventions, background_intervention)
        
        Num_Infected <- NA
        Num_Infected <- unlist(lapply(children_list, sum))
        
        cat('\nGeneration',g, ': n=', nrow(Pop_next), 'Effective Reproductive Number:', mean(Num_Infected))
        output[g,"n"] <- nrow(Pop_next)
        output[g,"R"] <- mean(Num_Infected)
        output[g,"obs_to_iso"] <- mean(Pop_next$t_iso - Pop_next$t_obs)
        output[g,"prop_lat_before_obs"] <- mean((Pop_next$T_lat < Pop_next$t_obs) == 1)
        if (parms_serial_interval$dist != "unknown" & subseq_interventions == "u"){
          output[g,"ks"] <- serial_interval_fcn(Pop_prev, Pop_next, parms_serial_interval, plot="False")
        }
        
        Pop_prev <- NA
        Pop_prev <- Pop_next
        
      } else {cat('\nPopulation size dropped below 20')}
    }
  }
  cat('\n\n')
  
  output <- output[is.na(output$n)==0,]
  In_Out <- list(input, output)
  names(In_Out) <- c("input","output")
  return(In_Out)
}

##### serial_interval_fcn ####
# Generate Serial Interval Distribution
serial_interval_fcn <- function(Pop1, Pop2, parms_serial_interval, plot=TRUE){
  output <- rep(NA, length(Pop2$ID))
  for (i in Pop2$ID){
    output[i] <- jitter(Pop2[i,"T_inc"] - Pop1[Pop1$ID == Pop2[i,"infector"], "T_inc"], amount = 1) #jitter to prevent ties
    if (output[i] < 24){output[i] <- runif(1, min=1, max=24)} # If the serial interval is less than one day, we're really unsure of what hour it was
  }
  if (sum(is.na(output)==1)>0){cat("Error 1: serial_interval_fcn")}
  
  if (plot=="True"){
    #Need a multiplier for the curves to plot on the same axes
    multipliers <- (hist(output/24, breaks=seq(min(output/24), max(output/24)+1), plot=FALSE)$counts / hist(output/24,breaks=seq(min(output/24), max(output/24)+1), plot=FALSE)$density)
    multiplier <- multipliers[is.na(multipliers)==0][1]
    
    hist(output/24, breaks=seq(min(output/24), max(output/24)+1), density=10, 
         xlab = "Serial Interval (days)",
         main = "")
    
    if (parms_serial_interval$dist == "gamma"){
      fit <- fitdistr(output/24, "gamma", lower=0.00001, start = list(shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))
      curve(multiplier*dgamma(x, shape = as.numeric(fit$estimate[1]), rate = as.numeric(fit$estimate[2])), 
            from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="darkblue", lwd=2)
      curve(multiplier*dgamma(x, shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2),
            from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="green", lwd=2)
      ks <- ks.test(x=output/24, "pgamma", shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2)
    }
    if (parms_serial_interval$dist == "weibull"){
      fit <- fitdistr(output/24, "weibull", lower=0.00001, start = list(shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2))
      curve(multiplier*dweibull(x, shape = as.numeric(fit$estimate[1]), scale = as.numeric(fit$estimate[2])), 
            from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="darkblue", lwd=2)
      curve(multiplier*dweibull(x, shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2),
            from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="green", lwd=2)
      ks <- ks.test(x=output/24, "pweibull", shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2)
    }
    if (parms_serial_interval$dist == "lognormal"){
      fit <- fitdistr(output/24, "lognormal", lower=0.00001)
      curve(multiplier*dlnorm(x, meanlog = as.numeric(fit$estimate[1]), sdlog = as.numeric(fit$estimate[2])), 
            from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="darkblue", lwd=2)
      curve(multiplier*dlnorm(x, meanlog= log(parms_serial_interval$parm1), sdlog=log(parms_serial_interval$parm2)),
            from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="green", lwd=2)
      ks <- ks.test(x=output/24, "plnorm", meanlog=log(parms_serial_interval$parm1), sdlog=log(parms_serial_interval$parm2))
    }
    
  } else if (plot=="Add"){
    #Need a multiplier for the curves to plot on the same axes
    multipliers <- (hist(output/24, breaks=seq(min(output/24), max(output/24)+1), plot=FALSE)$counts / hist(output/24,breaks=seq(min(output/24), max(output/24)+1), plot=FALSE)$density)
    multiplier <- multipliers[is.na(multipliers)==0][1]
    
    if (parms_serial_interval$dist == "gamma"){
      fit <- fitdistr(output/24, "gamma", lower=0.00001, start = list(shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2))
      suppressWarnings(curve(multiplier*dgamma(x, shape = as.numeric(fit$estimate[1]), rate = as.numeric(fit$estimate[2])), 
                             from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="darkblue", lwd=0.5, 
                             ylab="Empirical Distribution", xlab="Days", main="Serial Intervals"))
      ks <- ks.test(x=output/24, "pgamma", shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2)
    }
    if (parms_serial_interval$dist == "weibull"){
      fit <- fitdistr(output/24, "weibull", lower=0.00001, start = list(shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2))
      suppressWarnings(curve(multiplier*dweibull(x, shape = as.numeric(fit$estimate[1]), scale = as.numeric(fit$estimate[2])), 
                             from=0, to=max(output/24)+1, add=TRUE, yaxt="n", col="darkblue", lwd=0.5, 
                             ylab="Empirical Distribution", xlab="Days", main="Serial Intervals"))
      ks <- ks.test(x=output/24, "pweibull", shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2)
    }
    
  } else if (plot == "False"){
    #Need a multiplier for the curves to plot on the same axes    
    multipliers <- (hist(output/24, breaks=seq(min(output/24), max(output/24)+1), plot=FALSE)$counts / hist(output/24,breaks=seq(min(output/24), max(output/24)+1), plot=FALSE)$density)
    multiplier <- multipliers[is.na(multipliers)==0][1]
    
    if (parms_serial_interval$dist == "gamma"){
      ks <- ks.test(x=output/24, "pgamma", shape=parms_serial_interval$parm1, rate=parms_serial_interval$parm2, exact=FALSE)
    }
    if (parms_serial_interval$dist == "weibull"){
      ks <- ks.test(x=output/24, "pweibull", shape=parms_serial_interval$parm1, scale=parms_serial_interval$parm2, exact=FALSE)
    }
  }
  return(as.numeric(ks$statistic))
}

#### calculate_R_0_hsb ####
# Calculate R_O_hsb from a desired R_0 distribution in a setting where isolation occurs
# calculate_R_0_hsb <- function(Pop, desired_R_0,  
#                               parms_T_inc, 
#                               parms_T_lat, 
#                               parms_d_inf, 
#                               parms_d_symp,
#                               parms_epsilon, 
#                               num_generations,
#                               pi_t_distribution,
#                               gamma,
#                               prob_CT,
#                               parms_CT_delay,
#                               times,
#                               num_generations = 5){
#   
#   n_pop = nrow(Pop)
#   names <- c("R_0", "R_hsb", "R_s", "R_q", "Abs_Benefit","Rel_Benefit","NNQ","obs_to_iso_q","Abs_Benefit_per_Qday")
#   data <- data.frame(matrix(rep(NA, length(names)*times), nrow=times))
#   names(data) <- names
#   dimensions <- c("R_0")
#   
#   lhs <- maximinLHS(times, length(dimensions))
#   
#   R_0.min <- desired_R_0
#   R_0.max <- 40
#   
#   params.set <- cbind(
#     gamma = lhs[,1]*(gamma.max - gamma.min) + gamma.min,
#     prob_CT = lhs[,2]*(prob_CT.max - prob_CT.min) + prob_CT.min,
#     CT_delay = lhs[,3]*(CT_delay.max - CT_delay.min) + CT_delay.min,
#     epsilon = lhs[,4]*(epsilon.max - epsilon.min) + epsilon.min,
#     R_0 = lhs[,5]*(R_0.max - R_0.min) + R_0.min,
#     pi_t_distribution = lhs[,6]*(pi_t_distribution.max - pi_t_distribution.min) + pi_t_distribution.min)
#   
#   for (i in 1:times){
#     cat('\nIteration',i, '\n')
#     
#     gamma <- params.set[i,"gamma"]
#     prob_CT <- params.set[i,"prob_CT"]
#     parms_CT_delay$parm2 <- params.set[i,"CT_delay"]
#     parms_epsilon$parm2 <- params.set[i,"epsilon"]
#     
#     # Spread of R_0 is always 0.5. The mean is changing
#     parms_R_0[c("parm1","parm2","parm3")] <- c(params.set[i,"R_0"] - 0.25, params.set[i,"R_0"] + 0.25, params.set[i,"R_0"])
#     
#     if (params.set[i,"pi_t_distribution"] == 1){pi_t_distribution <- "linear_increase"}
#     if (params.set[i,"pi_t_distribution"] == 0){pi_t_distribution <- "triangle"}
#     
#     for (subseq_interventions in c(background_intervention, "hsb", "s","q")){      
#       In_Out <- repeat_call_fcn(n_pop=n_pop, 
#                                 parms_T_inc, 
#                                 parms_T_lat, 
#                                 parms_d_inf, 
#                                 parms_d_symp, 
#                                 parms_R_0, 
#                                 parms_epsilon, 
#                                 num_generations,
#                                 background_intervention,
#                                 subseq_interventions,
#                                 pi_t_distribution,
#                                 gamma,
#                                 prob_CT,
#                                 parms_CT_delay)
#       if (subseq_interventions == background_intervention){
#         data[i,"R_0"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
#       }
#       if (subseq_interventions == "hsb"){
#         data[i,"R_hsb"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
#       }
#       if (subseq_interventions == "s"){
#         data[i,"R_s"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
#       }
#       if (subseq_interventions == "q"){
#         data[i,"R_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"R"], w=In_Out$output[3:nrow(In_Out$output),"n"])
#         data[i,"obs_to_iso_q"] <- weighted.mean(x=In_Out$output[3:nrow(In_Out$output),"obs_to_iso"], w=In_Out$output[3:nrow(In_Out$output),"n"]) / 24
#       }
#     }
#   }
# }

#### decile_plot_fcn ####
decile_plot_fcn <- function(data, params.set){
  require(Hmisc)
  params <- names(data.frame(params.set))
  if (is.element(paste(params[1], "_quantile",sep=""), names(data))){
    for (col in params){
      col_name <- paste(col, "_quantile",sep="")
      data[,col_name] <- as.numeric(cut2(data[,col], g=10, levels.mean=TRUE))
    }
  } else {
    for (col in params){
      newcol <- as.numeric(cut2(data[,col], g=10, levels.mean=TRUE))
      data <- cbind(data, newcol)
      col_name <- paste(col, "_quantile",sep="")
      names <- c(names(data)[1:ncol(data)-1], col_name)
      names(data) <- names
    }
  }
  quant_vars <- names(data)[grep("quantile", names(data))]
  data_melt <- melt.data.frame(data, id.vars = c("ks"), measure.vars = quant_vars)
  ggplot(data_melt, aes(x=as.factor(value), y=ks)) + geom_boxplot() + facet_grid(.~variable) + theme_bw() + xlab("Decile")
}

#### multiplot ####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

