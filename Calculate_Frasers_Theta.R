#### Header ####
# Calculate Fraser's Theta for a given set of disease characteristics
# Borrowed heavily from "Figure_Fraster2004.R" on Dec 29, 2015

#### Load Workspace ####
root <- "20151022_SARS"
root <- "20151024_Ebola"
root <- "20151026_HepatitisA"
root <- "20151026_Pertussis"
root <- "20151027_MERS"
root <- "20151028_InfluenzaA"
root <- "20151028_Smallpox"
root <- "20151104_Ebola"
root <- "20151104_InfluenzaA"

load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "/", root, "_SMC.RData", sep=""))

head(data)

#### Calculate Fraser's Theta ####
# Theta is the proprtion of infections that occur prior to symptom onset
data$theta <- NA

for (i in 1:nrow(data)){
  a <- 0
  b <- data[i, "d_inf"]
  c <- data[i, "pi_t_triangle_center"] * data[i, "d_inf"]
  x <- -(data[i, "T_lat_offset"])   # Time of symptom onset
  
  if (data[i, "T_lat_offset"] >= 0){   # If infectious after symptom onset
    data[i, "theta"] = 0
  } else if (b <= x){   # If end of infectiousness is before symptom onset
    data[i, "theta"] = 1
  } else {
    if (c <= x){   # If the peak occurs before onset of symptoms
      # Calculate the right hand side triangle
      height <- 2*(b-x) / ((b-a)*(b-c))
      base <- b-x
      rhs <- 0.5*base*height
      data[i, "theta"] <- 1-rhs
    }
    if (c > x){
      # Calculate the left hand side triangle
      height <- 2*(x-a) / ((b-a)*(c-a))
      base <- x
      lhs <- 0.5*base*height
      data[i, "theta"] <- lhs
    }
  }
}

#### Summary statistics for Fraser Theta ####
cat("The median theta is", median(data$theta))
cat("The 95% CI for theta is [", sort(data$theta)[length(data$theta)*0.025], ", ", sort(data$theta)[length(data$theta)*0.975], "]", sep = "")
cat("The mean theta is", mean(data$theta))

hist(data$theta, breaks = 20, xlim = c(0,1))

#### Plots with Theta as Y ####
plot(data$T_lat_offset, data$theta)
plot(data$d_inf, data$theta)
plot(data$pi_t_triangle_center, data$theta)
plot(data$ks,data$theta)

#### Save Workspace ####
save.image(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", root, "/", root, "_SMC.RData", sep=""))
 

