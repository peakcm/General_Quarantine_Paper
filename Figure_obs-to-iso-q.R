#### Header ####
# Plot the obs_to_iso_q for each disease to see how long we are following them before they develop symptoms

#### Load Libraries ####
library(ggplot2)

#### Load Workspaces ####
desired_date <- "20151111"
# load(paste("~/Dropbox/Ebola/General_Quarantine_Paper/General_Quarantine_Paper/", desired_date, "_PRCC.RData", sep=""))

#### Plot ####
ggplot(df.prcc, aes(x=obs_to_iso_q, fill = disease)) + facet_grid(disease~.) + geom_bar()

ggplot(df.prcc[df.prcc$prob_CT < .6,], aes(x=obs_to_iso_q, fill = disease)) + facet_grid(disease~.) + geom_bar()

ggplot(df.prcc[df.prcc$prob_CT > .6,], aes(x=obs_to_iso_q, fill = disease)) + facet_grid(disease~.) + geom_bar()

ggplot(df.prcc, aes(x=obs_to_iso_q, fill = (prob_CT>0.6))) + facet_grid(disease~.) + geom_bar(position = "dodge")
