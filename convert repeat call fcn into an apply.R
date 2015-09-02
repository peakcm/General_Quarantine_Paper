params.set.old <- params.set
params.set <- data.frame(params.set)
names.old <- names(params.set)
names.new <- c("anchor_value", "parm1", "parm2", "triangle_center")
names(params.set) <- names.new

params.set[,"parm2"] <- params.set[,"parm1"] + params.set[,"parm2"]

start.time <- proc.time()
out <- apply(params.set, 1, function(x) repeat_call_fcn(n_pop, 
                                                 parms_T_inc, 
                                                 parms_T_lat = c(parms_T_lat[1:4], x["anchor_value"], parms_T_lat[6]), 
                                                 parms_d_inf, 
                                                 parms_d_symp = c(parms_d_symp[1], x["parm1"], x["parm2"],parms_d_symp[4:6]), 
                                                 parms_R_0, 
                                                 parms_epsilon, 
                                                 parms_pi_t = c(parms_pi_t[1], x["triangle_center"]),
                                                 num_generations,
                                                 background_intervention,
                                                 subseq_interventions,
                                                 gamma,
                                                 prob_CT,
                                                 parms_CT_delay,
                                                 parms_serial_interval,
                                                 printing = printing))
proc.time() - start.time
out[[1]]$output
