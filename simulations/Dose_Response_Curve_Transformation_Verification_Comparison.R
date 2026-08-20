################################################################################################################################
#                                                                                                                              #
#      Verification of the design-sensitivity formula with the kappa transformation                                            #
#      and the computation of the design sensitivities under different settings                                                #
#      for Section C.3 of the supplement, Table S.7 and Table S.8                                                              #
#                                                                                                                              #
#      DGP: Square treatment effect, and log DGP                                                                               #
#      Kappa choices: identity, z^2 , log(z)                                                                                   #
#                                                                                                                              #                                                               #                                                                      #
#                                                                                                                              #
################################################################################################################################



library(sensitivitymv)
library(mvtnorm)
source("matched_sensitivity_analysis_functions.R")




################################################################################
#                                                                              #
#               Define the three kappa choices                                 #
#                                                                              #
################################################################################


kappa.identity = function(z) z
kappa.square   = function(z) z^2
kappa.log      = function(z) log(z)




################################################################################
#                                                                              #
#       Step 1: Design sensitivity under each kappa, Square DGP                #
#                                                                              #
################################################################################


set.seed(100)
population_data = generate_polynomial_effect(I=500000, power=2)


## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.square, design.log)
rownames(design.kappa.tb) = c("kappa=identity","kappa=z^2","kappa=log")

design.kappa.tb




################################################################################
#                                                                              #
#       Step 2: Simulated power under each kappa, square DGP                   #
#                                                                              #
################################################################################


## simulated_power_polynomial_2 is defined in
## Continuous_Outcome_Generalized_Design_Bahadur_Efficiency.R
## and accepts the kappa argument (defaults to identity). We copy-paste it here

simulated_power_polynomial_2 = function(method=c("wilcoxon","dose.weighted","U","dose.weighted.U"),rank.order = 2, m=2,m1=2,m2=2,Gamma_seq, alpha=0.01, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_polynomial_effect(I=5000,seed=seed, power = 2)
      p_value = sendose(x=sample_data$dose,y=sample_data$response, method=method,alternative = "greater than",rank.order=rank.order,m=m,m1=m1,m2=m2,
                        Gamma_ave = Gamma_seq[i], kappa = kappa)$pval
      p.value[seed] = p_value
    }
    simulated.power = mean(p.value<alpha)
    power_seq[i] = simulated.power
    cat("current Gamma", Gamma_seq[i],"with power",simulated.power,"\n")
  }
  return(list(method=method,Gamma_seq=Gamma_seq,simulated.power = power_seq, alpha=alpha))
}



Gamma_seq = c(1.5,2,2.5,3,3.5,4)


## kappa = identity (baseline; matches existing Table 2 numbers)
wilcoxon.power.id      = simulated_power_polynomial_2(method="wilcoxon",        Gamma_seq=Gamma_seq, kappa=kappa.identity)
dose.weighted.power.id = simulated_power_polynomial_2(method="dose.weighted",   Gamma_seq=Gamma_seq, kappa=kappa.identity)
U878.power.id          = simulated_power_polynomial_2(method="U",               m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.identity)
dose.U878.power.id     = simulated_power_polynomial_2(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.identity)


## and compute the power for the adaptive (wilcoxon, dose-weighted wilcoxon)

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.1.id = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    poly.data = generate_polynomial_effect(I = 5000, power = 2, seed=run)
    res = sendose.adaptive(
      x = poly.data$dose, 
      y = poly.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.identity, kappa2 = kappa.identity
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.1.id[index] = power
}



## and compute the power for the adaptive (U, dose-weighted U)
Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.2.id = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    poly.data = generate_polynomial_effect(I = 5000, power = 2, seed=run)
    res = sendose.adaptive(
      x = poly.data$dose, 
      y = poly.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.identity, kappa2 = kappa.identity
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.2.id[index] = power
}






power.identity = rbind(Gamma_seq,
                       wilcoxon.power.id$simulated.power,
                       dose.weighted.power.id$simulated.power, power.adaptive.1.id,
                       U878.power.id$simulated.power,
                       dose.U878.power.id$simulated.power, power.adaptive.2.id)
rownames(power.identity) = c("Gamma","Wilcoxon","dose-wt","adaptive-1","(8,7,8)","dose.(8,7,8)","adaptive-2")


## Reproduced Table S.7 in Section C.3 of the supplementary material, the second 
## block of Dose-response Square and \kappa(z) = z
power.identity





## kappa = z^2
wilcoxon.power.sq      = simulated_power_polynomial_2(method="wilcoxon",        Gamma_seq=Gamma_seq, kappa=kappa.square)
dose.weighted.power.sq = simulated_power_polynomial_2(method="dose.weighted",   Gamma_seq=Gamma_seq, kappa=kappa.square)
U878.power.sq          = simulated_power_polynomial_2(method="U",               m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.square)
dose.U878.power.sq     = simulated_power_polynomial_2(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.square)


## and compute the power for the adaptive (wilcoxon, dose-weighted wilcoxon)

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.1.sq = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    poly.data = generate_polynomial_effect(I = 5000, power = 2, seed=run)
    res = sendose.adaptive(
      x = poly.data$dose, 
      y = poly.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.square, kappa2 = kappa.square
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.1.sq[index] = power
}



## and compute the power for the adaptive (U, dose-weighted U)
Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.2.sq = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    poly.data = generate_polynomial_effect(I = 5000, power = 2, seed=run)
    res = sendose.adaptive(
      x = poly.data$dose, 
      y = poly.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.square, kappa2 = kappa.square
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.2.sq[index] = power
}





power.square = rbind(Gamma_seq,
                     wilcoxon.power.sq$simulated.power,
                     dose.weighted.power.sq$simulated.power, power.adaptive.1.sq,
                     U878.power.sq$simulated.power,
                     dose.U878.power.sq$simulated.power, power.adaptive.2.sq)
rownames(power.square) = c("Gamma","Wilcoxon","dose-wt","adaptive-1","(8,7,8)","dose.(8,7,8)","adaptive-2")


## Reproduced Table S.7 Square dose-response curve with \kappa(z) = z^2
power.square


## kappa = log(z)
wilcoxon.power.log      = simulated_power_polynomial_2(method="wilcoxon",        Gamma_seq=Gamma_seq, kappa=kappa.log)
dose.weighted.power.log = simulated_power_polynomial_2(method="dose.weighted",   Gamma_seq=Gamma_seq, kappa=kappa.log)
U878.power.log          = simulated_power_polynomial_2(method="U",               m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.log)
dose.U878.power.log     = simulated_power_polynomial_2(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.log)


## and compute the power for the adaptive (wilcoxon, dose-weighted wilcoxon)

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.1.log = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    poly.data = generate_polynomial_effect(I = 5000, power = 2, seed=run)
    res = sendose.adaptive(
      x = poly.data$dose, 
      y = poly.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.log, kappa2 = kappa.log
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.1.log[index] = power
}



## and compute the power for the adaptive (U, dose-weighted U)
Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.2.log = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    poly.data = generate_polynomial_effect(I = 5000, power = 2, seed=run)
    res = sendose.adaptive(
      x = poly.data$dose, 
      y = poly.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.log, kappa2 = kappa.log
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.2.log[index] = power
}




power.log = rbind(Gamma_seq,
                  wilcoxon.power.log$simulated.power,
                  dose.weighted.power.log$simulated.power,
                  power.adaptive.1.log,
                  U878.power.log$simulated.power,
                  dose.U878.power.log$simulated.power, power.adaptive.2.log)
rownames(power.log) = c("Gamma","Wilcoxon","dose-wt","adaptive-1","(8,7,8)","dose.(8,7,8)", "adaptive-2")

## Reproduced Table S.7 Square dose-response curve and \kappa(z) = log(z)
power.log


###################################################################
#    log DGP
###################################################################


set.seed(100)
population_data = generate_log_effect(I=500000)


## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb.log.dgp = rbind(design.identity, design.square, design.log)
rownames(design.kappa.tb.log.dgp) = c("kappa=identity","kappa=z^2","kappa=log")


## Table S.7, design sensitivities
design.kappa.tb.log.dgp




################################################################################
#                                                                              #
#       Step 2: Simulated power under each kappa, log DGP                      #
#                                                                              #
################################################################################


simulated_power_log = function(method=c("wilcoxon","dose.weighted","U","dose.weighted.U"),rank.order = 2, m=2,m1=2,m2=2,Gamma_seq, alpha=0.01, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_log_effect(I=5000,seed=seed)
      p_value = sendose(x=sample_data$dose,y=sample_data$response, method=method,alternative = "greater than",rank.order=rank.order,m=m,m1=m1,m2=m2,
                        Gamma_ave = Gamma_seq[i], kappa = kappa)$pval
      p.value[seed] = p_value
    }
    simulated.power = mean(p.value<alpha)
    power_seq[i] = simulated.power
    cat("current Gamma", Gamma_seq[i],"with power",simulated.power,"\n")
  }
  return(list(method=method,Gamma_seq=Gamma_seq,simulated.power = power_seq, alpha=alpha))
}



Gamma_seq = c(1.5,2,2.5,3,3.5,4)


## kappa = identity
wilcoxon.power.id      = simulated_power_log(method="wilcoxon",        Gamma_seq=Gamma_seq, kappa=kappa.identity)
dose.weighted.power.id = simulated_power_log(method="dose.weighted",   Gamma_seq=Gamma_seq, kappa=kappa.identity)
U878.power.id          = simulated_power_log(method="U",               m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.identity)
dose.U878.power.id     = simulated_power_log(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.identity)


## and compute the power for the adaptive (wilcoxon, dose-weighted wilcoxon)

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.1.id = rep(NA, length(Gamma_vec))

for(index in 1:length(Gamma_vec)){
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data
    log.data = generate_log_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = log.data$dose,
      y = log.data$response,
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.identity, kappa2 = kappa.identity
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.1.id[index] = power
}



## and compute the power for the adaptive (U, dose-weighted U)
Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.2.id = rep(NA, length(Gamma_vec))

for(index in 1:length(Gamma_vec)){
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data
    log.data = generate_log_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = log.data$dose,
      y = log.data$response,
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.identity, kappa2 = kappa.identity
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.2.id[index] = power
}




power.identity.log = rbind(Gamma_seq,
                           wilcoxon.power.id$simulated.power,
                           dose.weighted.power.id$simulated.power, power.adaptive.1.id,
                           U878.power.id$simulated.power,
                           dose.U878.power.id$simulated.power, power.adaptive.2.id)
rownames(power.identity.log) = c("Gamma","Wilcoxon","dose-wt","adaptive-1","(8,7,8)","dose.(8,7,8)","adaptive-2")



## Table S.7 Log dose-response curve, \kappa(z) = z
power.identity.log





## kappa = z^2
wilcoxon.power.sq      = simulated_power_log(method="wilcoxon",        Gamma_seq=Gamma_seq, kappa=kappa.square)
dose.weighted.power.sq = simulated_power_log(method="dose.weighted",   Gamma_seq=Gamma_seq, kappa=kappa.square)
U878.power.sq          = simulated_power_log(method="U",               m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.square)
dose.U878.power.sq     = simulated_power_log(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.square)


## and compute the power for the adaptive (wilcoxon, dose-weighted wilcoxon)

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.1.sq = rep(NA, length(Gamma_vec))

for(index in 1:length(Gamma_vec)){
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data
    log.data = generate_log_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = log.data$dose,
      y = log.data$response,
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.square, kappa2 = kappa.square
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.1.sq[index] = power
}



## and compute the power for the adaptive (U, dose-weighted U)
Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.2.sq = rep(NA, length(Gamma_vec))

for(index in 1:length(Gamma_vec)){
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data
    log.data = generate_log_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = log.data$dose,
      y = log.data$response,
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.square, kappa2 = kappa.square
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.2.sq[index] = power
}





power.square.log = rbind(Gamma_seq,
                         wilcoxon.power.sq$simulated.power,
                         dose.weighted.power.sq$simulated.power, power.adaptive.1.sq,
                         U878.power.sq$simulated.power,
                         dose.U878.power.sq$simulated.power, power.adaptive.2.sq)
rownames(power.square.log) = c("Gamma","Wilcoxon","dose-wt","adaptive-1","(8,7,8)","dose.(8,7,8)","adaptive-2")


## Reproduced Table S.7 log dose-response curve, \kappa(z) = z^2
power.square.log


## kappa = log(z)
wilcoxon.power.log      = simulated_power_log(method="wilcoxon",        Gamma_seq=Gamma_seq, kappa=kappa.log)
dose.weighted.power.log = simulated_power_log(method="dose.weighted",   Gamma_seq=Gamma_seq, kappa=kappa.log)
U878.power.log          = simulated_power_log(method="U",               m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.log)
dose.U878.power.log     = simulated_power_log(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=Gamma_seq, kappa=kappa.log)


## and compute the power for the adaptive (wilcoxon, dose-weighted wilcoxon)

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.1.log = rep(NA, length(Gamma_vec))

for(index in 1:length(Gamma_vec)){
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data
    log.data = generate_log_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = log.data$dose,
      y = log.data$response,
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.log, kappa2 = kappa.log
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.1.log[index] = power
}



## and compute the power for the adaptive (U, dose-weighted U)
Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.adaptive.2.log = rep(NA, length(Gamma_vec))

for(index in 1:length(Gamma_vec)){
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data
    log.data = generate_log_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = log.data$dose,
      y = log.data$response,
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave, kappa1 = kappa.log, kappa2 = kappa.log
    )
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 100 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.adaptive.2.log[index] = power
}




power.log.log = rbind(Gamma_seq,
                      wilcoxon.power.log$simulated.power,
                      dose.weighted.power.log$simulated.power,
                      power.adaptive.1.log,
                      U878.power.log$simulated.power,
                      dose.U878.power.log$simulated.power, power.adaptive.2.log)
rownames(power.log.log) = c("Gamma","Wilcoxon","dose-wt","adaptive-1","(8,7,8)","dose.(8,7,8)", "adaptive-2")


## Reproduced Table S.7 log dose-response curve, \kappa(z) = log(z)
power.log.log




#######################################################################################################################
#
#
#   Compute the design sensitivity across different choice of kappas (dose-transformation) and different tests
#   under six data generating processes, reproduced Table S.8 in the supplementary material
#
#######################################################################################################################



kappa.identity = function(z) z
kappa.log      = function(z) log(z)
kappa.square   = function(z) z^2
kappa.sqrt     = function(z) sqrt(z)



###### The design sensitivity under the Square DGP


set.seed(100)
population_data = generate_polynomial_effect(I=500000, power=2)


## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")





## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = sqrt(z)

wilcoxon.design.sqrt      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.sqrt)
dose.weighted.design.sqrt = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.sqrt)
U878.design.sqrt          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.sqrt)
dose.U878.design.sqrt     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.sqrt)
design.sqrt = c(wilcoxon.design.sqrt, dose.weighted.design.sqrt, U878.design.sqrt, dose.U878.design.sqrt)
names(design.sqrt) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")



## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.log, design.square, design.sqrt)
rownames(design.kappa.tb) = c("kappa=identity","kappa=log","kappa=z^2", "kappa = sqrt(z)")



## Table S.8 the Square dose-response curve
design.kappa.tb




############ Kink DGP


set.seed(100)
population_data = generate_kink_effect(I=500000, linear_coefficient = 1.5, kink=0.8)


## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")





## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = sqrt(z)

wilcoxon.design.sqrt      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.sqrt)
dose.weighted.design.sqrt = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.sqrt)
U878.design.sqrt          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.sqrt)
dose.U878.design.sqrt     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.sqrt)
design.sqrt = c(wilcoxon.design.sqrt, dose.weighted.design.sqrt, U878.design.sqrt, dose.U878.design.sqrt)
names(design.sqrt) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")



## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.log, design.square, design.sqrt)
rownames(design.kappa.tb) = c("kappa=identity","kappa=log","kappa=z^2", "kappa = sqrt(z)")



## Table S.8, Kink dose-response curve
design.kappa.tb




############ Linear DGP


set.seed(100)
population_data = generate_linear_effect(I=500000, linear_coefficient =1)


## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")





## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = sqrt(z)

wilcoxon.design.sqrt      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.sqrt)
dose.weighted.design.sqrt = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.sqrt)
U878.design.sqrt          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.sqrt)
dose.U878.design.sqrt     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.sqrt)
design.sqrt = c(wilcoxon.design.sqrt, dose.weighted.design.sqrt, U878.design.sqrt, dose.U878.design.sqrt)
names(design.sqrt) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")



## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.log, design.square, design.sqrt)
rownames(design.kappa.tb) = c("kappa=identity","kappa=log","kappa=z^2", "kappa = sqrt(z)")



## Table S.8, Linear dose-response curve
design.kappa.tb




############ Squared Root DGP


set.seed(100)
population_data = generate_sqrt_effect(I=500000)

## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")





## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = sqrt(z)

wilcoxon.design.sqrt      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.sqrt)
dose.weighted.design.sqrt = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.sqrt)
U878.design.sqrt          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.sqrt)
dose.U878.design.sqrt     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.sqrt)
design.sqrt = c(wilcoxon.design.sqrt, dose.weighted.design.sqrt, U878.design.sqrt, dose.U878.design.sqrt)
names(design.sqrt) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")



## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.log, design.square, design.sqrt)
rownames(design.kappa.tb) = c("kappa=identity","kappa=log","kappa=z^2", "kappa = sqrt(z)")



## Table S.8, Square root dose-response curve
design.kappa.tb






############ Flat DGP


set.seed(100)
population_data = generate_flat_effect(I=500000)

## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")





## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = sqrt(z)

wilcoxon.design.sqrt      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.sqrt)
dose.weighted.design.sqrt = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.sqrt)
U878.design.sqrt          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.sqrt)
dose.U878.design.sqrt     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.sqrt)
design.sqrt = c(wilcoxon.design.sqrt, dose.weighted.design.sqrt, U878.design.sqrt, dose.U878.design.sqrt)
names(design.sqrt) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")



## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.log, design.square, design.sqrt)
rownames(design.kappa.tb) = c("kappa=identity","kappa=log","kappa=z^2", "kappa = sqrt(z)")



## Table S.8, Flat dose-response curve
design.kappa.tb







############ Log DGP


set.seed(100)
population_data = generate_log_effect(I=500000)


## kappa = identity
wilcoxon.design.id      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.identity)
dose.weighted.design.id = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.identity)
U878.design.id          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.identity)
dose.U878.design.id     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.identity)
design.identity = c(wilcoxon.design.id, dose.weighted.design.id, U878.design.id, dose.U878.design.id)
names(design.identity) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")





## kappa = log(z)
wilcoxon.design.log      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.log)
dose.weighted.design.log = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.log)
U878.design.log          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.log)
dose.U878.design.log     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.log)
design.log = c(wilcoxon.design.log, dose.weighted.design.log, U878.design.log, dose.U878.design.log)
names(design.log) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = z^2
wilcoxon.design.sq      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.square)
dose.weighted.design.sq = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.square)
U878.design.sq          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.square)
dose.U878.design.sq     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.square)
design.square = c(wilcoxon.design.sq, dose.weighted.design.sq, U878.design.sq, dose.U878.design.sq)
names(design.square) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


## kappa = sqrt(z)

wilcoxon.design.sqrt      = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon",        kappa=kappa.sqrt)
dose.weighted.design.sqrt = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted",   kappa=kappa.sqrt)
U878.design.sqrt          = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",               m=8,m1=7,m2=8, kappa=kappa.sqrt)
dose.U878.design.sqrt     = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U", m=8,m1=7,m2=8, kappa=kappa.sqrt)
design.sqrt = c(wilcoxon.design.sqrt, dose.weighted.design.sqrt, U878.design.sqrt, dose.U878.design.sqrt)
names(design.sqrt) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")



## Stack into a design-sensitivity table, rows = kappa, cols = method
design.kappa.tb = rbind(design.identity, design.log, design.square, design.sqrt)
rownames(design.kappa.tb) = c("kappa=identity","kappa=log","kappa=z^2", "kappa = sqrt(z)")

## Table S.8, Log dose-response curve
design.kappa.tb












