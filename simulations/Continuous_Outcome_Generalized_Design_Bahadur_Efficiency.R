#####################################################################################################################################
#                                                                                                                                   #
#      This file contains the simulations for                                                                                       #
#      "Towards Robust Matched Observational Studies with General Treatment Types: Consistency, Efficiency, and Adaptivity"         #
#      This file reproduces Table 2, Table 3, Table 4, and Fig 2 in the main text, and Table S.1 and Table S.2 with binary outcome  #
#      in the supplementary material, which takes the \kappa(z) = z in the text                                                     #
#      design stage, and fix \phi(z) = z for the sensitivity analysis. See Section C.3 for discussion on applying different         #
#      \kappa(z), the simulations for different \kappa(z) choice can be seen in                                                     #
#       Dose_Response_Curve_Transformation_Verification_Comparison.R and Bahadur_Efficiency_Simulations_with_Transformation.R       #                       
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################


source("matched_sensitivity_analysis_functions.R")






####################################################################################################
#########             Here starts our experiment                                            ########
####################################################################################################



###################################################################################
#            The plot of various dose-response curves in the main text, Figure 2, #
#             with their design sensitivities                                     #
###################################################################################


###################################################################################
library(ggplot2)

pop1      <- generate_linear_effect(I = 500000, linear_coefficient = 1)
pop_kink  <- generate_kink_effect(I = 500000, linear_coefficient = 1.5, kink = 0.8)
## The users may use this function to specify different polynomial dose-response
## curves, in our manuscript, we specify the square dose-response curve.
pop_poly  <- generate_polynomial_effect(I = 500000, power = 2)
pop_flat  <- generate_flat_effect(I = 500000)
pop_sqrt  <- generate_sqrt_effect(I = 500000)
pop_log   <- generate_log_effect(I = 500000)

make_df <- function(population_data, label) {
  data.frame(
    dose     = c(population_data$dose[,1], population_data$dose[,2]),
    response = c(population_data$response[,1], population_data$response[,2]),
    effect   = label
  )
}

all_data <- rbind(
  make_df(pop_poly, "Square"),
  make_df(pop_kink, "Kink"),
  make_df(pop1,     "Linear"),
  make_df(pop_sqrt, "Square Root"),
  make_df(pop_flat, "Flat"),
  make_df(pop_log,  "Log")
)

# -------------------------
# Step 3: Build ggplot
# -------------------------
# Keep original short names in factor
all_data$effect <- factor(all_data$effect,
                          levels = c("Square", "Kink", "Linear", "Square Root", "Flat", "Log")
)

# Legend labels with just the numbers
labels_numbers_only <- c(
  "2.32, 2.44, 3.32, 4.16",
  "2.26, 2.43, 3.17, 4.28", 
  "2.68, 2.68, 4.16, 4.43",
  "2.45, 2.38, 3.67, 3.50",
  "2.60, 2.45, 4.02, 3.37",
  "2.83, 2.63, 4.65, 3.84"
)

dose.response.curve.p <- ggplot(all_data, aes(x = dose, y = response, color = effect, linetype = effect)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, size = 1.2) +
  
  # More contrasting colors
  scale_color_manual(
    values = c(
      "Square"   = "#E31A1C",  # Bright red
      "Kink"         = "#1F78B4",  # Blue  
      "Linear"       = "#33A02C",  # Green
      "Square Root" = "#FF7F00",  # Orange
      "Flat"         = "#6A3D9A",  # Purple
      "Log"          = "#A6761D"   # Brown
    ),
    labels = labels_numbers_only
  ) +
  
  # Group line types: solid for Square/Kink/Linear, dotdash for others
  scale_linetype_manual(
    values = c(
      "Square"   = "solid",
      "Kink"         = "solid", 
      "Linear"       = "solid",
      "Square Root" = "longdash",
      "Flat"         = "longdash",
      "Log"          = "longdash"
    ), 
    labels = labels_numbers_only
  ) +
  
  scale_x_continuous(limits = c(0, 2.0), expand = c(0, 0.05)) +
  labs(x = "", y = "", 
       color = "", 
       linetype = "") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, family="Times New Roman"),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.key.spacing.y = unit(0.6, "cm"),
    legend.box.margin = margin(0, 0, 0, -10),
    legend.margin = margin(0, 0, 0, -5),
    plot.margin = margin(5, 5, 5, 5),
    panel.grid.major = element_line(color = "gray85", size = 0.5),
    panel.grid.minor = element_line(color = "gray85", size = 0.25),
    axis.text = element_text(size = 18, family="Times New Roman"),
    axis.text.y = element_text(margin = margin(r = 2)),
    axis.text.x = element_text(margin = margin(t = 2)),
    aspect.ratio = 0.6
  ) +
  coord_cartesian(clip = "on")  # Keep this to restrict curve drawing

# Save as PNG
ggsave("DoseResponseCurve.png", plot = dose.response.curve.p,
       width = 9, height = 6.5, dpi = 600)



#################################################################################
#                                                                               # 
#       Main text, Table 2, Design Sensitivity and simulated power with         #
#.          Square Treatment Effect                                             #
#                                                                               #
#################################################################################



set.seed(100)
population_data = generate_polynomial_effect(I=500000, power=2)

wilcoxon.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon")
dose.weighted.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted")
U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",m=8,m1=7,m2=8)
dose.U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U",m=8,m1=7,m2=8)
design.power.2.effect = c(wilcoxon.design,dose.weighted.design,U878.design, dose.U878.design)
names(design.power.2.effect) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")

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

wilcoxon.simulated.power      = simulated_power_polynomial_2(method="wilcoxon",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.weighted.simulated.power = simulated_power_polynomial_2(method="dose.weighted",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
U878.simulated.power = simulated_power_polynomial_2(method="U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.U878.simulated.power = simulated_power_polynomial_2(method="dose.weighted.U", m=8,m1=7,m2=8, Gamma_seq=c(1.5,2,2.5,3,3.5,4))

polynomial_2_simulated_power = rbind(wilcoxon.simulated.power$Gamma_seq, wilcoxon.simulated.power$simulated.power,dose.weighted.simulated.power$simulated.power, 
                                     U878.simulated.power$simulated.power,dose.U878.simulated.power$simulated.power)
rownames(polynomial_2_simulated_power) = c("Gamma","Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


design.power.2.effect

round(t(polynomial_2_simulated_power),2)










#################################################################################
#                                                                               # 
#             Main text, Table 2, Design Sensitivity and simulated power        # 
#                             with Kink Treatment Effect                        #
#                                                                               #
#################################################################################


set.seed(100)
population_data = generate_kink_effect(I=500000, linear_coefficient = 1.5, kink = 0.8)

wilcoxon.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon")


dose.weighted.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted")


U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",m=8,m1=7,m2=8)

dose.U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U",m=8,m1=7,m2=8)


design.kink.effect = c(wilcoxon.design,dose.weighted.design,U878.design, dose.U878.design)
names(design.kink.effect) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


simulated_power_kink = function(method=c("wilcoxon","dose.weighted","U", "dose.weighted.U"),rank.order = 2, m=2,m1=2,m2=2,Gamma_seq, alpha=0.01,linear_coefficient= 1.5, kink  = 0.8, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_kink_effect(I=5000,seed=seed, linear_coefficient = linear_coefficient, kink = kink)
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

wilcoxon.simulated.power      = simulated_power_kink(method="wilcoxon",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.weighted.simulated.power = simulated_power_kink(method="dose.weighted",Gamma_seq=c(1.5,2,2.5,3,3.5,4))

U878.simulated.power = simulated_power_kink(method="U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))


dose.U878.simulated.power = simulated_power_kink(method="dose.weighted.U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))

kink_simulated_power = rbind(wilcoxon.simulated.power$Gamma_seq, wilcoxon.simulated.power$simulated.power, 
                             dose.weighted.simulated.power$simulated.power, 
                             U878.simulated.power$simulated.power,dose.U878.simulated.power$simulated.power)
rownames(kink_simulated_power) = c("Gamma","Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


design.kink.effect

round(t(kink_simulated_power),2)






#################################################################################
#           Main text, Table, Design Sensitivity and simulated powers           #
#                              with Linear Treatment Effect                     #
#                                                                               #
#################################################################################


set.seed(100)
population_data = generate_linear_effect(I=500000, linear_coefficient =1)

wilcoxon.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon")


dose.weighted.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted")


U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",m=8,m1=7,m2=8)

U878.dose.design = design_sensitivity_fun(x=population_data$dose, y=population_data$response, method="dose.weighted.U", m=8,m1=7,m2=8)
# Report the design sensitivities

design.linear.1.effect = c(wilcoxon.design,dose.weighted.design,U878.design,U878.dose.design)
names(design.linear.1.effect) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


simulated_power_linear_1 = function(method=c("wilcoxon","dose.weighted","U","dose.weighted.U"),rank.order = 2, m=2,m1=2,m2=2,Gamma_seq, alpha=0.01,linear_coefficient = 1, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_linear_effect(I=5000,seed=seed, linear_coefficient = linear_coefficient)
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

wilcoxon.simulated.power      = simulated_power_linear_1(method="wilcoxon",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.weighted.simulated.power = simulated_power_linear_1(method="dose.weighted",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
U878.simulated.power = simulated_power_linear_1(method="U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.U878.simulated.power = simulated_power_linear_1(method = "dose.weighted.U",m=8,m1=7,m2=8, Gamma_seq = c(1.5,2,2.5,3,3.5,4))

linear_1_simulated_power = rbind(wilcoxon.simulated.power$Gamma_seq, wilcoxon.simulated.power$simulated.power, 
                                 dose.weighted.simulated.power$simulated.power, 
                                 U878.simulated.power$simulated.power, dose.U878.simulated.power$simulated.power)
rownames(linear_1_simulated_power) = c("Gamma","Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


design.linear.1.effect

round(t(linear_1_simulated_power),2)






#################################################################################
#                                                                               # 
#             Main text, Table 2, Design Sensitivity and simulated powers       #
#                          with Squared Root Treatment Effect                   #
#                                                                               #
#################################################################################


set.seed(100)
population_data = generate_sqrt_effect(I=500000)

wilcoxon.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon")


dose.weighted.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted")


U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",m=8,m1=7,m2=8)

U878.dose.design = design_sensitivity_fun(x=population_data$dose, y=population_data$response, method="dose.weighted.U", m=8,m1=7,m2=8)
# Report the design sensitivities

design.sqrt.effect = c(wilcoxon.design,dose.weighted.design,U878.design,U878.dose.design)
names(design.sqrt.effect) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")
round(design.sqrt.effect,2)



simulated_power_sqrt = function(method=c("wilcoxon","dose.weighted","U","dose.weighted.U"),rank.order = 2, m=2,m1=2,m2=2,Gamma_seq, alpha=0.01, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_sqrt_effect(I=5000,seed=seed)
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

wilcoxon.simulated.power      = simulated_power_sqrt(method="wilcoxon",Gamma_seq=c(1.5,2,2.5,3,3.5,4.0))
dose.weighted.simulated.power = simulated_power_sqrt(method="dose.weighted",Gamma_seq=c(1.5,2,2.5,3,3.5,4))

U878.simulated.power = simulated_power_sqrt(method="U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.U878.simulated.power = simulated_power_sqrt(method="dose.weighted.U", m=8, m1=7,m2=8,Gamma_seq = c(1.5,2,2.5,3,3.5,4))
sqrt_simulated_power = rbind(wilcoxon.simulated.power$Gamma_seq, wilcoxon.simulated.power$simulated.power, 
                             dose.weighted.simulated.power$simulated.power,
                             U878.simulated.power$simulated.power,dose.U878.simulated.power$simulated.power)
rownames(sqrt_simulated_power) = c("Gamma","Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


design.sqrt.effect

round(t(sqrt_simulated_power),2)




#################################################################################
#                                                                               # 
#             Main text, Table 2, Design Sensitivity and simulated powers       # 
#                       with Flat Treatment Effect                              #
#                                                                               #
#################################################################################


set.seed(100)
population_data = generate_flat_effect(I=500000)
wilcoxon.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon")
dose.weighted.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted")
U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",m=8,m1=7,m2=8)
dose.U878.design =design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U",m=8,m1=7,m2=8)
design.flat.effect = c(wilcoxon.design,dose.weighted.design,U878.design,dose.U878.design)
names(design.flat.effect) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


simulated_power_flat = function(method=c("wilcoxon","dose.weighted","U","dose.weighted.U"),rank.order = 2, m=2,m1=2,m2=2,Gamma_seq, alpha=0.01, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_flat_effect(I=5000,seed=seed)
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

wilcoxon.simulated.power      = simulated_power_flat(method="wilcoxon",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.weighted.simulated.power = simulated_power_flat(method="dose.weighted",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
U878.simulated.power = simulated_power_flat(method="U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.U878.simulated.power = simulated_power_flat(method="dose.weighted.U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))


flat_simulated_power = rbind(wilcoxon.simulated.power$Gamma_seq, wilcoxon.simulated.power$simulated.power, 
                             dose.weighted.simulated.power$simulated.power, 
                             U878.simulated.power$simulated.power, dose.U878.simulated.power$simulated.power)
rownames(flat_simulated_power) = c("Gamma","Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


design.flat.effect

round(t(flat_simulated_power),2)





#################################################################################
#                                                                               # 
#             Main text, Table 2, Design Sensitivity and simulated powers 
#                              with Log Treatment Effect                        #
#                                                                               #
#################################################################################


set.seed(100)
population_data = generate_log_effect(I=500000)

wilcoxon.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="wilcoxon")
dose.weighted.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted")
U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="U",m=8,m1=7,m2=8)
dose.U878.design = design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose.weighted.U",m=8,m1=7,m2=8)

design.log.effect = c(wilcoxon.design,dose.weighted.design,U878.design,dose.U878.design)
names(design.log.effect) = c("Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


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

wilcoxon.simulated.power      = simulated_power_log(method="wilcoxon",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.weighted.simulated.power = simulated_power_log(method="dose.weighted",Gamma_seq=c(1.5,2,2.5,3,3.5,4))
U878.simulated.power = simulated_power_log(method="U",m=8,m1=7,m2=8,Gamma_seq=c(1.5,2,2.5,3,3.5,4))
dose.U878.simulated.power = simulated_power_log(method="dose.weighted.U", m=8, m1=7,m2=8,Gamma_seq = c(1.5,2,2.5,3,3.5,4))
log_simulated_power = rbind(wilcoxon.simulated.power$Gamma_seq, wilcoxon.simulated.power$simulated.power, 
                            dose.weighted.simulated.power$simulated.power,
                            U878.simulated.power$simulated.power,dose.U878.simulated.power$simulated.power)
rownames(log_simulated_power) = c("Gamma","Wilcoxon","dose-wt","(8,7,8)","dose.(8,7,8)")


design.log.effect

round(t(log_simulated_power),2)








########################################################################################
####### Here starts the experiment of Bahadur-Rosenbaum relative efficiency ############
########################################################################################



## this is a helper function to tell the column giving the row-wise maximum
max_col_per_row <- function(mat) {
  # Check if input is a matrix
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }
  
  # Find the column index of the max value in each row
  max_cols <- apply(mat, 1, which.max)
  
  return(max_cols)
}





########################################################################################
#                                                                                      #
#     Main text, Table 3                                                               #
#    Bahadur-Rosenbaum relative efficiency, Square treatment effect                    #
#                                                                                      #
########################################################################################

set.seed(100)
population_data = generate_polynomial_effect(I=500000,power = 2)

wilcoxon_slope_vec = rep(NA,4)
dose.weighted_slope_vec = rep(NA,4)
U878_slope_vec = rep(NA,4)
dose.U878_slope_vec = rep(NA,4)
Gamma_ave_vec = c(1.0,1.1,1.3,1.5)
for(index in 1:4){
  wilcoxon_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="wilcoxon")
  dose.weighted_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted")
  U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="U",m=8,m1=7,m2=8)
  dose.U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted.U",m=8,m1=7,m2=8)
  
}




polynomial.effect.Bahadur.relative.efficiency = cbind(wilcoxon_slope_vec/dose.weighted_slope_vec,U878_slope_vec/dose.U878_slope_vec)


sample_size_polynomial_fun = function(alpha,method = c("wilcoxon","dose.weighted","U","dose.weighted.U"),power,Gamma_ave, upper.bound=1000,rank.order=2,m1=2,m2=2,m=2, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      
      sample_data = generate_polynomial_effect(I=I,seed=mc, power = 2)
      
      sample.data_x =sample_data$dose
      sample.data_y =sample_data$response
      
      # Compute the p-value
      p.value.vec[mc] = sendose(x = sample.data_x,y = sample.data_y,method = method,  
                                alternative = "greater than",  Gamma_ave = Gamma_ave,rank.order = rank.order,m=m,m1=m1,m2=m2, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}



dose.to.wilcoxon.sample.size.ratio.0.01 = rep(NA,4)
dose.U878.to.U878.sample.size.ratio.0.01 = rep(NA,4)
wilcoxon.sample.size.0.01 = rep(NA,4)
U878.sample.size.0.01 =rep(NA,4)


set.seed(100)

for(index in 1:4){
  wilcoxon.size = sample_size_polynomial_fun(alpha=0.01,method="wilcoxon", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  wilcoxon.sample.size.0.01[index] = wilcoxon.size
  dose.weighted.rank.size = sample_size_polynomial_fun(alpha=0.01,method="dose.weighted", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  U878.size = sample_size_polynomial_fun(alpha=0.01,method="U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  dose.U878.size = sample_size_polynomial_fun(alpha=0.01,method="dose.weighted.U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  U878.sample.size.0.01[index] = U878.size
  dose.to.wilcoxon.sample.size.ratio.0.01[index] = dose.weighted.rank.size/wilcoxon.size
  dose.U878.to.U878.sample.size.ratio.0.01[index] = dose.U878.size/U878.size
  
  
  
}
alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01,dose.U878.to.U878.sample.size.ratio.0.01)

wilcoxon.sample.size.0.01
U878.sample.size.0.01

## The minimal sample size ratio from bisection method 
round(alpha.0.01.size.ratio,2)

## The Bahadur relatie efficiency by taking the ratio of the two Bahadur slopes 
round(polynomial.effect.Bahadur.relative.efficiency,2)



## examine the Bahadur slope altogether 
Bahadur.mat = cbind(wilcoxon_slope_vec, dose.weighted_slope_vec, U878_slope_vec, dose.U878_slope_vec)
max_col_per_row(Bahadur.mat)









########################################################################################
#                                                                                      #
#          Main text, Table 3
#       Bahadur-Rosenbaum relative efficiency, Kink treatment effect                   #
#                                                                                      #
########################################################################################


set.seed(100)
population_data = generate_kink_effect(I=500000, linear_coefficient = 1.5, kink=0.8)

wilcoxon_slope_vec = rep(NA,4)
dose.weighted_slope_vec = rep(NA,4)
U878_slope_vec = rep(NA,4)
dose.U878_slope_vec = rep(NA,4)
Gamma_ave_vec = c(1.0,1.1,1.3,1.5)
for(index in 1:4){
  wilcoxon_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="wilcoxon")
  dose.weighted_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted")
  U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="U",m=8,m1=7,m2=8)
  dose.U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted.U",m=8,m1=7,m2=8)
  
}

kink.effect.Bahadur.relative.efficiency = cbind(wilcoxon_slope_vec/dose.weighted_slope_vec,U878_slope_vec/dose.U878_slope_vec)


sample_size_kink_fun = function(alpha,method = c("wilcoxon","dose.weighted","U","dose.weighted.U"),power,Gamma_ave, upper.bound=1500,rank.order=2,m1=2,m2=2,m=2, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      
      sample_data = generate_kink_effect(I=I,seed=mc, linear_coefficient = 1.5, kink = 0.8)
      
      sample.data_x =sample_data$dose
      sample.data_y =sample_data$response
      
      # Compute the p-value
      p.value.vec[mc] = sendose(x = sample.data_x,y = sample.data_y,method = method,  
                                alternative = "greater than",  Gamma_ave = Gamma_ave,rank.order = rank.order,m=m,m1=m1,m2=m2, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}


dose.to.wilcoxon.sample.size.ratio.0.01 = rep(NA,4)
dose.U878.to.U878.sample.size.ratio.0.01 = rep(NA,4)
wilcoxon.sample.size.0.01 = rep(NA,4)
U878.sample.size.0.01 =rep(NA,4)


set.seed(100)

for(index in 1:4){
  wilcoxon.size = sample_size_kink_fun(alpha=0.01,method="wilcoxon", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  wilcoxon.sample.size.0.01[index] = wilcoxon.size
  dose.weighted.rank.size = sample_size_kink_fun(alpha=0.01,method="dose.weighted", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  
  U878.size = sample_size_kink_fun(alpha=0.01,method="U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  dose.U878.size = sample_size_kink_fun(alpha=0.01,method="dose.weighted.U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  U878.sample.size.0.01[index] = U878.size
  
  dose.to.wilcoxon.sample.size.ratio.0.01[index] = dose.weighted.rank.size/wilcoxon.size
  dose.U878.to.U878.sample.size.ratio.0.01[index] = dose.U878.size/U878.size
  
  # alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01 ,dose.U878.to.U878.sample.size.ratio.0.01)
  
  
  
  # print(round(alpha.0.01.size.ratio,2))
  # print(round(kink.effect.Bahadur.relative.efficiency,2))
  
  
  
}

wilcoxon.sample.size.0.01

U878.sample.size.0.01

alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01 ,dose.U878.to.U878.sample.size.ratio.0.01)

round(alpha.0.01.size.ratio,2)

round(kink.effect.Bahadur.relative.efficiency,2)

## examine the Bahadur slope altogether 
Bahadur.mat = cbind(wilcoxon_slope_vec, dose.weighted_slope_vec, U878_slope_vec, dose.U878_slope_vec)
max_col_per_row(Bahadur.mat)









########################################################################################
#                                                                                      #
#         Main text, Table 3
#        Bahadur-Rosenbaum relative efficiency, Linear treatment effect                #
#                                                                                      #
########################################################################################


set.seed(100)
population_data = generate_linear_effect(I=500000, linear_coefficient =1)

wilcoxon_slope_vec = rep(NA,4)
dose.weighted_slope_vec = rep(NA,4)
U878_slope_vec = rep(NA,4)
dose.U878_slope_vec = rep(NA,4)
Gamma_ave_vec = c(1.0,1.1,1.3,1.5)
for(index in 1:4){
  wilcoxon_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="wilcoxon")
  dose.weighted_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted")
  U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="U",m=8,m1=7,m2=8)
  dose.U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted.U",m=8,m1=7,m2=8)
  
}


linear.1.effect.Bahadur.relative.efficiency = cbind(wilcoxon_slope_vec/dose.weighted_slope_vec,U878_slope_vec/dose.U878_slope_vec)

sample_size_linear_fun = function(alpha,method = c("wilcoxon","dose.weighted","U","dose.weighted.U"),power,Gamma_ave, upper.bound=1000,rank.order=2,m1=2,m2=2,m=2, linear_coefficient, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      
      sample_data = generate_linear_effect(I=I,seed=mc, linear_coefficient = linear_coefficient)
      
      sample.data_x =sample_data$dose
      sample.data_y =sample_data$response
      
      # Compute the p-value
      p.value.vec[mc] = sendose(x = sample.data_x,y = sample.data_y,method = method,  
                                alternative = "greater than",  Gamma_ave = Gamma_ave,rank.order = rank.order,m=m,m1=m1,m2=m2, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}


dose.to.wilcoxon.sample.size.ratio.0.01 = rep(NA,4)
dose.U878.to.U878.sample.size.ratio.0.01 = rep(NA,4)
wilcoxon.sample.size.0.01 = rep(NA,4)
U878.sample.size.0.01 =rep(NA,4)


set.seed(100)

for(index in 1:4){
  wilcoxon.size = sample_size_linear_fun(alpha=0.01,method="wilcoxon", power = 0.95, Gamma_ave = Gamma_ave_vec[index], linear_coefficient = 1)
  wilcoxon.sample.size.0.01[index] = wilcoxon.size
  dose.weighted.rank.size = sample_size_linear_fun(alpha=0.01,method="dose.weighted", power = 0.95, Gamma_ave = Gamma_ave_vec[index], linear_coefficient = 1)
  U878.size = sample_size_linear_fun(alpha=0.01,method="U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8, linear_coefficient = 1)
  dose.U878.size = sample_size_linear_fun(alpha=0.01,method="dose.weighted.U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8, linear_coefficient = 1)
  U878.sample.size.0.01[index] = U878.size
  
  dose.to.wilcoxon.sample.size.ratio.0.01[index] = dose.weighted.rank.size/wilcoxon.size
  dose.U878.to.U878.sample.size.ratio.0.01[index] = dose.U878.size/U878.size
  alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01, dose.U878.to.U878.sample.size.ratio.0.01)
  print(round(alpha.0.01.size.ratio,2))
  print(round(linear.1.effect.Bahadur.relative.efficiency,2))
  
  
  
  
}
alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01, dose.U878.to.U878.sample.size.ratio.0.01)

wilcoxon.sample.size.0.01
U878.sample.size.0.01
round(alpha.0.01.size.ratio,2)

round(linear.1.effect.Bahadur.relative.efficiency,2)


## examine the Bahadur slope altogether 
Bahadur.mat = cbind(wilcoxon_slope_vec, dose.weighted_slope_vec, U878_slope_vec, dose.U878_slope_vec)
max_col_per_row(Bahadur.mat)






########################################################################################
#                                                                                      #
# Main text, Table 3, 
# Bahadur-Rosenbaum relative efficiency, Squared Root treatment effect                 #
#                                                                                      #
########################################################################################


set.seed(100)
population_data = generate_sqrt_effect(I=500000)


wilcoxon_slope_vec = rep(NA,4)
dose.weighted_slope_vec = rep(NA,4)
U878_slope_vec = rep(NA,4)
dose.U878_slope_vec = rep(NA,4)
Gamma_ave_vec = c(1.0,1.1,1.3,1.5)
for(index in 1:4){
  wilcoxon_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="wilcoxon")
  dose.weighted_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted")
  U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="U",m=8,m1=7,m2=8)
  dose.U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted.U",m=8,m1=7,m2=8)
  
}


sqrt.effect.Bahadur.relative.efficiency = cbind(wilcoxon_slope_vec/dose.weighted_slope_vec,U878_slope_vec/dose.U878_slope_vec)

sample_size_sqrt_fun = function(alpha,method = c("wilcoxon","dose.weighted","U","dose.weighted.U"),power,Gamma_ave, upper.bound=1000,rank.order=2,m1=2,m2=2,m=2, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      
      sample_data = generate_sqrt_effect(I=I,seed=mc)
      
      sample.data_x =sample_data$dose
      sample.data_y =sample_data$response
      
      # Compute the p-value
      p.value.vec[mc] = sendose(x = sample.data_x,y = sample.data_y,method = method,  
                                alternative = "greater than",  Gamma_ave = Gamma_ave,rank.order = rank.order,m=m,m1=m1,m2=m2, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}



dose.to.wilcoxon.sample.size.ratio.0.01 = rep(NA,4)
dose.U878.to.U878.sample.size.ratio.0.01 = rep(NA,4)
wilcoxon.sample.size.0.01 = rep(NA,4)
U878.sample.size.0.01 =rep(NA,4)


set.seed(100)


for(index in 1:4){
  wilcoxon.size = sample_size_sqrt_fun(alpha=0.01,method="wilcoxon", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  wilcoxon.sample.size.0.01[index] = wilcoxon.size
  dose.weighted.rank.size = sample_size_sqrt_fun(alpha=0.01,method="dose.weighted", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  U878.size = sample_size_sqrt_fun(alpha=0.01,method="U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  dose.U878.size = sample_size_sqrt_fun(alpha=0.01,method="dose.weighted.U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  U878.sample.size.0.01[index] = U878.size
  
  dose.to.wilcoxon.sample.size.ratio.0.01[index] = dose.weighted.rank.size/wilcoxon.size
  dose.U878.to.U878.sample.size.ratio.0.01[index] = dose.U878.size/U878.size
  alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01,
                                dose.U878.to.U878.sample.size.ratio.0.01)
  
  # print(round(alpha.0.01.size.ratio,2))
  # print(round(sqrt.effect.Bahadur.relative.efficiency,2))
  
  
  
  
  
  
}
alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01,dose.U878.to.U878.sample.size.ratio.0.01)

U878.sample.size.0.01
wilcoxon.sample.size.0.01

round(alpha.0.01.size.ratio,2)

round(sqrt.effect.Bahadur.relative.efficiency,2)




## examine the Bahadur slope altogether 
Bahadur.mat = cbind(wilcoxon_slope_vec, dose.weighted_slope_vec, U878_slope_vec, dose.U878_slope_vec)
max_col_per_row(Bahadur.mat)








########################################################################################
#                                                                                      #
#   Main text, Table 3
#   Bahadur-Rosenbaum  relative efficiency, Flat treatment effect                      #
#                                                                                      #
########################################################################################


set.seed(100)
population_data = generate_flat_effect(I=500000)

wilcoxon_slope_vec = rep(NA,4)
dose.weighted_slope_vec = rep(NA,4)
U878_slope_vec = rep(NA,4)
dose.U878_slope_vec = rep(NA,4)
Gamma_ave_vec = c(1.0,1.1,1.3,1.5)
for(index in 1:4){
  wilcoxon_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="wilcoxon")
  dose.weighted_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted")
  U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="U",m=8,m1=7,m2=8)
  dose.U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted.U",m=8,m1=7,m2=8)
  
}


flat.effect.Bahadur.relative.efficiency = cbind(wilcoxon_slope_vec/dose.weighted_slope_vec,U878_slope_vec/dose.U878_slope_vec)



sample_size_flat_fun = function(alpha,method = c("wilcoxon","dose.weighted","U","dose.weighted.U"),power,Gamma_ave, upper.bound=1000,rank.order=2,m1=2,m2=2,m=2, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      
      sample_data = generate_flat_effect(I=I,seed=mc)
      
      sample.data_x =sample_data$dose
      sample.data_y =sample_data$response
      
      # Compute the p-value
      p.value.vec[mc] = sendose(x = sample.data_x,y = sample.data_y,method = method,  
                                alternative = "greater than",  Gamma_ave = Gamma_ave,rank.order = rank.order,m=m,m1=m1,m2=m2, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}


dose.to.wilcoxon.sample.size.ratio.0.01 = rep(NA,4)
dose.U878.to.U878.sample.size.ratio.0.01 = rep(NA,4)
wilcoxon.sample.size.0.01 = rep(NA,4)
U878.sample.size.0.01 =rep(NA,4)


set.seed(100)

for(index in 1:4){
  wilcoxon.size = sample_size_flat_fun(alpha=0.01,method="wilcoxon", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  wilcoxon.sample.size.0.01[index] = wilcoxon.size
  dose.weighted.rank.size = sample_size_flat_fun(alpha=0.01,method="dose.weighted", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  U878.size = sample_size_flat_fun(alpha=0.01,method="U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  dose.U878.size = sample_size_flat_fun(alpha=0.01,method="dose.weighted.U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  U878.sample.size.0.01[index] = U878.size
  
  dose.to.wilcoxon.sample.size.ratio.0.01[index] = dose.weighted.rank.size/wilcoxon.size
  dose.U878.to.U878.sample.size.ratio.0.01[index] = dose.U878.size/U878.size
  
  
  
  
}
alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01,dose.U878.to.U878.sample.size.ratio.0.01)

wilcoxon.sample.size.0.01
U878.sample.size.0.01

round(alpha.0.01.size.ratio,2)

round(flat.effect.Bahadur.relative.efficiency,2)





## examine the Bahadur slope altogether 
Bahadur.mat = cbind(wilcoxon_slope_vec, dose.weighted_slope_vec, U878_slope_vec, dose.U878_slope_vec)
max_col_per_row(Bahadur.mat)




########################################################################################
#                                                                                      #
#          Main text, Table 3 
#       Bahadur-Rosebaum relative efficiency, Log treatment effect                     #
#                                                                                      #
########################################################################################


set.seed(100)
population_data = generate_log_effect(I=500000)

wilcoxon_slope_vec = rep(NA,4)
dose.weighted_slope_vec = rep(NA,4)
U878_slope_vec = rep(NA,4)
dose.U878_slope_vec = rep(NA,4)
Gamma_ave_vec = c(1.0,1.1,1.3,1.5)
for(index in 1:4){
  wilcoxon_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="wilcoxon")
  dose.weighted_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted")
  U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="U",m=8,m1=7,m2=8)
  dose.U878_slope_vec[index] = Bahadur_slope_fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose.weighted.U",m=8,m1=7,m2=8)
  
}

log.effect.Bahadur.relative.efficiency =  cbind(wilcoxon_slope_vec/dose.weighted_slope_vec,U878_slope_vec/dose.U878_slope_vec)


sample_size_log_fun = function(alpha,method = c("wilcoxon","dose.weighted","U","dose.weighted.U"),power,Gamma_ave, upper.bound=2000,rank.order=2,m1=2,m2=2,m=2, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      
      sample_data = generate_log_effect(I=I,seed=mc)
      
      sample.data_x =sample_data$dose
      sample.data_y =sample_data$response
      
      # Compute the p-value
      p.value.vec[mc] = sendose(x = sample.data_x,y = sample.data_y,method = method,  
                                alternative = "greater than",  Gamma_ave = Gamma_ave,rank.order = rank.order,m=m,m1=m1,m2=m2, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}


dose.to.wilcoxon.sample.size.ratio.0.01 = rep(NA,4)
dose.U878.to.U878.sample.size.ratio.0.01 = rep(NA,4)
wilcoxon.sample.size.0.01 = rep(NA,4)
U878.sample.size.0.01 =rep(NA,4)


set.seed(100)


for(index in 1:4){
  wilcoxon.size = sample_size_log_fun(alpha=0.01,method="wilcoxon", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  wilcoxon.sample.size.0.01[index] = wilcoxon.size
  dose.weighted.rank.size = sample_size_log_fun(alpha=0.01,method="dose.weighted", power = 0.95, Gamma_ave = Gamma_ave_vec[index])
  U878.size = sample_size_log_fun(alpha=0.01,method="U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  dose.U878.size = sample_size_log_fun(alpha=0.01,method="dose.weighted.U",power=0.95,Gamma_ave = Gamma_ave_vec[index],m=8,m1=7,m2=8)
  U878.sample.size.0.01[index] = U878.size
  dose.to.wilcoxon.sample.size.ratio.0.01[index] = dose.weighted.rank.size/wilcoxon.size
  dose.U878.to.U878.sample.size.ratio.0.01[index] = dose.U878.size/U878.size
  alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01,dose.U878.to.U878.sample.size.ratio.0.01)
  print(round(alpha.0.01.size.ratio,2))
  
  print(round(log.effect.Bahadur.relative.efficiency,2))
  
  
}
alpha.0.01.size.ratio = cbind(dose.to.wilcoxon.sample.size.ratio.0.01,dose.U878.to.U878.sample.size.ratio.0.01)

wilcoxon.sample.size.0.01
U878.sample.size.0.01

round(alpha.0.01.size.ratio,2)

round(log.effect.Bahadur.relative.efficiency,2)



## examine the Bahadur slope altogether 
Bahadur.mat = cbind(wilcoxon_slope_vec, dose.weighted_slope_vec, U878_slope_vec, dose.U878_slope_vec)
max_col_per_row(Bahadur.mat)










####################################################################################################################################
########################## Here starts the simulations for adaptive test under continuous treatments ###############################
####################################################################################################################################



################################################################################
#                                                                              #
#                  Main text, Table 4,                                         #
#  Square effect, Wilcoxon and dose-weighted signed rank, adaptive test        #
# The adaptive test function is called "sendose.adaptive" in                   #
#              matched_sensitivity_analysis_functions.R                        #
#                                                                              #
################################################################################



Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

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
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec





##############################################################################
#
#  Main text, Table 4,
#  Square effect, U and dose-weighted U test, adaptive test
#
##############################################################################



Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

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
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec




######################################################################
#   Main text, Table 4
#   kink model, wilcoxon and dose-weighted signed rank adaptive test
#
######################################################################

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    kink.data = generate_kink_effect(I = 5000,linear_coefficient = 1.5, kink =0.8, seed=run)
    res = sendose.adaptive(
      x = kink.data$dose, 
      y = kink.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec




######################################################################
#   Main text, Table 4 
#   kink model, U test and dose-weighted U test, adaptive test
#
######################################################################

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    kink.data = generate_kink_effect(I = 5000,linear_coefficient = 1.5, kink =0.8, seed=run)
    res = sendose.adaptive(
      x = kink.data$dose, 
      y = kink.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec




###############################################################################
#
#  Main text, Table 4,
#  Linear model, Wilcoxon and dose-weighted signed rank, adaptive test 
#
###############################################################################

Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    linear.data = generate_linear_effect(I = 5000, linear_coefficient = 1, seed=run)
    res = sendose.adaptive(
      x = linear.data$dose, 
      y = linear.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec



################################################################################
#
#   Main text, Table 4, 
#   Linear model, U and dose-weighted U test, adaptive test
#
#
###############################################################################


Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    linear.data = generate_linear_effect(I = 5000, linear_coefficient = 1, seed=run)
    res = sendose.adaptive(
      x = linear.data$dose, 
      y = linear.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec








##############################################################################
#
#  Main text, Table 4,
#  Squared root effect, Wilcoxon and dose-weighted signed rank, adaptive test
#
##############################################################################



Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    sqrt.data = generate_sqrt_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = sqrt.data$dose, 
      y = sqrt.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec






##############################################################################
#
#  Main text, Table 4,
#  Squared root effect, U and dose-weighted U test, adaptive test
#
##############################################################################



Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    sqrt.data = generate_sqrt_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = sqrt.data$dose, 
      y = sqrt.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec




##############################################################################
#
#  Main text, Table 4,
#  Flat effect, Wilcoxon and dose-weighted signed rank, adaptive test
#
##############################################################################



Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    flat.data = generate_flat_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = flat.data$dose, 
      y = flat.data$response, 
      method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec






##############################################################################
#
#  Main text, Table 4,
#  Flat effect, U and dose-weighted U test, adaptive test
#
##############################################################################



Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    flat.data = generate_flat_effect(I = 5000, seed=run)
    res = sendose.adaptive(
      x = flat.data$dose, 
      y = flat.data$response, 
      method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec



##############################################################################
#
#   Main text, Table 4
#   log model, wilcoxon and dose weighted signed rank, adaptive test 
#
##############################################################################


Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

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
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec



################################################################################
#      Main text, Table 4 
#      Log model, U test and dose-weighted U test, adaptive test 
#
################################################################################


Gamma_vec = c(1.5,2.0,2.5,3.0,3.5,4.0)
power.vec = rep(NA, length(Gamma_vec))  

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
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
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
  power.vec[index] = power
}

Gamma_vec
power.vec












####################################################################################################################################
####################################################################################################################################
############ Simulation for Binary outcome Design Sensitivity , Bahadur-Rosenbaum Relative Efficiency and adaptive testing, ######## 
############                        Section C in the Supplementary Material                                                 ########







## Here first writes functions for binary outcome simulations
## mcnemar.Test.stat.fun, mcnemar.q_I.fun, mcnemar.sendose, 
## mcnemar.design_sensitivity_fun, mcnemar.Bhadur.fun, generate_binary_effect
## binary_sample_size_fun, sendose.mcnemar.adaptive


mcnemar.Test.stat.fun <- function(x,y, method = c("mcnemar","dose-mcnemar"), kappa=function(z) z) {
  first.treatment = x[,1]
  second.treatment= x[,2]
  first.response  = y[,1]
  second.response = y[,2]
  treatment.diff = first.treatment - second.treatment
  I <- length(treatment.diff)
  response.diff = first.response  - second.response
  scaled.rank.treatment = rank(abs(kappa(first.treatment) - kappa(second.treatment)))/I
  scaled.rank.response  = rank(abs(response.diff))/I
  
  if (method == "mcnemar") {
    scores  = rep(1,I)
  }
  
  if(method=="dose-mcnemar"){
    scores <- scaled.rank.treatment
  }
  
  sign <- ifelse(treatment.diff*response.diff > 0, 1, 0)
  test.stat <- as.numeric(sign %*% scores)
  return(test.stat)
}


mcnemar.q_I.fun <- function(x,y, method=c("mcnemar","dose-mcnemar"), kappa=function(z) z) {
  first.treatment = x[,1]
  second.treatment= x[,2]
  first.response  = y[,1]
  second.response = y[,2]
  treatment.diff = first.treatment - second.treatment
  I <- length(treatment.diff)
  response.diff = first.response  - second.response
  scaled.rank.treatment = rank(abs(kappa(first.treatment) - kappa(second.treatment)))/I
  scaled.rank.response  = rank(abs(response.diff))/I
  
  if (method == "mcnemar") {
    scores  = rep(1,I)
  }
  
  if(method=="dose-mcnemar"){
    scores <- scaled.rank.treatment
    
  }
  
  
  return(scores)
}




mcnemar.sendose <-function(x, y, method=c("mcnemar", "dose-mcnemar"),alternative=c("less than","greater than","two.sided"),Gamma_ave=1,exact=FALSE,mc.iteration=1000,log.p=FALSE, kappa=function(z) z){
  ### Basic check on input
  # Check on data
  stopifnot(is.matrix(x),is.matrix(y),
            is.numeric(x),(is.numeric(y) || is.logical(y)))
  stopifnot(ncol(x) == 2,ncol(y) == 2,all(!is.na(x)), all(!is.na(y)))
  stopifnot(all(abs(x[,2] - x[,1]) > .Machine$double.eps)) # Checks whether any dose difference are essentially zero.
  
  # Check on hypothesis testing and Gamma_bar
  stopifnot(is.numeric(Gamma_ave) && Gamma_ave> 0)
  alternative = match.arg(alternative,c("greater than","less than","two.sided"))
  method = match.arg(method,c("mcnemar","dose-mcnemar"))
  
  
  
  ## extract the Gamma_i and then p_i
  first.treatment = x[,1]
  second.treatment= x[,2]
  first.response = y[,1]
  second.response = y[,2]
  I = length(first.treatment)
  treatment.diff = first.treatment - second.treatment
  tilde.gamma = Gamma_to_tilde_gamma_fun(x=x, generalized.Gamma = Gamma_ave)
  Gamma_i = exp(tilde.gamma*abs(treatment.diff))
  p_i_plus = Gamma_i/(1+Gamma_i)
  p_i_product = p_i_plus*(1-p_i_plus)
  q_I = mcnemar.q_I.fun(x=x,y=y,method=method, kappa = kappa)
  sense.variance = as.numeric(q_I^2%*%p_i_product)
  ## 
  test.stat = mcnemar.Test.stat.fun(x=x,y=y,method=method, kappa = kappa)
  if(exact==FALSE){
    ## compute the normal p value 
    if(alternative=="greater than"){
      sense.mean = as.numeric(q_I%*%p_i_plus)
      sense.z.score = (test.stat - sense.mean)/sqrt(sense.variance)
      pval = pnorm(q=sense.z.score,lower.tail = FALSE,log.p = log.p)
    }
    else if(alternative=="less than"){
      sense.mean = as.numeric(q_I%*%(1-p_i_plus))
      sense.z.score = (test.stat - sense.mean)/sqrt(sense.variance)
      pval = pnorm(q=sense.z.score,lower.tail = TRUE,log.p = log.p)
    }else{
      ## two tails, need to compute two z scores and two one-tailed p value
      ## and take twice of the mininum one as the final p value 
      upper.tailed.sense.mean = as.numeric(q_I%*%p_i_plus)
      lower.tailed.sense.mean = as.numermic(q_I%*%(1-p_i_plus))
      upper.z.score  = (test.stat - upper.tailed.sense.mean)/sqrt(sense.variance)
      lower.z.score = (test.stat - lower.tailed.sense.mean)/(sqrt(sense.variance))
      
      pval = 2*min(pnorm(q=upper.z.score,lower.tail = F,log.p = log.p),pnorm(q=lower.z.score,lower.tail = T,log.p=log.p))
    }
    
    
  }
  if(exact==TRUE){
    mat = replicate(mc.iteration,rbinom(I,1,p_i_plus))
    pval = switch(alternative,
                  "greater than" = mean((q_I %*% mat) >= test.stat),
                  "less than"= mean((q_I %*% (1-mat)) <= test.stat),
                  "two.sided" = 2*min(mean((q_I %*% mat) >= test.stat),
                                      mean((q_I %*% (1-mat)) <= test.stat)))
    
  }
  return(list(pval=pval, sense.mean=sense.mean, sense.variance = sense.variance))
}


## x and y shall be obtained by simulating a population data 
mcnemar.design_sensitivity_fun = function(x,y,method = c("mcnemar","dose-mcnemar"), kappa=function(z) z){
  ## the number of matched pairs is I
  I = length(x)*0.5
  if(method=="mcnemar"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    
    right.hand.side.expectation = mean(dose.response.product.sign)
    f = function(x){
      mean(exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation
    }
    
    gamma_star = uniroot(f,lower = -10,upper=10,tol=.Machine$double.eps^0.9)$root
    
  }
  
  if(method=="dose-mcnemar"){
    dose.diff = x[,2] - x[,1]
    kappa.dose.diff = kappa(x[,2]) - kappa(x[,1])
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    F_z = ecdf(abs(kappa.dose.diff))
    right.hand.side.expectation = mean(F_z(abs(kappa.dose.diff))*dose.response.product.sign)
    f = function(x){
      mean(F_z(abs(kappa.dose.diff))*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f,lower = -10,upper=20,tol=.Machine$double.eps^0.9)$root
  }
  
  ## and obtain the design sensitivity 
  Gamma_star = mean(exp(gamma_star*abs(x[,2]-x[,1])))
  return(Gamma_star)
  
  
}




mcnemar.Bhadur.fun =function(x,y,Gamma_ave,method=c("mcnemar","dose-mcnemar"), upper_bound = 10, kappa=function(z) z){
  tilde_gamma = Gamma_to_tilde_gamma_fun(x=x, generalized.Gamma = Gamma_ave)
  Gamma_i = exp(tilde_gamma*abs(x[,1]-x[,2]))
  ## from Gamma_i to p_i 
  pi_plus  = Gamma_i/(1+Gamma_i)
  response.diff = y[,1] - y[,2]
  dose.diff = x[,1] - x[,2]
  kappa.dose.diff = kappa(x[,1]) - kappa(x[,2])
  diff.product = response.diff*dose.diff
  F_z = ecdf(abs(kappa.dose.diff))
  I = length(x)/2
  
  
  if(method=="mcnemar"){
    psi_star = rep(1,I)
    
  }
  if(method=="dose-mcnemar"){
    psi_star = F_z(abs(kappa.dose.diff))
  }
  
  
  mu = mean(psi_star*ifelse(diff.product>0,1,0))
  
  ## solve on the log scale 
  log.mu = log(mu)
  
  log.f = function(t){
    log(mean(pi_plus*psi_star*exp(t*psi_star)/(pi_plus*exp(t*psi_star)+(1-pi_plus)))) - log.mu
  }
  
  ## Solve on the log scale
  ## Discuss a proper bound on tilde_t with Advisors later!!
  lower_bound <- .Machine$double.eps^0.9
  tolerance <- .Machine$double.eps
  
  ## Function to adjust the bounds if needed
  while (sign(log.f(lower_bound)) == sign(log.f(upper_bound))) {
    if (sign(log.f(lower_bound)) > 0) {
      lower_bound <- lower_bound-10  ## Decrease lower bound
      warning("Function value at lower bound is positive; decreasing lower bound to ", 
              lower_bound)
    } else if (sign(log.f(upper_bound)) < 0) {
      upper_bound <- upper_bound+10  ## Increase upper bound
      warning("Function value at upper bound is negative; increasing upper bound to ", 
              upper_bound)
    }
  }
  
  ## Solve for the root
  tilde_t <- uniroot(log.f, lower = lower_bound, upper = upper_bound, tol = tolerance)$root
  Bhadur = tilde_t*mu - mean(log(pi_plus*exp(tilde_t*psi_star)+(1-pi_plus)))
  return(Bhadur)
  
}

generate_binary_effect = function(I, linear_coefficient,seed=100){
  ## set.seed 
  set.seed(seed)
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = rnorm(n=1,mean=-1,sd=1)
    dose.diff = runif(n=1,min=1,max = 5)
    high.dose = low.dose + dose.diff
    ## generate the matched pair random effect, delta_i
    matched.random.effect = rnorm(n=1,mean=0,sd=1)
    
    low.dose.error = rnorm(n=1,mean=0,sd=1)
    low.dose.response =  ifelse(linear_coefficient*low.dose + matched.random.effect > low.dose.error,1,0)
    
    high.dose.error = rnorm(n=1,mean=0,sd=1)
    
    high.dose.response =  ifelse(linear_coefficient*high.dose + matched.random.effect > high.dose.error,1,0)
    
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}




binary_sample_size_fun = function(alpha,method = c("mcnemar", "dose-mcnemar"),power,Gamma_ave,linear_coefficient,upper.bound=5000, kappa=function(z) z){
  f = function(I){
    p.value.vec = rep(NA, 1000)  # Vector to store p-values for each iteration
    for(mc in 1:1000){
      ## sample from the population 
      set.seed(mc)
      sample.data = generate_binary_effect(I=I,linear_coefficient = linear_coefficient, seed=mc)
      sample.data_x = sample.data$dose
      sample.data_y = sample.data$response
      
      # Compute the p-value
      p.value.vec[mc] = mcnemar.sendose(x = sample.data_x,y = sample.data_y,method = method,alternative = "greater than",  Gamma_ave = Gamma_ave, kappa = kappa)$pval
      
    }
    power.gap = mean(p.value.vec<alpha)-power
    return(power.gap)
  }
  
  sample.size = uniroot(f, lower=10, upper=upper.bound,tol=1)$root
  return(ceiling(sample.size))
  
}



simulated_power_binary = function(method=c("mcnemar","dose-mcnemar"),Gamma_seq,linear_coefficient, alpha=0.01, kappa=function(z) z){
  Gamma_leng = length(Gamma_seq)
  power_seq = rep(NA,Gamma_leng)
  for(i in 1:Gamma_leng){
    p.value = rep(NA,1000)
    for(seed in 1:1000){
      sample_data = generate_binary_effect(I=5000, linear_coefficient = linear_coefficient,seed=seed)
      p_value = mcnemar.sendose(x=sample_data$dose,y=sample_data$response, method=method,alternative = "greater than",Gamma_ave = Gamma_seq[i], kappa = kappa)$pval
      p.value[seed] = p_value
    }
    simulated.power = mean(p.value<alpha)
    power_seq[i] = simulated.power
  }
  return(list(method=method,Gamma_seq=Gamma_seq,simulated.power = power_seq, alpha=alpha))
}






sendose.mcnemar.adaptive <- function(
    x, y,
    method1 = c("mcnemar", "dose-mcnemar"),
    method2 = c("mcnemar", "dose-mcnemar"),
    log.p = FALSE, alternative = "greater than",
    Gamma_ave = 1, alpha = 0.05, kappa1 = function(z) z, kappa2 = function(z) z) {
  
  stopifnot(is.matrix(x), is.matrix(y), is.numeric(x), (is.numeric(y) || is.logical(y)))
  stopifnot(ncol(x) == 2, ncol(y) == 2, all(!is.na(x)), all(!is.na(y)))
  stopifnot(all(abs(x[,2] - x[,1]) > .Machine$double.eps))
  stopifnot(is.numeric(Gamma_ave) && Gamma_ave > 0)
  
  alternative <- match.arg(alternative, c("greater than", "less than"))
  
  
  
  
  find_Q_given_correlation <- function(alpha, correlation){
    finderFS <- function(Q, correlation, alpha)
      1 - mvtnorm::pmvnorm(upper = rep(Q, 2), corr = matrix(c(1, correlation, correlation, 1), 2))[1] - alpha
    uniroot(finderFS, c(0, qnorm(1 - alpha / 4)), correlation, alpha)$root
  }
  
  
  
  find.worst.case.cor.obj.p <- function(dose.score.mat, tilde.gamma){
    ## Worst-case correlation rho* of the adaptive procedure, equation (7)
    ## (binary-treatment / McNemar analogue).
    ## Under the Rosenbaum sensitivity bound (3), the least favourable within-pair
    ## treatment-assignment probability is attained at the upper endpoint
    ##     p_plus_i = Gamma_i / (1 + Gamma_i),   Gamma_i = exp(tilde.gamma * |D_i|).
    ## rho* is therefore obtained in closed form by evaluating the correlation of the
    ## two score statistics at this single p_plus; no optimization is performed. The
    ## same p_plus underlies the marginal standardizations returned by
    ## mcnemar.sendose(), so the correlation and the two standardized deviates are
    ## evaluated at one common worst-case configuration, as in equation (7).
    cor.matched.pair.obj.p <- function(dose.score.mat, pi){
      q_i <- dose.score.mat[,2]; s_i <- dose.score.mat[,3]
      cov <- sum(q_i * s_i * pi * (1-pi))
      sd1 <- sqrt(sum(q_i^2 * pi * (1-pi)))
      sd2 <- sqrt(sum(s_i^2 * pi * (1-pi)))
      return(cov / sd1 / sd2)
    }
    abs.dose.diff <- dose.score.mat[,1]
    Gamma_i <- exp(tilde.gamma * abs.dose.diff)
    p.plus  <- Gamma_i / (1 + Gamma_i)                 ## upper Rosenbaum bound, eq. (7)
    list(worst.case.pi          = p.plus,
         worst.case.correlation = cor.matched.pair.obj.p(dose.score.mat, p.plus))
  }
  
  
  first.treatment <- x[,1]; second.treatment <- x[,2]
  first.response <- y[,1]; second.response <- y[,2]
  abs.dose.diff <- abs(first.treatment - second.treatment)
  
  q_i <- mcnemar.q_I.fun(x=x,y=y, method=method1, kappa = kappa1)
  s_i <- mcnemar.q_I.fun(x=x,y=y, method=method2, kappa = kappa2)
  
  # if (alternative == "less than") {
  #  q_i <- -q_i
  #  s_i <- -s_i
  # }
  
  tilde.gamma <- Gamma_to_tilde_gamma_fun(x=x, generalized.Gamma = Gamma_ave)
  
  
  dose.score.mat <- cbind(abs.dose.diff, q_i, s_i)
  
  
  worst.case.res <- find.worst.case.cor.obj.p(dose.score.mat, tilde.gamma)
  worst.case.pi <- worst.case.res$worst.case.pi
  worst.case.correlation <- worst.case.res$worst.case.correlation
  
  test.stat.1 <- mcnemar.Test.stat.fun(x=x, y=y, method=method1, kappa = kappa1)
  test.stat.2 <- mcnemar.Test.stat.fun(x=x, y=y, method=method2, kappa = kappa2)
  
  # if (alternative == "less than") {
  #  test.stat.1 <- -test.stat.1
  #  test.stat.2 <- -test.stat.2
  # }
  
  critical.value <- find_Q_given_correlation(alpha = alpha, correlation = worst.case.correlation)
  
  
  ## use the sendose function to compute the mean and standard deviation
  test.1.res= mcnemar.sendose(x=x,y=y, method = method1, alternative = alternative, Gamma_ave = Gamma_ave, kappa = kappa1)
  test.2.res = mcnemar.sendose(x=x,y=y, method = method2, alternative = alternative, Gamma_ave = Gamma_ave, kappa = kappa2) 
  
  test.1.mean = test.1.res$sense.mean
  test.2.mean = test.2.res$sense.mean
  test.1.sd = sqrt(test.1.res$sense.variance)
  test.2.sd = sqrt(test.2.res$sense.variance)
  
  
  
  std.1 = (test.stat.1 - test.1.mean)/test.1.sd
  std.2 = (test.stat.2 - test.2.mean)/test.2.sd
  ## compute the p value based on the worst-case bivariate normal distribution and whether rejection or not
  if(alternative=="greater than"){
    rejection = ifelse(std.1>=critical.value || std.2 >= critical.value, 1,0)
  }
  if(alternative=="less than"){
    rejection = ifelse(std.1<=-critical.value || std.2 <=-critical.value,1,0)
  }
  
  
  
  return(list(
    rejection = rejection,
    worst.case.correlation = worst.case.correlation,
    worst.case.pi = worst.case.pi,
    critical.value = critical.value,
    std.1 = std.1,
    std.2 = std.2,
    alternative = alternative,
    tilde.gamma = tilde.gamma
  ))
}







############################################################################################
#             Section C.1 in the Supplementary Material, Table S.1                         #
#         Design Sensitivity and simulated powers, Mcnemar and dose-weighted Mcnemar,      #
#                                                                                          #
############################################################################################


## When beta = 3, the design sensitivities 
set.seed(100)
population_data = generate_binary_effect(I=500000, linear_coefficient =3)
mcnemar.design = mcnemar.design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="mcnemar")
dose.weighted.mcnemar.design = mcnemar.design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose-mcnemar")
design.3.effect = c(mcnemar.design,dose.weighted.mcnemar.design)
names(design.3.effect) = c("mcnemar","dose-mcnemar")

## When beta = 5, the design sensitivities 

set.seed(100)
population_data = generate_binary_effect(I=500000, linear_coefficient =5)
mcnemar.design = mcnemar.design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="mcnemar")
dose.weighted.mcnemar.design = mcnemar.design_sensitivity_fun(x=population_data$dose,y=population_data$response,method="dose-mcnemar")
design.5.effect = c(mcnemar.design,dose.weighted.mcnemar.design)
names(design.5.effect) = c("mcnemar","dose-mcnemar")

## Powers

set.seed(100)
mcnemar.3.power      = simulated_power_binary(method="mcnemar",Gamma_seq=c(1,1.5,2,2.5,3,3.5),linear_coefficient = 3)
dose.mcnemar.3.power = simulated_power_binary(method="dose-mcnemar",Gamma_seq=c(1,1.5,2,2.5,3,3.5),linear_coefficient = 3)
mcnemar.5.power      = simulated_power_binary(method="mcnemar",Gamma_seq=c(1,1.5,2,2.5,3,3.5),linear_coefficient = 5)
dose.mcnemar.5.power = simulated_power_binary(method="dose-mcnemar",Gamma_seq=c(1,1.5,2,2.5,3,3.5),linear_coefficient = 5)


design.sensitivity = c(as.numeric(design.3.effect), as.numeric(design.5.effect))
powers = cbind(mcnemar.3.power$simulated.power, dose.mcnemar.3.power$simulated.power, mcnemar.5.power$simulated.power, 
               dose.mcnemar.5.power$simulated.power)

mcnemar.design.power.tb = rbind(design.sensitivity, powers)

mcnemar.design.power.tb







#############################################################################################
# Bahadur-Rosenbaum relative efficiency and sample size, Mcnemar and dose-weighted Mcnemar, # 
#                        Section C in the Supplementary Material ,Table S.2, for \beta = 3  #     
#############################################################################################


set.seed(100)
population_data = generate_binary_effect(I=500000, linear_coefficient =3)
effect.3.mcnemar.to.dose.mcnemar.baha.ratio = rep(NA,4)
effect.3.dose.mcnemar.to.mcnemar.sample.size.ratio.0.01 = rep(NA,4)
effect.3.mcnemar.size.vec = rep(NA,4)
effect.3.dose.mcnemar.size.vec = rep(NA,4)

Gamma_ave_vec = seq(from=1.125,to=1.50,length.out=4)

for(index in 1:4){
  mcnemar_slope = mcnemar.Bhadur.fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="mcnemar")
  dose.mcnemar_slope = mcnemar.Bhadur.fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose-mcnemar")
  effect.3.mcnemar.to.dose.mcnemar.baha.ratio[index] = mcnemar_slope/dose.mcnemar_slope
}

for(index in 1:4){
  ## 
  mcnemar.size = binary_sample_size_fun(alpha=0.01,method="mcnemar", power = 0.95, Gamma_ave = Gamma_ave_vec[index], linear_coefficient = 3)
  effect.3.mcnemar.size.vec[index] = mcnemar.size
  dose.mcnemar.size = binary_sample_size_fun(alpha=0.01,method="dose-mcnemar", power = 0.95, Gamma_ave = Gamma_ave_vec[index],linear_coefficient = 3)
  effect.3.dose.mcnemar.size.vec[index] = dose.mcnemar.size
  effect.3.dose.mcnemar.to.mcnemar.sample.size.ratio.0.01[index] = dose.mcnemar.size/mcnemar.size
  
}

effect.3.table = cbind(effect.3.mcnemar.size.vec, effect.3.dose.mcnemar.size.vec, 
                       effect.3.mcnemar.to.dose.mcnemar.baha.ratio, effect.3.dose.mcnemar.to.mcnemar.sample.size.ratio.0.01)

effect.3.table





#############################################################################################
# Bahadur-Rosenbaum relative efficiency and sample size, Mcnemar and dose-weighted Mcnemar, #
#                     Section C in the Supplementary Material, Table S.2, for \beta = 5     #
#############################################################################################

set.seed(100)
population_data = generate_binary_effect(I=500000, linear_coefficient =5)
effect.5.mcnemar.to.dose.mcnemar.baha.ratio = rep(NA,4)
effect.5.dose.mcnemar.to.mcnemar.sample.size.ratio.0.01 = rep(NA,4)
effect.5.mcnemar.size.vec = rep(NA,4)
effect.5.dose.mcnemar.size.vec = rep(NA,4)

Gamma_ave_vec = seq(from=1.125,to=1.50,length.out=4)

for(index in 1:4){
  mcnemar_slope = mcnemar.Bhadur.fun(x=population_data$dose, y = population_data$response, Gamma_ave = Gamma_ave_vec[index], method="mcnemar")
  dose.mcnemar_slope = mcnemar.Bhadur.fun(x=population_data$dose,y = population_data$response,Gamma_ave = Gamma_ave_vec[index], method="dose-mcnemar")
  effect.5.mcnemar.to.dose.mcnemar.baha.ratio[index] = mcnemar_slope/dose.mcnemar_slope
}

for(index in 1:4){
  
  mcnemar.size = binary_sample_size_fun(alpha=0.01,method="mcnemar", power = 0.95, Gamma_ave = Gamma_ave_vec[index], linear_coefficient = 5)
  effect.5.mcnemar.size.vec[index] = mcnemar.size
  dose.mcnemar.size = binary_sample_size_fun(alpha=0.01,method="dose-mcnemar", power = 0.95, Gamma_ave = Gamma_ave_vec[index],linear_coefficient = 5)
  effect.5.dose.mcnemar.size.vec[index] = dose.mcnemar.size
  effect.5.dose.mcnemar.to.mcnemar.sample.size.ratio.0.01[index] = dose.mcnemar.size/mcnemar.size
  
}

effect.5.table = cbind(effect.5.mcnemar.size.vec, effect.5.dose.mcnemar.size.vec, 
                       effect.5.mcnemar.to.dose.mcnemar.baha.ratio, effect.5.dose.mcnemar.to.mcnemar.sample.size.ratio.0.01)

effect.5.table






###############################################################################
# 
#  Mcnemar and dose-mcnemar, Adaptive testing, effect = 3, Table S.1 in 
#  Section C.1 in the supplementary material, the last column             
#
#
###############################################################################





Gamma_vec = c(1,1.5,2,2.5,3.0,3.5)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    binary.data = generate_binary_effect(I = 5000, linear_coefficient = 3,seed=run)
    res = sendose.mcnemar.adaptive(
      x = binary.data$dose, 
      y = binary.data$response, 
      method1 = "mcnemar", method2 = "dose-mcnemar",alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 50 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.vec[index] = power
}

Gamma_vec
power.vec







###############################################################################
#
#  Mcnemar and dose-mcnemar, Adaptive testing, effect = 5, Table S.1 in 
#  Section C.1 in the supplementary material, the last column
#
#
###############################################################################





Gamma_vec = c(1,1.5,2,2.5,3.0,3.5)
power.vec = rep(NA, length(Gamma_vec))  

for(index in 1:length(Gamma_vec)){  
  rejection = c()
  gamma_start_time = Sys.time()
  Gamma_ave = Gamma_vec[index]
  for(run in 1:1000){
    ## generate the data 
    binary.data = generate_binary_effect(I = 5000, linear_coefficient = 5,seed=run)
    res = sendose.mcnemar.adaptive(
      x = binary.data$dose, 
      y = binary.data$response, 
      method1 = "mcnemar", method2 = "dose-mcnemar",alpha = 0.01,
      Gamma_ave = Gamma_ave
    )
    # cat("std.1", res$std.1, "std.2", res$std.2,"\n", "critical value", res$critical.value,"worst.case.correlation", res$worst.case.correlation,"\n")
    
    rejection = c(rejection, res$rejection)
    
    
    if(run %% 50 == 0){
      cat("current rejection rate", mean(rejection), "current run", run, "at Gamma",Gamma_ave,"\n")
    }
  }
  gamma_end_time = Sys.time()
  gamma_total_time = difftime(gamma_end_time, gamma_start_time, units = "secs")
  
  power = mean(rejection)
  cat("Gamma_ave",Gamma_ave,"rejection rate",power,"\n")
  cat("total_time",gamma_total_time,"secs","\n")
  power.vec[index] = power
}

Gamma_vec
power.vec