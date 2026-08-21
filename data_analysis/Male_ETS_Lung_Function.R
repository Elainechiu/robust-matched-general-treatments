#####################################################################################################
#                                                                                                   #  
#            This file conducts the data analysis for the paper:                                    #
# Towards Robust Matched Observational Studies with General Treatment Types:                        #
#                       Consistency, Efficiency, and Adaptivity                                     #
#                                                                                                   # 
# The codes are divided into two parts: matching and data analysis.                                 #
# This file recovers Table 5 in the main text and Table S.13 in the supplementary                   #
#                                                                                                   #
#                                                                                                   #
# The Data Cleaning and Matching section assumes that the following NHANES 2011-2012                #
# data files are present in the working directory:                                                  #
#                                                                                                   #
#     DEMO_G.xpt, BMX_G.xpt, COTNAL_G.xpt, SPX_G.xpt, MCQ_G.xpt, SMQ_G.xpt                          #
#                                                                                                   #
# They are available, without registration or a data use agreement, from                            #
#  https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2011                    #
# under the sections Demographics Data (DEMO_G), Examination Data (BMX_G, SPX_G),                   #
# Laboratory Data (COTNAL_G), and Questionnaire Data (MCQ_G, SMQ_G).                                #
#                                                                                                   #
# NHANES public-release files are in the public domain and carry no privacy                         #
# restriction. The post-matching data set derived below, male_40_matched_data.csv,                  #
# is included in this repository, so Table 5, Table S.13, and Figure 3 can be                       #
# reproduced without downloading any NHANES file or repeating the matching step.                    #
#####################################################################################################




########################################################################################
#                                                                                      #
#                                                                                      #  
#                          Data Cleaning and Matching                                  #
#                                                                                      #
#                                                                                      #
########################################################################################

# Data Codes from NHANES 2011 - 2012 we verify the distribution of these variables before matching  

# SEQN - Respondent sequence number  

# LBXCOT - Cotinine, Serum (ng/mL)

# RIAGENDR - Gender: 1 male; 2 Female  

# RIDAGEYR - Age in years at screening

# SPXNFEV1 - Baseline FEV 1 (mL)

# SPXNFVC - Baseline FVC (mL)

# BMXWT - weight in (kg)  

# BMXBMI - bmi 

# BMXWAIST - Waist Circumference (cm)

# RIDRETH1 - Race, 3 means white 

# INDFMIN2 - Annual family income 

# INDFMPIR - Income to poverty ratio 

# DMDHREDU - HH ref person's education level

# MCQ010 - Ever been told you have asthma



####### Load Libraries ##############################################################################################

library(MASS)
library(dplyr)
library(moments)
library(designmatch)
library(haven)


##### Load the Data from NHANES organize them into a data set ######################################################

DEMO_G <- read_xpt("DEMO_G.xpt")
BMX_G <- read_xpt("BMX_G.xpt")
COT_G<- read_xpt("COTNAL_G.xpt")
SPX_G <- read_xpt("SPX_G.xpt")
MCQ_G <- read_xpt("MCQ_G.xpt")
SMQ_G <- read_xpt("SMQ_G.xpt")

DEMO <- bind_rows(DEMO_G)
COT = bind_rows(COT_G)
SPX = bind_rows(SPX_G)
BMX = bind_rows(BMX_G)
MCQ = bind_rows(MCQ_G)
SMQ = bind_rows(SMQ_G)


Data_2011_2012<- left_join(DEMO, COT, by="SEQN")
Data_2011_2012<- left_join(Data_2011_2012, SPX, by="SEQN")
Data_2011_2012 = left_join(Data_2011_2012, BMX, by="SEQN")
Data_2011_2012 = left_join(Data_2011_2012, MCQ, by="SEQN")
Data_2011_2012 = left_join(Data_2011_2012, SMQ, by="SEQN")


VarsToKeep = c("LBXCOT","SEQN","SPXNFEV1","RIAGENDR","RIDAGEYR","SPXNFVC",
               "BMXBMI","BMXWAIST","RIDRETH1","INDFMIN2", "INDFMPIR","DMDHREDU","MCQ010")
health_data_2011_2012<- Data_2011_2012[,VarsToKeep]





### Here first checks the variable distribution is the same as the NHANES documentation and rename them  ################### 

## RIAGENDR - Gender
## 1	Male	    4856	
## 2	Female	  4900	
## .	Missing	  0
table(health_data_2011_2012$RIAGENDR)


## RIDAGEYR - Age in years at screening
## 0 to 79	Range of Values	          9393		
## 80	      80 years of age and over	363		
##  .	      Missing	                  0

table(health_data_2011_2012$RIDAGEYR)


## LBXCOT - Cotinine (ng/mL)
## 0.011 to 1700    7378		
## .	Missing	     1057
table(ifelse(health_data_2011_2012$LBXCOT<=1700,1,0))


## SPXNFVC - Baseline FVC (mL)
## 592 to 7863		6652		
##  .	Missing	    843	
summary(health_data_2011_2012$SPXNFVC)
sum(health_data_2011_2012$SPXNFVC>=592 & health_data_2011_2012$SPXNFVC<=7863, na.rm=TRUE)


##  SPXNFEV1 - Baseline FEV 1 (mL)
## 409 to 6923		6652	
## .	Missing	    843	
summary(health_data_2011_2012$SPXNFEV1)
sum(health_data_2011_2012$SPXNFEV1>=409 & health_data_2011_2012$SPXNFEV1<=6923, na.rm=TRUE)



## BMXBMI - Body Mass Index (kg/m**2)

## BMXWAIST - Waist Circumference (cm)

## RIDRETH1 - Race/Hispanic origin
# 1	Mexican American	                  1355	
# 2	Other Hispanic	                    1076	
# 3	Non-Hispanic White	                2973		
# 4	Non-Hispanic Black	                2683		
# 5	Other Race - Including Multi-Racial	1669		
# .	Missing	                             0	

table(health_data_2011_2012$RIDRETH1)



## INDFMPIR - Ratio of family income to poverty
# 0 to 4.99	Range of Values	            7656		
# 5	Value greater than or equal to 5.00	1260		
# .	Missing	                             840
summary(health_data_2011_2012$INDFMPIR)
sum(health_data_2011_2012$INDFMPIR<=4.99,na.rm=TRUE)


## DMDHREDU - HH ref person's education level
## 1	Less Than 9th Grade	                                 966		
## 2	9-11th Grade (Includes 12th grade with no diploma)	1430		
## 3	High School Grad/GED or Equivalent	                2003		
## 4	Some College or AA degree	                          2665		
## 5	College Graduate or above	                          2294		
## 7	Refused	                                               3		
## 9	Don't Know	                                          33	
## .	Missing	                                             362


table(health_data_2011_2012$DMDHREDU)



## MCQ010 Have been told to have asthma 
# 1	Yes	                                                1448		
# 2	No	                                                7906	  	
# 7	Refused	                                             1	
# 9	Don't know	                                         8		
# .	Missing	                                             1	

table(health_data_2011_2012$MCQ010)

health_data_2011_2012$sample = complete.cases(health_data_2011_2012)




complete_data = subset(health_data_2011_2012, health_data_2011_2012$sample==T)



complete_data <- complete_data %>%
  rename(
    ID = SEQN,
    gender = RIAGENDR,
    age = RIDAGEYR,
    smoke = LBXCOT,
    FEV1 = SPXNFEV1,
    FVC = SPXNFVC,
    RACE = RIDRETH1,
    familyincome = INDFMIN2,
    incometopoverty = INDFMPIR,
    hheducation = DMDHREDU,
    asthma = MCQ010,
    bmi = BMXBMI,
    waist_circumference = BMXWAIST
  )




## only focus on male data older than 40, also, if there is any random variable >7, means the answer is do not know or refuse
## also get rid of them
male.data = subset(complete_data, complete_data$gender==1 & complete_data$age>=40 & complete_data$hheducation<7 & complete_data$asthma <7)
male.data$ratio = male.data$FEV1/male.data$FVC
male.data$asthma = ifelse(male.data$asthma==1,1,0)
male.data$white = ifelse(male.data$RACE==3,1,0)
male.data$some.college = ifelse(male.data$hheducation>=4,1,0)


# Box-Cox Transformation, this transformation is monotonic, does not change the direction of analysis
bc <- boxcox(male.data$smoke ~ 1, plotit = T)

lambda <- bc$x[which.max(bc$y)]
# lambda -0.02
print(lambda)
male.data$bc_smoke <- (male.data$smoke^lambda - 1) / lambda
skewness(male.data$bc_smoke, na.rm = TRUE)
hist(male.data$bc_smoke, main = "Box-Cox-transformed smoke")





################################ conduct the matching ############################################################


# Here use the linear regression to define a balancing score and define the distance
balancing_smoke_male =  lm(bc_smoke~age+incometopoverty, data = male.data)
summary(balancing_smoke_male)

balancing_score_male = fitted.values(balancing_smoke_male)

non_binary_matching <- function(balancing.scores, treatment.dose, ID, 
                                max.penalty=99999, multiplier=100, epsilon = 1e-6) {
  # Check inputs
  if (!is.numeric(balancing.scores)) {
    stop("The balancing scores must be numeric.")
  }
  if (!is.numeric(treatment.dose)) {
    stop("The treatment.dose must be numeric.")
  }
  if (length(balancing.scores) != length(treatment.dose)) {
    stop("The balancing scores and treatment dose must have the same length.")
  }
  if (!is.numeric(max.penalty) || max.penalty <= 0) {
    stop("max.penalty should be a positive numeric value.")
  }
  if (!is.numeric(multiplier) || multiplier <= 0) {
    stop("multiplier should be a positive numeric value.")
  }
  if (length(ID) != length(treatment.dose)) {
    stop("The length of ID must equal the length of treatment.dose.")
  }
  
  data_length <- length(balancing.scores)
  
  distance.mat <- matrix(0, nrow = data_length, ncol = data_length,
                         dimnames = list(ID, ID))
  
  
  for (i in 1:data_length) {
    for (j in i:data_length) {
      if (i == j) {
        # Set diagonal entries to max.penalty to avoid self-matching.
        distance.mat[i, j] <- max.penalty
      } else {
        # Get treatment doses and balancing scores for observations i and j.
        zi <- treatment.dose[i]
        zj <- treatment.dose[j]
        linear.xi <- balancing.scores[i]
        linear.xj <- balancing.scores[j]
        
        # If the treatment doses are equal, assign max.penalty.
        if (zi == zj) {
          pair.distance <- max.penalty
        } else {
          pair.distance <- multiplier * (((linear.xi - linear.xj)^2 + epsilon) / ((zi - zj)^2))
        }
        # Fill both [i,j] and [j,i] to ensure symmetry.
        distance.mat[i, j] <- pair.distance
        distance.mat[j, i] <- pair.distance
      }
    }
  }
  
  return(distance.mat)
}


distance = non_binary_matching(balancing.scores = balancing_score_male, ID=male.data$ID,treatment.dose = male.data$bc_smoke)
distance_mat = as.matrix(distance)



total_pairs = 500
subset_weight = 1 
t_max = 60*5 
solver = "highs"
approximate = 1 


solver = list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, 
              trace_cplex = 0)

## Match                  
out = nmatch(dist_mat = distance_mat, subset_weight = subset_weight, total_pairs = total_pairs, solver = solver)    



## Indices of the treated units and matched controls
id_1 = out$id_1
id_2 = out$id_2


group1.id = id_1
group2.id = id_2



health.matched.pair.df = data.frame(pair.id = 1:length(group1.id),
                                    group1.dose=male.data$bc_smoke[group1.id], 
                                    group2.dose=male.data$bc_smoke[group2.id],
                                    group1.response=male.data$ratio[group1.id], 
                                    group2.response=male.data$ratio[group2.id])


high.dose = c()
high.dose.id = c()
low.dose  = c()
low.dose.id = c()
high.response = c()
low.response = c()

for(i in 1:length(group1.id)){
  high = max(health.matched.pair.df$group1.dose[i],
             health.matched.pair.df$group2.dose[i])
  low = min(health.matched.pair.df$group1.dose[i],
            health.matched.pair.df$group2.dose[i])
  high.dose = c(high.dose, high)
  low.dose = c(low.dose, low)
  if(high ==health.matched.pair.df$group1.dose[i]){
    high.response = c(high.response,health.matched.pair.df$group1.response[i])
    low.response = c(low.response, health.matched.pair.df$group2.response[i])
    high.dose.id = c(high.dose.id,group1.id[i])
    low.dose.id = c(low.dose.id,group2.id[i])
  }else{
    
    high.response = c(high.response,health.matched.pair.df$group2.response[i])
    low.response = c(low.response, health.matched.pair.df$group1.response[i])
    low.dose.id = c(low.dose.id,group1.id[i])
    high.dose.id = c(high.dose.id,group2.id[i])
  }
}



health.matched.pair.df$low.dose = low.dose
health.matched.pair.df$high.dose = high.dose
health.matched.pair.df$low.response = low.response
health.matched.pair.df$high.response = high.response
health.matched.pair.df$low.dose.id = low.dose.id
health.matched.pair.df$high.dose.id = high.dose.id   

health.matched.pair.df$low.dose.age = male.data$age[low.dose.id]
health.matched.pair.df$high.dose.age = male.data$age[high.dose.id]
health.matched.pair.df$low.dose.bmi = male.data$bmi[low.dose.id]
health.matched.pair.df$high.dose.bmi=male.data$bmi[high.dose.id]
health.matched.pair.df$low.dose.waist_circumference=male.data$waist_circumference[low.dose.id]
health.matched.pair.df$high.dose.waist_circumference=male.data$waist_circumference[high.dose.id]

health.matched.pair.df$low.dose.incometopoverty=male.data$incometopoverty[low.dose.id]
health.matched.pair.df$high.dose.incometopoverty=male.data$incometopoverty[high.dose.id]

health.matched.pair.df$low.dose.some.college = male.data$some.college[low.dose.id]
health.matched.pair.df$high.dose.some.college = male.data$some.college[high.dose.id]

health.matched.pair.df$low.dose.asthma = male.data$asthma[low.dose.id]
health.matched.pair.df$high.dose.asthma = male.data$asthma[high.dose.id]


health.matched.pair.df$low.dose.white = male.data$white[low.dose.id]
health.matched.pair.df$high.dose.white = male.data$white[high.dose.id]


summary(health.matched.pair.df)


#### Assess balancing before storing the data 
mom_covs = male.data[,c("age","bmi","waist_circumference","incometopoverty","asthma","white","some.college")]
a = apply(mom_covs[low.dose.id, ], 2, mean)
b = apply(mom_covs[high.dose.id, ], 2, mean)

c = apply(mom_covs, 2, sd)
tab = round(cbind(a, b, (a-b)/c), 2)
colnames(tab) = c("Mean 1", "Mean 2", "Std.Diffs")
rownames(tab) = c("age","bmi","waist_cicumsference","incometopoverty","asthma","white","some.college")
tab


## store into a csv

write.csv(x=health.matched.pair.df, file="male_40_matched_data.csv")
















#############################################################################################################################
#                                                                                                                           #
#                                                                                                                           #
#                                           Data Analysis                                                                   #
#                                                                                                                           #
#                                                                                                                           #
#                                                                                                                           #
#                                                                                                                           #
#############################################################################################################################




################## Recover the post-matching balance, Table S.13 in Supplementary ##################################################

source("matched_sensitivity_analysis_functions.R")
health.matched.pair.df = read.csv("data_analysis/male_40_matched_data.csv")
x = cbind(health.matched.pair.df$group1.dose, health.matched.pair.df$group2.dose)
y = cbind(health.matched.pair.df$group1.response, health.matched.pair.df$group2.response)

## recheck the matching quality 

#### Assess balancing #########
variable = c("age","white","incometopoverty","some.college","bmi","waist_circumference","asthma")
low.dose.mean = rep(NA,7)
high.dose.mean = rep(NA,7)
sd = rep(NA,7)
Std.diff = rep(NA,7)
for(var in 1:7){
  var.name = variable[var]
  low.dose.column = health.matched.pair.df[,paste0("low.dose.",var.name)]
  high.dose.column = health.matched.pair.df[,paste0("high.dose.",var.name)]
  combine.column = c(low.dose.column, high.dose.column)
  low.dose.mean[var] = mean(low.dose.column)
  high.dose.mean[var] = mean(high.dose.column)
  sd[var] = sd(combine.column)
  Std.diff[var] = (mean(low.dose.column)-mean(high.dose.column))/sd(combine.column)
  
  
}

balance.table = round(cbind(low.dose.mean,high.dose.mean,Std.diff),2)
rownames(balance.table) = variable


### Reproduces Table S.13 in the supplementary material 
balance.table


# Plot the data points
plot(x, y,
     xlab = "Dose",
     ylab = "Response",
     pch  = 16, 
     cex  = 0.5,
     col  = rgb(0, 0, 0, 0.3))  



# Generate Figure 3 in the main text
fit_spline <- smooth.spline(x, y, spar = 1.5)
lines(fit_spline, col = "blue", lwd = 2)



################################# Sensitivity Analysis for the Data #########################################################

## We note that we set the alternative to be less than means we a smaller test statistic is an evidence against the null
## the alternative hypothesis posits that a higher dose (higher smoking) causes a decrease in lung function 
## therefore, the sign*(dose.diff*response.diff)*score is considered.


source("matched_sensitivity_analysis_functions.R")
health.matched.pair.df = read.csv("data_analysis/male_40_matched_data.csv")
x = cbind(health.matched.pair.df$group1.dose, health.matched.pair.df$group2.dose)
y = cbind(health.matched.pair.df$group1.response, health.matched.pair.df$group2.response)



Gamma_ave_vec = seq(1,1.8,length.out=9)
wilcoxon.p = rep(NA,9)
dose.weighted.p = rep(NA,9)
U878.p = rep(NA,9)
dose.weighted.U878.p = rep(NA,9)
for(index in 1:9){
  Gamma_ave = Gamma_ave_vec[index]
  wilcoxon.p[index] = sendose(x=x,y=y,method = "wilcoxon", alternative = "less than",Gamma_ave = Gamma_ave)$pval
  dose.weighted.p[index] = sendose(x=x,y=y,method = "dose.weighted", alternative = "less than",Gamma_ave = Gamma_ave)$pval
  U878.p[index] = sendose(x=x,y=y,method = "U", alternative = "less than",Gamma_ave = Gamma_ave,m=8,m1=7,m2=8)$pval
  dose.weighted.U878.p[index] = sendose(x=x,y=y, method="dose.weighted.U", alternative = "less than", Gamma_ave = Gamma_ave, m=8,m1=7,m2=8)$pval
  
}

p.value.tb = cbind(wilcoxon.p,dose.weighted.p,U878.p, dose.weighted.U878.p)
colnames(p.value.tb) = c("wilcoxon","dose-weighted","(8,7,8)","dose-weighted.(8,7,8)")
rownames(p.value.tb) = Gamma_ave_vec

p.value.tb = round(p.value.tb,3)



### Find the sensitivity value 
find_sensitivity_value <- function(x, y, method = c("U","wilcoxon","dose.weighted","response.polynomial","dose.weighted.U"), m, m1, m2, 
                                   alternative = "less than", rank.order = 3,
                                   target = 0.05, 
                                   lower = 0.01, upper = 10, init = 1.0, kappa = function(z) z) {
  objective_function <- function(Gamma_ave) {
    result <- sendose(x = x, y = y, method = method, m = m, m1 = m1, 
                      m2 = m2, alternative = alternative, Gamma_ave = Gamma_ave, rank.order=rank.order, kappa = kappa)$pval
    return((result - target)^2)
  }
  opt <- optim(par = init, fn = objective_function, method = "Brent", 
               lower = lower, upper = upper)
  
  # Return the estimated Gamma_ave value
  Gamma_ave = opt$par
  
  tilde.gamma = Gamma_to_tilde_gamma_fun(x=x, generalized.Gamma = Gamma_ave)
  abs.dose.diff = abs(x[,1]-x[,2])
  smallest_Gamma_i = exp(tilde.gamma*min(abs.dose.diff))
  largest_Gamma_i = exp(tilde.gamma*max(abs.dose.diff))
  return(list(Gamma_ave=Gamma_ave, smallest_Gamma_i=smallest_Gamma_i, largest_Gamma_i=largest_Gamma_i))
}


wilcoxon.sen.value     = find_sensitivity_value(x=x, y=y, method="wilcoxon", alternative = "less than")$Gamma_ave
dose.weighted.sen.value = find_sensitivity_value(x=x, y=y, method="dose.weighted", alternative = "less than")$Gamma_ave
U878.sen.value         = find_sensitivity_value(x=x, y=y, method="U", m=8,m1=7,m2=8, alternative = "less than")$Gamma_ave
dose.U878.sen.value    = find_sensitivity_value(x=x,y=y, method="dose.weighted.U", m=8,m1=7, m2=8, alternative = "less than")$Gamma_ave
Table.5 = rbind(p.value.tb,round(c(wilcoxon.sen.value,dose.weighted.sen.value,U878.sen.value,dose.U878.sen.value),2))


## First generate the sensitivity value of a single test in Table 5 in the main text 
Table.5





## Here conducts the sensitivity value calculation for the adaptive test 


sendose.adaptive <- function(
    x, y,
    method1 = c("wilcoxon", "dose.weighted", "response.polynomial", "U", "dose.weighted.U"),
    method.1.rank.order = 2, method.1.m = 2, method.1.m1 = 2, method.1.m2 = 2,
    method2 = c("wilcoxon", "dose.weighted", "response.polynomial", "U", "dose.weighted.U"),
    method.2.rank.order = 2, method.2.m = 2, method.2.m1 = 2, method.2.m2 = 2,
    log.p = FALSE, alternative = "greater than",
    Gamma_ave = 1, alpha = 0.05, kappa1 = function(z) z, kappa2 = function(z) z) {
  
  stopifnot(is.matrix(x), is.matrix(y), is.numeric(x), (is.numeric(y) || is.logical(y)))
  stopifnot(ncol(x) == 2, ncol(y) == 2, all(!is.na(x)), all(!is.na(y)))
  stopifnot(all(abs(x[,2] - x[,1]) > .Machine$double.eps))
  stopifnot(is.numeric(Gamma_ave) && Gamma_ave > 0)
  
  alternative <- match.arg(alternative, c("greater than", "less than"))
  
  required_functions <- c("q_I.fun", "Test.stat.fun")
  missing_functions <- c()
  
  for(func in required_functions) {
    if (!exists(func, envir = globalenv()) && !exists(func, envir = parent.frame())) {
      missing_functions <- c(missing_functions, func)
    }
  }
  
  if(length(missing_functions) > 0) {
    stop("Missing required functions: ", paste(missing_functions, collapse = ", "))
  }
  
  find_Q_given_correlation <- function(alpha, correlation){
    finderFS <- function(Q, correlation, alpha)
      1 - mvtnorm::pmvnorm(upper = rep(Q, 2), corr = matrix(c(1, correlation, correlation, 1), 2))[1] - alpha
    uniroot(finderFS, c(0, qnorm(1 - alpha / 4)), correlation, alpha)$root
  }
  
  
  
  find.worst.case.cor.obj.p <- function(dose.score.mat, tilde.gamma){
    ## Worst-case correlation rho* of the adaptive procedure, equation (7) in the main text.
    ## Under the Rosenbaum sensitivity bound (3), the least favourable within-pair
    ## treatment-assignment probability is attained at the upper endpoint
    ##     p_plus_i = Gamma_i / (1 + Gamma_i),   Gamma_i = exp(tilde.gamma * |D_i|).
    ## rho* is therefore obtained in closed form by evaluating the correlation of the
    ## two score statistics at this single p_plus; no optimization is performed. The
    ## same p_plus underlies the marginal standardizations returned by sendose(), so
    ## the correlation and the two standardized deviates are evaluated at one common
    ## worst-case configuration, as in equation (7).
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
  
  q_i <- q_I.fun(x, y, method1, method.1.rank.order, method.1.m, method.1.m1, method.1.m2, kappa = kappa1)
  s_i <- q_I.fun(x, y, method2, method.2.rank.order, method.2.m, method.2.m1, method.2.m2, kappa = kappa2)
  
  # if (alternative == "less than") {
  #  q_i <- -q_i
  #  s_i <- -s_i
  # }
  
  tilde.gamma <- Gamma_to_tilde_gamma_fun(x=x, generalized.Gamma = Gamma_ave)
  
  
  dose.score.mat <- cbind(abs.dose.diff, q_i, s_i)
  
  
  worst.case.res <- find.worst.case.cor.obj.p(dose.score.mat, tilde.gamma)
  worst.case.pi <- worst.case.res$worst.case.pi
  worst.case.correlation <- worst.case.res$worst.case.correlation
  
  test.stat.1 <- Test.stat.fun(x, y, method1, method.1.rank.order, method.1.m, method.1.m1, method.1.m2, kappa = kappa1)
  test.stat.2 <- Test.stat.fun(x, y, method2, method.2.rank.order, method.2.m, method.2.m1, method.2.m2, kappa = kappa2)
  
  # if (alternative == "less than") {
  #  test.stat.1 <- -test.stat.1
  #  test.stat.2 <- -test.stat.2
  # }
  
  critical.value <- find_Q_given_correlation(alpha = alpha, correlation = worst.case.correlation)
  
  
  ## use the sendose function to compute the mean and standard deviation
  test.1.res= sendose(x=x,y=y, method = method1, rank.order = method.1.rank.order, m = method.1.m,  m1 = method.1.m1, m2 =method.1.m2, alternative = alternative, Gamma_ave = Gamma_ave, kappa = kappa1)
  test.2.res = sendose(x=x,y=y, method = method2, rank.order = method.2.rank.order, m = method.2.m, m1 = method.2.m1, m2 =method.2.m2, alternative = alternative, Gamma_ave = Gamma_ave, kappa = kappa2) 
  
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

## conduct the adaptive test at alpha = 0.05, and we remember to set the alternative to "less than"
## The first pair is Wilcoxon and the dose weighted wilcoxon signed rank test 
Gamma_ave = c(1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8)
rejection = numeric(9)
for(i in 1:9){
  res = sendose.adaptive(
    x = x, 
    y = y, 
    method1 = "wilcoxon", method2 = "dose.weighted", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.05,
    Gamma_ave = Gamma_ave[i], alternative="less than"
  )
  rejection[i] = res$rejection
}

## Table 5 in the main text, the adaptive test for the wilcoxon and dose-weighted wilcoxon column
rejection



## conduct the adaptive test at alpha = 0.05 for the second pair U test and dose-weighted U-test

for(i in 1:9){
  res = sendose.adaptive(
    x = x, 
    y = y, 
    method1 = "U", method2 = "dose.weighted.U", method.1.m = 8, method.1.m1=7, method.1.m2 = 8, method.2.m=8, method.2.m1=7, method.2.m2=8,alpha = 0.05,
    Gamma_ave = Gamma_ave[i], alternative="less than"
  )
  rejection[i] = res$rejection
}


## Table 5 in the main text, the adaptive test for the U and dose-weighted U column
rejection




################## Sensitivity Value of the Adaptive Tests ##################################################
##
## The adaptive test rejects H_0 at level alpha if either standardized deviate falls beyond the
## critical value Q(alpha, rho*), where rho* is the worst-case correlation of the two component
## statistics evaluated at the common worst-case treatment-assignment probabilities. For the
## "less than" alternative, the test rejects if and only if min(std.1, std.2) <= -Q(alpha, rho*).
## The attained significance level (worst-case p-value) of the adaptive test is therefore
##     p(Gamma_bar) = 1 - Phi_2(Q*, Q*; rho*),   with   Q* = -min(std.1, std.2),
## where Phi_2 denotes the bivariate normal distribution function with correlation rho*.
## The sensitivity value of the adaptive test is the value of Gamma_bar at which p(Gamma_bar)
## equals the significance level 0.05; since p(Gamma_bar) is monotonically increasing in
## Gamma_bar, the sensitivity value is obtained by root finding/bisection method 
###########################################################################################################


## Compute the worst-case p-value of the adaptive test at a given Gamma_bar
sendose.adaptive.pval <- function(
    x, y,
    method1 = c("wilcoxon", "dose.weighted", "response.polynomial", "U", "dose.weighted.U"),
    method.1.rank.order = 2, method.1.m = 2, method.1.m1 = 2, method.1.m2 = 2,
    method2 = c("wilcoxon", "dose.weighted", "response.polynomial", "U", "dose.weighted.U"),
    method.2.rank.order = 2, method.2.m = 2, method.2.m1 = 2, method.2.m2 = 2,
    alternative = "greater than",
    Gamma_ave = 1, alpha = 0.05,
    kappa1 = function(z) z, kappa2 = function(z) z) {
  
  res <- sendose.adaptive(
    x = x, y = y,
    method1 = method1,
    method.1.rank.order = method.1.rank.order,
    method.1.m = method.1.m, method.1.m1 = method.1.m1, method.1.m2 = method.1.m2,
    method2 = method2,
    method.2.rank.order = method.2.rank.order,
    method.2.m = method.2.m, method.2.m1 = method.2.m1, method.2.m2 = method.2.m2,
    alternative = alternative,
    Gamma_ave = Gamma_ave, alpha = alpha,
    kappa1 = kappa1, kappa2 = kappa2
  )
  
  ## The attained level corresponds to the more extreme of the two standardized deviates
  if (res$alternative == "less than") {
    Q.star <- -min(res$std.1, res$std.2)
  } else {
    Q.star <- max(res$std.1, res$std.2)
  }
  
  ## Evaluate the worst-case bivariate normal at the attained deviate
  rho <- res$worst.case.correlation
  pval <- 1 - mvtnorm::pmvnorm(upper = rep(Q.star, 2),
                               corr = matrix(c(1, rho, rho, 1), 2))[1]
  return(as.numeric(pval))
}


## Find the sensitivity value of the adaptive test at the target significance level
find_adaptive_sensitivity_value <- function(
    x, y,
    method1, method2,
    method.1.rank.order = 2, method.1.m = 2, method.1.m1 = 2, method.1.m2 = 2,
    method.2.rank.order = 2, method.2.m = 2, method.2.m1 = 2, method.2.m2 = 2,
    alternative = "less than",
    target = 0.05,
    lower = 1, upper = 10,
    kappa1 = function(z) z, kappa2 = function(z) z) {
  
  objective_function <- function(Gamma_ave) {
    sendose.adaptive.pval(
      x = x, y = y,
      method1 = method1,
      method.1.rank.order = method.1.rank.order,
      method.1.m = method.1.m, method.1.m1 = method.1.m1, method.1.m2 = method.1.m2,
      method2 = method2,
      method.2.rank.order = method.2.rank.order,
      method.2.m = method.2.m, method.2.m1 = method.2.m1, method.2.m2 = method.2.m2,
      alternative = alternative,
      Gamma_ave = Gamma_ave,
      kappa1 = kappa1, kappa2 = kappa2
    ) - target
  }
  
  Gamma_ave <- uniroot(objective_function, lower = lower, upper = upper,
                       tol = .Machine$double.eps^0.5)$root
  
  ## Report the corresponding gamma and the implied range of the per-pair Gamma_i
  tilde.gamma <- Gamma_to_tilde_gamma_fun(x = x, generalized.Gamma = Gamma_ave)
  abs.dose.diff <- abs(x[, 1] - x[, 2])
  smallest_Gamma_i <- exp(tilde.gamma * min(abs.dose.diff))
  largest_Gamma_i <- exp(tilde.gamma * max(abs.dose.diff))
  
  return(list(Gamma_ave = Gamma_ave,
              tilde.gamma = tilde.gamma,
              smallest_Gamma_i = smallest_Gamma_i,
              largest_Gamma_i = largest_Gamma_i))
}


## Sensitivity value of the adaptive test combining the Wilcoxon and dose-weighted Wilcoxon tests
adaptive.wilcoxon.sen.value <- find_adaptive_sensitivity_value(
  x = x, y = y,
  method1 = "wilcoxon", method2 = "dose.weighted",
  method.1.m = 8, method.1.m1 = 7, method.1.m2 = 8,
  method.2.m = 8, method.2.m1 = 7, method.2.m2 = 8,
  alternative = "less than"
)$Gamma_ave

## Sensitivity value of the adaptive test combining the U and dose-weighted U tests
adaptive.U.sen.value <- find_adaptive_sensitivity_value(
  x = x, y = y,
  method1 = "U", method2 = "dose.weighted.U",
  method.1.m = 8, method.1.m1 = 7, method.1.m2 = 8,
  method.2.m = 8, method.2.m1 = 7, method.2.m2 = 8,
  alternative = "less than"
)$Gamma_ave

## Sensitivity values of the two adaptive tests reported in Table 5 in the main text
round(c(adaptive.wilcoxon.sen.value, adaptive.U.sen.value), 2)


## Convert the sensitivity values of the component and adaptive tests to the gamma scale
sensitivity.values <- c(wilcoxon.sen.value, dose.weighted.sen.value,
                        U878.sen.value, dose.U878.sen.value,
                        adaptive.wilcoxon.sen.value, adaptive.U.sen.value)
gamma.sensitivity.values <- sapply(sensitivity.values, function(G) {
  Gamma_to_tilde_gamma_fun(x = x, generalized.Gamma = G)
})
sensitivity.conversion.tb <- rbind(sensitivity.values, gamma.sensitivity.values)
rownames(sensitivity.conversion.tb) <- c("Gamma_bar", "gamma")
colnames(sensitivity.conversion.tb) <- c("wilcoxon", "dose-weighted", "(8,7,8)",
                                         "dose-weighted.(8,7,8)",
                                         "adaptive.wilcoxon", "adaptive.U")
round(sensitivity.conversion.tb, 3)



################## Conversion of Generalized Gamma to the Rate Parameter gamma ##############################
##
## For each generalized sensitivity parameter Gamma_bar used in the data analysis,
## this section recovers the underlying rate parameter gamma satisfying
##     (1/I) * sum_i exp(gamma * |D_i|) = Gamma_bar,
## where D_i is the within-pair difference in Box-Cox-transformed cotinine.
## The implied per-pair sensitivity parameters Gamma_i = exp(gamma * |D_i|) are
## also reported at their minimum and maximum over the I = 500 matched pairs.
## See Proposition 3.2 in the main text for the one-to-one correspondence between 
## Gamma and gamma 

source("matched_sensitivity_analysis_functions.R")
health.matched.pair.df <- read.csv("data_analysis/male_40_matched_data.csv")
x <- cbind(health.matched.pair.df$group1.dose, health.matched.pair.df$group2.dose)

abs.dose.diff <- abs(x[, 1] - x[, 2])

## The grid of generalized Gamma values used in the sensitivity analysis
Gamma_ave_vec <- seq(1, 1.8, length.out = 9)

## Convert each Gamma_bar to gamma and the implied range of Gamma_i
gamma.vec <- sapply(Gamma_ave_vec, function(G) {
  Gamma_to_tilde_gamma_fun(x = x, generalized.Gamma = G)
})

Gamma.conversion.tb <- cbind(
  Gamma_ave_vec,
  gamma.vec,
  exp(gamma.vec * min(abs.dose.diff)),
  exp(gamma.vec * max(abs.dose.diff))
)
colnames(Gamma.conversion.tb) <- c("Gamma_bar", "gamma", "min Gamma_i", "max Gamma_i")

round(Gamma.conversion.tb, 3)




