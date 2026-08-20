# Load necessary library, this library is for the U test statistic defined in Rosenbaum 2011
require(sensitivitymv)
require(mvtnorm)
require("ggplot2")


#' Convert Generalized Gamma to Tilde Gamma
#'
#' This function converts a generalized Gamma value to a tilde gamma, taking 
#' into account the treatment dose information. For the relationship between 
#' generalized Gamma and tilde gamma, see Heng (2019).
#' @param x this is a I by 2 matrix having the treatment doses for I pairs  
#' @param generalized.Gamma The generalized Gamma value.
#' @param tol Tolerance level for the root finding. Defaults to the 0.9 power
#' of the machine tolerance.
#' @keywords Gamma conversion
Gamma_to_tilde_gamma_fun <- function(x, generalized.Gamma, tol=.Machine$double.eps^0.9) {
  if (generalized.Gamma == 1) {
    return(0)
  }
  
  first.treatment  = x[,1]
  second.treatment = x[,2]
  treatment.diff = first.treatment - second.treatment
  I <- length(treatment.diff)
  M <- max(abs(treatment.diff))
  m <- min(abs(treatment.diff))
  
  lower.bound <- min(log(generalized.Gamma) / M, log(generalized.Gamma) / m)
  upper.bound <- max(log(generalized.Gamma) / m, log(generalized.Gamma) / M)
  
  f <- function(x) I * generalized.Gamma - sum(exp(x * abs(treatment.diff)))
  if(lower.bound==upper.bound){
    root = lower.bound
  }else{
    root <- suppressWarnings(uniroot(f, lower = lower.bound, upper = upper.bound, tol = tol))$root
  }
  return(root)
}




#' Compute Test Statistic
#'
#' This function computes the test statistic based on the ranking option. This 
#' function can compute the usual signed rank test statistic and weighted rank 
#' test statistic. The scores based on the ranks are scaled so that the sum of 
#' scores will always be one for a given data.
#' To see the scores q_i for a given data set, check the function q_I.fun. 
#' The U test statistic with (m,m1,m2)=(2,2,2) is a constant multiple of 
#' the wilcoxon test.
#' @param x the I by 2 matrix containing the treatment dose information
#' @param y the I by 2 matrix containing the response information
#' @param m The m in the u test statistic specification.
#' @param m1 The m1 in the u test statistic specification
#' @param m2 The m2 in the u test statistic specification
#' @param rank.order The order of the polynomial rank
#' @param method The ranking option to use. 
#' @param kappa Monotone dose transformation applied before ranking. Defaults to identity.
#' @keywords test statistic


Test.stat.fun <- function(x,y, method = c("wilcoxon","dose.weighted","response.polynomial","U", "dose.weighted.U"),rank.order=2,m=2,m1=2,m2=2,kappa=function(z) z) {
  first.treatment = x[,1]
  second.treatment= x[,2]
  first.response  = y[,1]
  second.response = y[,2]
  treatment.diff = first.treatment - second.treatment
  I <- length(treatment.diff)
  response.diff = first.response  - second.response
  scaled.rank.treatment = rank(abs(kappa(first.treatment) - kappa(second.treatment)))/I
  scaled.rank.response  = rank(abs(response.diff))/I
  
  if (method == "wilcoxon") {
    scores  = scaled.rank.response
  }
  
  if(method=="dose.weighted"){
    scores <- scaled.rank.response*scaled.rank.treatment
  }
  if(method=="response.polynomial"){
    scores <- scaled.rank.response^rank.order
  }
  
  if(method=="U"){
    scores = multrnks(rk=scaled.rank.response*I,m1=m1,m2=m2,m=m)
  }
  if(method=="dose.weighted.U"){
    scores = multrnks(rk=scaled.rank.response*scaled.rank.treatment*I, m1=m1,m2=m2,m=m)
  }
  
  sign <- ifelse(treatment.diff*response.diff > 0, 1, 0)
  test.stat <- as.numeric(sign %*% scores)
  return(test.stat)
}


q_I.fun <- function(x,y, method = c("wilcoxon","dose.weighted","response.polynomial","U","dose.weighted.U"),rank.order=2,m=2,m1=2,m2=2,kappa=function(z) z) {
  first.treatment = x[,1]
  second.treatment= x[,2]
  first.response  = y[,1]
  second.response = y[,2]
  treatment.diff = first.treatment - second.treatment
  I <- length(treatment.diff)
  response.diff = first.response  - second.response
  scaled.rank.treatment = rank(abs(kappa(first.treatment) - kappa(second.treatment)))/I
  scaled.rank.response  = rank(abs(response.diff))/I
  
  if (method == "wilcoxon") {
    scores  = scaled.rank.response
  }
  
  if(method=="dose.weighted"){
    scores <- scaled.rank.response*scaled.rank.treatment
    
  }
  if(method=="response.polynomial"){
    scores <- scaled.rank.response^rank.order
  }
  
  if(method=="U"){
    scores <- multrnks(rk=scaled.rank.response*I,m=m,m1=m1,m2=m2)
  }
  if(method=="dose.weighted.U"){
    scores <- multrnks(rk=scaled.rank.response*scaled.rank.treatment*I, m=m,m1=m1,m2=m2)
  }
  
  return(scores)
}



#' Compute the upper bound of p value considering Gamma_i, either with
#' normal approximation or with permutation. If the exact option is set FALSE
#' expectation, variance, p value are obtained with 
#' monte carlo simulations; if the exact option is set TRUE, expectation, variance
#' ,p value and are obtained with normal approximation.
#' @param x A I by 2 matrix containing the treatment dose information.
#' @param y A I by 2 matrix containing the response information
#' @param alternative specifies the direction of alternative hypothesis, if the 
#' alternative is greater than, it tests if higher treatment dose increases
#' the response; if the alternative is set "less than", it tests if higher dose
#' decreases responses, for two-sided hypothesis, the algorithm will compute 
#' two p values corresponding to two one-sided hypothesis, and the resulting 
#' p value is twice of the minimum. 
#' @param Gamma_ave is the generalized Gamma, default to one
#' @param exact is false for normal approximated p value; true for permutation p. 
#' default to false. 
#' @param mc.iteration how many times of permutation is conducted for estimating
#' the p value without normal approximation. Default to 1,000.  
#' @param kappa Monotone dose transformation applied before ranking. Defaults to identity.

#' @export
sendose <-function(x, y,method=c("wilcoxon","dose.weighted","response.polynomial","U","dose.weighted.U"),rank.order=2,m=2,m1=2,m2=2,alternative=c("less than","greater than","two.sided"),Gamma_ave=1,exact=FALSE,mc.iteration=1000,log.p=FALSE,kappa=function(z) z){
  ### Basic check on input
  # Check on data
  stopifnot(is.matrix(x),is.matrix(y),
            is.numeric(x),(is.numeric(y) || is.logical(y)))
  stopifnot(ncol(x) == 2,ncol(y) == 2,all(!is.na(x)), all(!is.na(y)))
  stopifnot(all(abs(x[,2] - x[,1]) > .Machine$double.eps)) # Checks whether any dose difference are essentially zero.
  
  # Check on hypothesis testing and Gamma_bar
  stopifnot(is.numeric(Gamma_ave) && Gamma_ave> 0)
  alternative = match.arg(alternative,c("greater than","less than","two.sided"))
  method = match.arg(method,c("wilcoxon","dose.weighted","response.polynomial","U","dose.weighted.U"))
  
  
  
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
  q_I = q_I.fun(x=x,y=y,method=method,rank.order=rank.order,m=m,m1=m1,m2=m2,kappa=kappa)
  sense.variance = as.numeric(q_I^2%*%p_i_product)
  ## 
  test.stat = Test.stat.fun(x=x,y=y,method=method,rank.order=rank.order,m=m,m1=m1,m2=m2,kappa=kappa)
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
  return(list(pval = pval, sense.mean = sense.mean, sense.variance = sense.variance))
}


## x and y shall be obtained by simulating a population data 
design_sensitivity_fun = function(x,y,method = c("wilcoxon","dose.weighted","response.polynomial","U","dose.weighted.U"),rank.order=2,m=m,m1=m1,m2=m2,kappa=function(z) z){
  ## the number of matched pairs is I
  I = length(x)*0.5
  if(method=="wilcoxon"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    ## first obtain the empirical based cumulative distribution function F_y
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(F_y(abs(response.diff))*dose.response.product.sign)
    f = function(x){
      mean(F_y(abs(response.diff))*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation
    }
    gamma_star = uniroot(f,lower = 0.000001,upper=200,tol=.Machine$double.eps^0.9)$root
    
  }
  
  if(method=="dose.weighted"){
    dose.diff = x[,2] - x[,1]
    kappa.dose.diff = kappa(x[,2]) - kappa(x[,1])
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    ## first obtain the empirical based cumulative distribution function F_y
    F_y = ecdf(abs(response.diff))
    F_z = ecdf(abs(kappa.dose.diff))
    right.hand.side.expectation = mean(F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*dose.response.product.sign)
    f = function(x){
      mean(F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f,lower = 0.000001,upper=200,tol=.Machine$double.eps^0.9)$root
    
    
  }
  
  if(method=="response.polynomial"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    ## first obtain the empirical based cumulative distribution function F_y
    F_z = ecdf(abs(dose.diff))
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(F_y(abs(response.diff))^rank.order*dose.response.product.sign)
    f = function(x){
      mean(F_y(abs(response.diff))^rank.order*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f,lower = 0.000001,upper=200,tol=.Machine$double.eps^0.9)$root
    
  }
  
  if(method=="U"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    ## first obtain the empirical based cumulative distribution function F_y
    F_z = ecdf(abs(dose.diff))
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(multrnks(rk= F_y(abs(response.diff))*I,m=m,m1=m1,m2=m2)*dose.response.product.sign)
    f = function(x){
      mean((multrnks(rk=F_y(abs(response.diff))*I,m=m,m1=m1,m2=m2)*exp(x*abs(dose.diff)))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f,lower = 0.000001,upper=200,tol=.Machine$double.eps^0.9)$root
    
  }
  
  if(method=="dose.weighted.U"){
    dose.diff = x[,2] - x[,1]
    kappa.dose.diff = kappa(x[,2]) - kappa(x[,1])
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    ## first obtain the empirical based cumulative distribution function F_y
    F_z = ecdf(abs(kappa.dose.diff))
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(multrnks(rk= F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*I,m=m,m1=m1,m2=m2)*dose.response.product.sign)
    f = function(x){
      mean((multrnks(rk=F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*I,m=m,m1=m1,m2=m2)*exp(x*abs(dose.diff)))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f,lower = 0.000001,upper=200,tol=.Machine$double.eps^0.9)$root
    
  }
  
  
  ## and obtain the design sensitivity 
  Gamma_star = mean(exp(gamma_star*abs(x[,2]-x[,1])))
  return(Gamma_star)
  
  
}



Bahadur_slope_fun = function(x,y,Gamma_ave,method=c("wilcoxon","dose.weighted","U","response.polynomial","dose.weighted.U"),rank.order=2,m=2,m1=2,m2=2, upper_bound = 10, kappa=function(z) z){
  tilde_gamma = Gamma_to_tilde_gamma_fun(x=x, generalized.Gamma = Gamma_ave)
  Gamma_i = exp(tilde_gamma*abs(x[,1]-x[,2]))
  ## from Gamma_i to p_i 
  pi_plus  = Gamma_i/(1+Gamma_i)
  response.diff = y[,1] - y[,2]
  dose.diff = x[,1] - x[,2]
  kappa.dose.diff = kappa(x[,1]) - kappa(x[,2])
  diff.product = response.diff*dose.diff
  F_y = ecdf(abs(response.diff))
  F_z = ecdf(abs(kappa.dose.diff))
  I = length(x)/2
  
  
  if(method=="wilcoxon"){
    psi_star = F_y(abs(response.diff))
    
  }
  if(method=="dose.weighted"){
    psi_star = F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))
  }
  
  
  if(method=="U"){
    psi_star = multrnks(rk=F_y(abs(response.diff))*I,m=m,m1=m1,m2=m2) 
  }
  
  
  if(method=="response.polynomial"){
    psi_star = F_y(abs(response.diff))^rank.order
  }
  
  if(method=="dose.weighted.U"){
    psi_star = multrnks(rk=F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*I, m=m,m1=m1,m2=m2)
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
  Bahadur = tilde_t*mu - mean(log(pi_plus*exp(tilde_t*psi_star)+(1-pi_plus)))
  return(Bahadur)
  
}



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
    ## Worst-case correlation rho* of the adaptive procedure, equation (7).
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






generate_linear_effect = function(I, linear_coefficient,seed=100){
  ## set.seed 
  set.seed(seed)
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = runif(n=1,min=0.1,max=1)
    ## make sure that the high.dose is always higher than low dose 
    dose.diff =runif(n=1,min=0.1,max=1)
    high.dose = low.dose + dose.diff
    low.dose.response =  linear_coefficient*low.dose  + rnorm(n=1,mean=0,sd=1) - 1.2
    ## the high.dose.response is generated with the low dose response for an individual plus a constant plus a random error
    
    high.dose.response = linear_coefficient*high.dose + rnorm(n=1,mean=0,sd=1) - 1.2
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}


generate_kink_effect = function(I, linear_coefficient=1,seed=100, kink = 0.5){
  ## set.seed 
  set.seed(seed)
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = runif(n=1,min=0.1,max =1)
    ## make sure that the high.dose is always higher than low dose 
    dose.diff =runif(n=1,min=0.1,max=1)
    high.dose = low.dose + dose.diff
    ## generate the low.dose.response by a normal 
    low.dose.response =  linear_coefficient*ifelse(low.dose>kink,1,0)*(low.dose-kink)  + rnorm(n=1,mean=0,sd=1)
    ## the high.dose.response is generated with the low dose response for an individual plus a constant plus a random error
    
    high.dose.response = linear_coefficient*ifelse(high.dose>=kink,1,0)*(high.dose-kink) + rnorm(n=1,mean=0,1)
    
    
    
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}


generate_polynomial_effect = function(I, power,seed=100){
  ## set.seed 
  set.seed(seed)
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = runif(n=1,min=0.1,max=1)
    ## make sure that the high.dose is always higher than low dose 
    dose.diff =runif(n=1,min=0.1,max=1)
    high.dose = low.dose + dose.diff
    low.dose.response =  0.5*low.dose^power  + rnorm(n=1,mean=0,sd=1)
    ## the high.dose.response is generated with the low dose response for an individual plus a constant plus a random error
    
    high.dose.response = 0.5*high.dose^power + rnorm(n=1,mean=0,sd=1)
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}



generate_log_effect = function(I,seed=100){
  ## set.seed 
  set.seed(seed)
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = runif(n=1,min=0.1,max=1)
    ## make sure that the high.dose is always higher than low dose 
    dose.diff =runif(n=1,min=0.1,max=1)
    high.dose = low.dose + dose.diff
    low.dose.response =  0.75*log(low.dose)  + rnorm(n=1,mean=0,sd=1)-0.8
    ## the high.dose.response is generated with the low dose response for an individual plus a constant plus a random error
    
    high.dose.response = 0.75*log(high.dose) + rnorm(n=1,mean=0,sd=1)-0.8
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}



generate_sqrt_effect = function(I,seed=100){
  ## set.seed 
  set.seed(seed)
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = runif(n=1,min=0.1,max=1)
    ## make sure that the high.dose is always higher than low dose 
    dose.diff =runif(n=1,min=0.1,max=1)
    high.dose = low.dose + dose.diff
    low.dose.response =  1.6*sqrt(low.dose)  + rnorm(n=1,mean=0,sd=1)-1.8
    ## the high.dose.response is generated with the low dose response for an individual plus a constant plus a random error
    
    high.dose.response = 1.6*sqrt(high.dose) + rnorm(n=1,mean=0,sd=1)-1.8
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}





generate_flat_effect = function(I,seed=100){
  ## set.seed 
  set.seed(seed)
  f <- function(x) {
    a <- 1.2
    b <- 1.2
    y <- ifelse(x <= b, a * x, a * b)
    return(y)
  }
  ## first generate the low dose
  dose = matrix(data=NA,nrow = I,ncol=2)
  response = matrix(data=NA,nrow = I, ncol=2)
  for(i in 1:I){
    low.dose = runif(n=1,min=0.1,max=1)
    ## make sure that the high.dose is always higher than low dose 
    dose.diff =runif(n=1,min=0.1,max=1)
    high.dose = low.dose + dose.diff
    low.dose.response =  f(low.dose)  + rnorm(n=1,mean=0,sd=1)-1.5
    ## the high.dose.response is generated with the low dose response for an individual plus a constant plus a random error
    
    high.dose.response = f(high.dose) + rnorm(n=1,mean=0,sd=1)-1.5
    dose[i,1] = low.dose
    dose[i,2] = high.dose
    response[i,1] = low.dose.response
    response[i,2] = high.dose.response
  }
  
  return(list(dose=dose,response=response))
}




