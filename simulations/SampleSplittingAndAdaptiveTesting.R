################################################################################################################################
#       The Design sensitivities and the simulated power in finite sample settings,                                            #
#       compares the adaptive test, and two sample splitting using either 10% or 20% as the planning sample                    #
#       and the rest as the testing sample                                                                                     #
#       Table S.3, Table S.4, Table S.5, and Table S.6 can be generated using the following code skeleton                      #
#       by setting different I (matched pair numbers for computing the simulated powers)                                       #
#                                                                                                                              #
################################################################################################################################

# Table S.3, Table S.4, Table S.5, Table S.6 runner
#
# For one data-generating model, one pair of candidate tests, and one sample
# size I, this reproduces a sub-table in the following layout:
#
#                 Test1   Test2   Adaptive   Sample split
#   Gamma*        <ds1>   <ds2>   max(ds1,ds2)   max(ds1,ds2)
#   1.50          power   power   power          power
#   2.00          ...
#   ...
#   4.00          ...
#
# Row "Gamma*" is the generalized design sensitivity of each test, computed
# from a large favorable population via design_sensitivity_fun. The adaptive
# and sample-splitting columns both take the larger of the two individual
# design sensitivities based on the theorems in our paper and the one from 
# "Heller, R., Rosenbaum, P. R., & Small, D. S. (2009). 
# Split Samples and Design Sensitivity in Observational Studies. Journal of the American Statistical Association". 

# The remaining rows are simulated power at each Gamma
# under the model at the finite I (100, 500, 1000, and 5000)
#
# run_model_table() places the Wilcox/D-Wilcox sub-table and the U/D-U
# sub-table side by side, matching one model block of the printed table.
#
# Core routines (Gamma_to_tilde_gamma_fun, q_I.fun, Test.stat.fun, sendose,
# multrnks) are copied-pasted from matched_sensitivity_analysis_functions.R.
# Here we present them with the sample splitting procedure for comparison
###############################################################################

source("matched_sensitivity_analysis_functions.R")
require(mvtnorm)


###############################################################################
# 1. Design sensitivity (verbatim from Continuous_Outcome_Generalized_Design_
#    Bahadur_Efficiency.R). x, y are a large favorable population sample.
###############################################################################

design_sensitivity_fun = function(x, y, method = c("wilcoxon","dose.weighted","response.polynomial","U","dose.weighted.U"), rank.order=2, m=m, m1=m1, m2=m2, kappa=function(z) z){
  I = length(x)*0.5
  if(method=="wilcoxon"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(F_y(abs(response.diff))*dose.response.product.sign)
    f = function(x){
      mean(F_y(abs(response.diff))*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation
    }
    gamma_star = uniroot(f, lower = 0.000001, upper = 200, tol = .Machine$double.eps^0.9)$root
  }
  if(method=="dose.weighted"){
    dose.diff = x[,2] - x[,1]
    kappa.dose.diff = kappa(x[,2]) - kappa(x[,1])
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    F_y = ecdf(abs(response.diff))
    F_z = ecdf(abs(kappa.dose.diff))
    right.hand.side.expectation = mean(F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*dose.response.product.sign)
    f = function(x){
      mean(F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f, lower = 0.000001, upper = 200, tol = .Machine$double.eps^0.9)$root
  }
  if(method=="response.polynomial"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    F_z = ecdf(abs(dose.diff))
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(F_y(abs(response.diff))^rank.order*dose.response.product.sign)
    f = function(x){
      mean(F_y(abs(response.diff))^rank.order*exp(x*abs(dose.diff))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f, lower = 0.000001, upper = 200, tol = .Machine$double.eps^0.9)$root
  }
  if(method=="U"){
    dose.diff = x[,2] - x[,1]
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    F_z = ecdf(abs(dose.diff))
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(multrnks(rk = F_y(abs(response.diff))*I, m=m, m1=m1, m2=m2)*dose.response.product.sign)
    f = function(x){
      mean((multrnks(rk = F_y(abs(response.diff))*I, m=m, m1=m1, m2=m2)*exp(x*abs(dose.diff)))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f, lower = 0.000001, upper = 200, tol = .Machine$double.eps^0.9)$root
  }
  if(method=="dose.weighted.U"){
    dose.diff = x[,2] - x[,1]
    kappa.dose.diff = kappa(x[,2]) - kappa(x[,1])
    response.diff = y[,2] - y[,1]
    dose.response.product.sign = ifelse(dose.diff*response.diff>0,1,0)
    F_z = ecdf(abs(kappa.dose.diff))
    F_y = ecdf(abs(response.diff))
    right.hand.side.expectation = mean(multrnks(rk = F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*I, m=m, m1=m1, m2=m2)*dose.response.product.sign)
    f = function(x){
      mean((multrnks(rk = F_y(abs(response.diff))*F_z(abs(kappa.dose.diff))*I, m=m, m1=m1, m2=m2)*exp(x*abs(dose.diff)))/(1+exp(x*abs(dose.diff))))-right.hand.side.expectation}
    gamma_star = uniroot(f, lower = 0.000001, upper = 200, tol = .Machine$double.eps^0.9)$root
  }
  Gamma_star = mean(exp(gamma_star*abs(x[,2]-x[,1])))
  return(Gamma_star)
}


###############################################################################
# 2. Adaptive test and sample-splitting helpers
###############################################################################

find_sensitivity_value <- function(x, y,
                                   method = c("U", "wilcoxon", "dose.weighted",
                                              "response.polynomial", "dose.weighted.U"),
                                   m, m1, m2,
                                   alternative = "less than", rank.order = 3,
                                   target = 0.05,
                                   lower = 0.01, upper = 10, init = 1.0) {
  objective_function <- function(Gamma_ave) {
    result <- sendose(x = x, y = y, method = method, m = m, m1 = m1,
                      m2 = m2, alternative = alternative, Gamma_ave = Gamma_ave,
                      rank.order = rank.order)$pval
    return((result - target)^2)
  }
  opt <- optim(par = init, fn = objective_function, method = "Brent",
               lower = lower, upper = upper)
  Gamma_ave <- opt$par
  tilde.gamma <- Gamma_to_tilde_gamma_fun(x = x, generalized.Gamma = Gamma_ave)
  abs.dose.diff <- abs(x[, 1] - x[, 2])
  list(Gamma_ave = Gamma_ave,
       smallest_Gamma_i = exp(tilde.gamma * min(abs.dose.diff)),
       largest_Gamma_i = exp(tilde.gamma * max(abs.dose.diff)))
}


sample_split_procedure <- function(
    x, y,
    method1, method.1.rank.order = 2, method.1.m = 2, method.1.m1 = 2, method.1.m2 = 2,
    method2, method.2.rank.order = 2, method.2.m = 2, method.2.m1 = 2, method.2.m2 = 2,
    Gamma_ave, alpha = 0.05, pilot.prop = 0.10,
    alternative = "less than", split.seed = NULL) {
  
  stopifnot(is.matrix(x), is.matrix(y), nrow(x) == nrow(y))
  stopifnot(all(abs(x[, 2] - x[, 1]) > .Machine$double.eps))
  stopifnot(pilot.prop > 0 && pilot.prop < 1)
  stopifnot(is.numeric(Gamma_ave) && Gamma_ave > 0)
  alternative <- match.arg(alternative, c("greater than", "less than"))
  
  I <- nrow(x)
  if (!is.null(split.seed)) set.seed(split.seed)
  
  n.pilot <- max(1, floor(pilot.prop * I))
  pilot.idx <- sample(1:I, size = n.pilot, replace = FALSE)
  analysis.idx <- setdiff(1:I, pilot.idx)
  
  x.pilot <- x[pilot.idx, , drop = FALSE]; y.pilot <- y[pilot.idx, , drop = FALSE]
  x.analysis <- x[analysis.idx, , drop = FALSE]; y.analysis <- y[analysis.idx, , drop = FALSE]
  
  sens1 <- find_sensitivity_value(x.pilot, y.pilot, method = method1,
                                  rank.order = method.1.rank.order, m = method.1.m,
                                  m1 = method.1.m1, m2 = method.1.m2,
                                  alternative = alternative, target = alpha)
  sens2 <- find_sensitivity_value(x.pilot, y.pilot, method = method2,
                                  rank.order = method.2.rank.order, m = method.2.m,
                                  m1 = method.2.m1, m2 = method.2.m2,
                                  alternative = alternative, target = alpha)
  
  if (sens1$Gamma_ave >= sens2$Gamma_ave) {
    sel <- list(method = method1, rank.order = method.1.rank.order,
                m = method.1.m, m1 = method.1.m1, m2 = method.1.m2,
                sensitivity.value = sens1$Gamma_ave)
  } else {
    sel <- list(method = method2, rank.order = method.2.rank.order,
                m = method.2.m, m1 = method.2.m1, m2 = method.2.m2,
                sensitivity.value = sens2$Gamma_ave)
  }
  
  final.res <- sendose(x = x.analysis, y = y.analysis, method = sel$method,
                       rank.order = sel$rank.order, m = sel$m, m1 = sel$m1, m2 = sel$m2,
                       alternative = alternative, Gamma_ave = Gamma_ave)
  
  list(selected.method = sel$method,
       selected.sensitivity.value = sel$sensitivity.value,
       p.value = final.res$pval,
       reject = as.integer(final.res$pval < alpha))
}


sendose.adaptive <- function(
    x, y,
    method1, method.1.rank.order = 2, method.1.m = 2, method.1.m1 = 2, method.1.m2 = 2,
    method2, method.2.rank.order = 2, method.2.m = 2, method.2.m1 = 2, method.2.m2 = 2,
    log.p = FALSE, alternative = "greater than", Gamma_ave = 1, alpha = 0.05) {
  
  stopifnot(all(abs(x[, 2] - x[, 1]) > .Machine$double.eps))
  stopifnot(is.numeric(Gamma_ave) && Gamma_ave > 0)
  alternative <- match.arg(alternative, c("greater than", "less than"))
  
  find_Q_given_correlation <- function(alpha, correlation) {
    finderFS <- function(Q, correlation, alpha)
      1 - mvtnorm::pmvnorm(upper = rep(Q, 2), corr = matrix(c(1, correlation, correlation, 1), 2))[1] - alpha
    uniroot(finderFS, c(0, qnorm(1 - alpha / 4)), correlation, alpha)$root
  }
  
  find.worst.case.cor.obj.p <- function(dose.score.mat, tilde.gamma) {
    ## Worst-case correlation rho* of the adaptive procedure, equation (7).
    ## Under the Rosenbaum sensitivity bound (3), the least favourable within-pair
    ## treatment-assignment probability is attained at the upper endpoint
    ##     p_plus_i = Gamma_i / (1 + Gamma_i),   Gamma_i = exp(tilde.gamma * |D_i|).
    ## rho* is therefore evaluated in closed form at this single p_plus; no
    ## optimization is performed. The same p_plus underlies the marginal
    ## standardizations returned by sendose(), so every component of the procedure
    ## is evaluated at one common worst-case configuration, as in equation (7).
    cor.matched.pair.obj.p <- function(dose.score.mat, pi) {
      q_i <- dose.score.mat[, 2]; s_i <- dose.score.mat[, 3]
      cov <- sum(q_i * s_i * pi * (1 - pi))
      sd1 <- sqrt(sum(q_i^2 * pi * (1 - pi))); sd2 <- sqrt(sum(s_i^2 * pi * (1 - pi)))
      cov / sd1 / sd2
    }
    abs.dose.diff <- dose.score.mat[, 1]
    Gamma_i <- exp(tilde.gamma * abs.dose.diff)
    p.plus  <- Gamma_i / (1 + Gamma_i)                 ## upper Rosenbaum bound, eq. (7)
    list(worst.case.pi          = p.plus,
         worst.case.correlation = cor.matched.pair.obj.p(dose.score.mat, p.plus))
  }
  
  abs.dose.diff <- abs(x[, 1] - x[, 2])
  q_i <- q_I.fun(x, y, method1, method.1.rank.order, method.1.m, method.1.m1, method.1.m2)
  s_i <- q_I.fun(x, y, method2, method.2.rank.order, method.2.m, method.2.m1, method.2.m2)
  tilde.gamma <- Gamma_to_tilde_gamma_fun(x = x, generalized.Gamma = Gamma_ave)
  
  dose.score.mat <- cbind(abs.dose.diff, q_i, s_i)
  worst.case.res <- find.worst.case.cor.obj.p(dose.score.mat, tilde.gamma)
  worst.case.correlation <- worst.case.res$worst.case.correlation
  
  test.stat.1 <- Test.stat.fun(x, y, method1, method.1.rank.order, method.1.m, method.1.m1, method.1.m2)
  test.stat.2 <- Test.stat.fun(x, y, method2, method.2.rank.order, method.2.m, method.2.m1, method.2.m2)
  
  critical.value <- find_Q_given_correlation(alpha = alpha, correlation = worst.case.correlation)
  
  test.1.res <- sendose(x = x, y = y, method = method1, rank.order = method.1.rank.order,
                        m = method.1.m, m1 = method.1.m1, m2 = method.1.m2,
                        alternative = alternative, Gamma_ave = Gamma_ave)
  test.2.res <- sendose(x = x, y = y, method = method2, rank.order = method.2.rank.order,
                        m = method.2.m, m1 = method.2.m1, m2 = method.2.m2,
                        alternative = alternative, Gamma_ave = Gamma_ave)
  
  std.1 <- (test.stat.1 - test.1.res$sense.mean) / sqrt(test.1.res$sense.variance)
  std.2 <- (test.stat.2 - test.2.res$sense.mean) / sqrt(test.2.res$sense.variance)
  
  if (alternative == "greater than") {
    rejection <- ifelse(std.1 >= critical.value || std.2 >= critical.value, 1, 0)
  } else {
    rejection <- ifelse(std.1 <= -critical.value || std.2 <= -critical.value, 1, 0)
  }
  
  list(rejection = rejection, critical.value = critical.value,
       worst.case.correlation = worst.case.correlation,
       std.1 = std.1, std.2 = std.2, alternative = alternative)
}


###############################################################################
# 3. Data-generating processes and model definitions
###############################################################################

generate_linear_effect <- function(I, linear_coefficient, seed = 100) {
  set.seed(seed)
  dose <- matrix(NA, I, 2); response <- matrix(NA, I, 2)
  for (i in 1:I) {
    low.dose <- runif(1, 0.1, 1); high.dose <- low.dose + runif(1, 0.1, 1)
    dose[i, ] <- c(low.dose, high.dose)
    response[i, 1] <- linear_coefficient * low.dose + rnorm(1) - 1.2
    response[i, 2] <- linear_coefficient * high.dose + rnorm(1) - 1.2
  }
  list(dose = dose, response = response)
}

generate_kink_effect <- function(I, linear_coefficient = 1, seed = 100, kink = 0.5) {
  set.seed(seed)
  dose <- matrix(NA, I, 2); response <- matrix(NA, I, 2)
  for (i in 1:I) {
    low.dose <- runif(1, 0.1, 1); high.dose <- low.dose + runif(1, 0.1, 1)
    dose[i, ] <- c(low.dose, high.dose)
    response[i, 1] <- linear_coefficient * ifelse(low.dose > kink, 1, 0) * (low.dose - kink) + rnorm(1)
    response[i, 2] <- linear_coefficient * ifelse(high.dose >= kink, 1, 0) * (high.dose - kink) + rnorm(n = 1, mean = 0, 1)
  }
  list(dose = dose, response = response)
}

generate_polynomial_effect <- function(I, power, seed = 100) {
  set.seed(seed)
  dose <- matrix(NA, I, 2); response <- matrix(NA, I, 2)
  for (i in 1:I) {
    low.dose <- runif(1, 0.1, 1); high.dose <- low.dose + runif(1, 0.1, 1)
    dose[i, ] <- c(low.dose, high.dose)
    response[i, 1] <- 0.5 * low.dose^power + rnorm(1)
    response[i, 2] <- 0.5 * high.dose^power + rnorm(1)
  }
  list(dose = dose, response = response)
}

generate_log_effect <- function(I, seed = 100) {
  set.seed(seed)
  dose <- matrix(NA, I, 2); response <- matrix(NA, I, 2)
  for (i in 1:I) {
    low.dose <- runif(1, 0.1, 1); high.dose <- low.dose + runif(1, 0.1, 1)
    dose[i, ] <- c(low.dose, high.dose)
    response[i, 1] <- 0.75 * log(low.dose) + rnorm(1) - 0.8
    response[i, 2] <- 0.75 * log(high.dose) + rnorm(1) - 0.8
  }
  list(dose = dose, response = response)
}

generate_sqrt_effect <- function(I, seed = 100) {
  set.seed(seed)
  dose <- matrix(NA, I, 2); response <- matrix(NA, I, 2)
  for (i in 1:I) {
    low.dose <- runif(1, 0.1, 1); high.dose <- low.dose + runif(1, 0.1, 1)
    dose[i, ] <- c(low.dose, high.dose)
    response[i, 1] <- 1.6 * sqrt(low.dose) + rnorm(1) - 1.8
    response[i, 2] <- 1.6 * sqrt(high.dose) + rnorm(1) - 1.8
  }
  list(dose = dose, response = response)
}

generate_flat_effect <- function(I, seed = 100) {
  set.seed(seed)
  f <- function(x) { a <- 1.2; b <- 1.2; ifelse(x <= b, a * x, a * b) }
  dose <- matrix(NA, I, 2); response <- matrix(NA, I, 2)
  for (i in 1:I) {
    low.dose <- runif(1, 0.1, 1); high.dose <- low.dose + runif(1, 0.1, 1)
    dose[i, ] <- c(low.dose, high.dose)
    response[i, 1] <- f(low.dose) + rnorm(1) - 1.5
    response[i, 2] <- f(high.dose) + rnorm(1) - 1.5
  }
  list(dose = dose, response = response)
}

# Each model is function(I, seed) -> list(dose, response).
models <- list(
  square     = function(I, seed) generate_polynomial_effect(I = I, power = 2, seed = seed),
  kink       = function(I, seed) generate_kink_effect(I = I, linear_coefficient = 1.5, kink = 0.8, seed = seed),
  linear     = function(I, seed) generate_linear_effect(I = I, linear_coefficient = 1, seed = seed),
  sqrt       = function(I, seed) generate_sqrt_effect(I = I, seed = seed),
  flat       = function(I, seed) generate_flat_effect(I = I, seed = seed),
  log        = function(I, seed) generate_log_effect(I = I, seed = seed)
)


###############################################################################
# 4. Sub-table runner
#
# Produces one Table S.3--S.6 style sub-table: a Gamma* row (design sensitivities, with
# the adaptive and sample-splitting columns both equal to the larger of the
# two individual design sensitivities) and power rows for the two individual
# tests, the adaptive test, and sample splitting. The single-test rejections
# are read from the same standardized deviates computed inside the adaptive
# test, so no separate single-test call is needed.
###############################################################################

run_subtable <- function(model, I, method1, method2,
                         label1 = method1, label2 = method2,
                         m1.set = list(rank.order = 2, m = 2, m1 = 2, m2 = 2),
                         m2.set = list(rank.order = 2, m = 2, m1 = 2, m2 = 2),
                         Gamma_vec = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
                         n.runs = 1000, alpha = 0.01,
                         alternative = "greater than", pilot.prop = 0.1,
                         design = NULL, pop.I = 5e5, pop.seed = 100,
                         verbose = TRUE) {
  
  ## ---- design sensitivities (computed once; independent of finite I) ----
  if (is.null(design)) {
    if (verbose) cat("Computing design sensitivities from population (pop.I =", pop.I, ") ...\n")
    pop <- model(pop.I, pop.seed)
    ds1 <- design_sensitivity_fun(pop$dose, pop$response, method = method1,
                                  rank.order = m1.set$rank.order, m = m1.set$m,
                                  m1 = m1.set$m1, m2 = m1.set$m2)
    ds2 <- design_sensitivity_fun(pop$dose, pop$response, method = method2,
                                  rank.order = m2.set$rank.order, m = m2.set$m,
                                  m1 = m2.set$m1, m2 = m2.set$m2)
  } else {
    ds1 <- design[1]; ds2 <- design[2]
  }
  ds.adaptive <- max(ds1, ds2)
  
  z_alpha <- qnorm(1 - alpha)
  nG <- length(Gamma_vec)
  p1 <- p2 <- pad <- psp <- rep(NA, nG)
  
  for (g in seq_len(nG)) {
    Gamma_ave <- Gamma_vec[g]
    r1 <- r2 <- rad <- rsp <- integer(n.runs)
    t0 <- Sys.time()
    
    for (run in seq_len(n.runs)) {
      d <- model(I, run); x <- d$dose; y <- d$response
      
      ad <- sendose.adaptive(
        x = x, y = y, method1 = method1, method2 = method2,
        method.1.rank.order = m1.set$rank.order, method.1.m = m1.set$m, method.1.m1 = m1.set$m1, method.1.m2 = m1.set$m2,
        method.2.rank.order = m2.set$rank.order, method.2.m = m2.set$m, method.2.m1 = m2.set$m1, method.2.m2 = m2.set$m2,
        alternative = alternative, Gamma_ave = Gamma_ave, alpha = alpha
      )
      if (alternative == "greater than") {
        r1[run] <- as.integer(ad$std.1 >= z_alpha)
        r2[run] <- as.integer(ad$std.2 >= z_alpha)
      } else {
        r1[run] <- as.integer(ad$std.1 <= -z_alpha)
        r2[run] <- as.integer(ad$std.2 <= -z_alpha)
      }
      rad[run] <- ad$rejection
      
      sp <- sample_split_procedure(
        x = x, y = y, method1 = method1, method2 = method2,
        method.1.rank.order = m1.set$rank.order, method.1.m = m1.set$m, method.1.m1 = m1.set$m1, method.1.m2 = m1.set$m2,
        method.2.rank.order = m2.set$rank.order, method.2.m = m2.set$m, method.2.m1 = m2.set$m1, method.2.m2 = m2.set$m2,
        Gamma_ave = Gamma_ave, alpha = alpha, pilot.prop = pilot.prop,
        alternative = alternative, split.seed = run
      )
      rsp[run] <- sp$reject
      
      if (verbose && run %% 100 == 0)
        cat("  Gamma", Gamma_ave, "run", run, "|",
            label1, round(mean(r1[1:run]), 2), label2, round(mean(r2[1:run]), 2),
            "Adapt", round(mean(rad[1:run]), 2), "Split", round(mean(rsp[1:run]), 2), "\n")
    }
    
    p1[g] <- mean(r1); p2[g] <- mean(r2); pad[g] <- mean(rad); psp[g] <- mean(rsp)
    if (verbose) cat("Gamma", Gamma_ave, "done in",
                     round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "secs\n")
  }
  
  tab <- rbind(c(ds1, ds2, ds.adaptive, ds.adaptive), cbind(p1, p2, pad, psp))
  rownames(tab) <- c("Gamma*", formatC(Gamma_vec, format = "f", digits = 2))
  colnames(tab) <- c(label1, label2, "Adaptive", "Sample split")
  round(tab, 2)
}


###############################################################################
# 5. Full model block (both test pairs side by side) and printing
###############################################################################

# Replace NA with blank for printing (matches the empty Sample split cell).
print_block <- function(mat) {
  out <- formatC(mat, format = "f", digits = 2)
  out[is.na(mat)] <- ""
  dim(out) <- dim(mat); dimnames(out) <- dimnames(mat)
  print(noquote(out))
  invisible(out)
}

run_model_table <- function(model, I, model.name = "",
                            Gamma_vec = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
                            n.runs = 1000, alpha = 0.01,
                            alternative = "greater than", pilot.prop = 0.1,
                            pop.I = 5e5, pop.seed = 100, verbose = TRUE) {
  
  left <- run_subtable(model, I, method1 = "wilcoxon", method2 = "dose.weighted",
                       label1 = "Wilcox", label2 = "D-Wilcox",
                       m1.set = list(rank.order = 2, m = 2, m1 = 2, m2 = 2),
                       m2.set = list(rank.order = 2, m = 2, m1 = 2, m2 = 2),
                       Gamma_vec = Gamma_vec, n.runs = n.runs, alpha = alpha,
                       alternative = alternative, pilot.prop = pilot.prop,
                       pop.I = pop.I, pop.seed = pop.seed, verbose = verbose)
  
  right <- run_subtable(model, I, method1 = "U", method2 = "dose.weighted.U",
                        label1 = "U", label2 = "D-U",
                        m1.set = list(rank.order = 2, m = 8, m1 = 7, m2 = 8),
                        m2.set = list(rank.order = 2, m = 8, m1 = 7, m2 = 8),
                        Gamma_vec = Gamma_vec, n.runs = n.runs, alpha = alpha,
                        alternative = alternative, pilot.prop = pilot.prop,
                        pop.I = pop.I, pop.seed = pop.seed, verbose = verbose)
  
  full <- cbind(left, right)
  cat("\n==== ", model.name, "  (I =", I, ", pilot.prop =", pilot.prop,
      ", alpha =", alpha, ") ====\n")
  print_block(full)
  invisible(list(left = left, right = right, full = full))
}


###############################################################################
# 6. Example: one model, one I
#
# This reproduces the Polynomial (Square) dose-response block of Table S.4 at I = 500.
# with sample splitting (0.10) column. The design
# sensitivity uses pop.I = 5e5 to compute the values based on our paper; 
# For simulating the power in finite sample setting, the following uses 
# I = 500, that is, to compute the empirical by counting the proportion of
# rejection with I = 500 matched pairs every time for 1000 draws 
# the default is pilot.prop = 0.10, corresponding to take 0.10 of data as planning samples
# One can then see in the following that we also specify pilot.prop = 0.20 to fully
# recover the simulations
###############################################################################

if (FALSE) {
  
  # full block (both sub-tables) for the square model at I = 500
  run_model_table(models$square, I = 500, model.name = "Polynomial")
  
  # a single sub-table only (Wilcox vs D-Wilcox, square, I = 500)
  run_subtable(models$square, I = 500,
               method1 = "wilcoxon", method2 = "dose.weighted",
               label1 = "Wilcox", label2 = "D-Wilcox")
  
}





#### Here starts the experiment ##############

#########################################################################################
#  Table S.3, I = 100, adaptive test with sample split (0.10) and sample split (0.20)
#
########################################################################################

### First run the model with square dose-response curve and I = 100, and sample splitting planning
### with planning sample 0.1 of the data, this 
### recovers Table S.3, the first four columns of square dose-response curve 
### with \kappa(z) = z, i.e., the columns of 
### (Wilcox, D-Wilcox, Adaptive, Split (0.10); U, D-U, Adaptive, Split (0.10)) 

run_model_table(models$square, I = 100, model.name = "Square")


### the model with kink and I = 100, and sample splitting planning sample 0.1

run_model_table(models$kink, I = 100, model.name = "Kink")


### the model with linear and I = 100, and sample splitting planning sample 0.1

run_model_table(models$linear, I = 100, model.name = "Linear", pilot.prop = 0.10)


## Square Root model I = 100, sample splitting planning sample 0.10

run_model_table(models$sqrt, I = 100, model.name = "Squared Root")


## Flat model I= 100, sample splitting planning sample 0.10

run_model_table(models$flat, I = 100, model.name = "Flat")


## Log model I = 100, sample splitting planning sample 0.10

run_model_table(models$log, I = 100, model.name = "Log")




### This batch runs the model with square and I = 100, and sample splitting planning
### with planning sample 0.20 of the data, this 
### recovers Table S.3, the first three columns of square dose-response curve 
### with \kappa(z) = z and the last column of a test X dose response curve combination, i.e.,
### the columns of 
### (Wilcox, D-Wilcox, Adaptive, Split (0.20); U, D-U, Adaptive, Split (0.20))

run_model_table(models$square, I = 100, model.name = "Square", pilot.prop = 0.20)


### the model with kink and I = 100, and sample splitting planning sample 0.20

run_model_table(models$kink, I = 100, model.name = "Kink", pilot.prop = 0.20)


### the model with linear and I = 100, and sample splitting planning sample 0.2

run_model_table(models$linear, I = 100, model.name = "Linear", pilot.prop = 0.20)


## Square Root model I = 100, sample splitting planning sample 0.20

run_model_table(models$sqrt, I = 100, model.name = "Squared Root", pilot.prop = 0.20)


## Flat model I= 100, sample splitting planning sample 0.20

run_model_table(models$flat, I = 100, model.name = "Flat", pilot.prop = 0.20)


## Log model I = 100, sample splitting planning sample 0.20

run_model_table(models$log, I = 100, model.name = "Log", pilot.prop = 0.20)






#########################################################################################
#  Table S.4, I = 500, adaptive test with sample split (0.10) and sample split (0.20)
#
########################################################################################


### First run the model with square and I = 500, and sample splitting planning
### sample 0.1, this batch will recover Table S.4, the first four column of each 
### Dose-response curve X test, i.e., the columns of 
### (Wilcox, D-Wilcox, Adaptive, Split (0.10); U, D-U, Adaptive, Split (0.10))

run_model_table(models$square, I = 500, model.name = "Square")


### the model with kink and I = 500, and sample splitting planning sample 0.1

run_model_table(models$kink, I = 500, model.name = "Kink")


### the model with linear and I = 500, and sample splitting planning sample 0.1

run_model_table(models$linear, I = 500, model.name = "Linear")


## Squared Root model I = 500, sample splitting planning sample 0.10

run_model_table(models$sqrt, I = 500, model.name = "Squared Root")


## Flat model I= 500, sample splitting planning sample 0.10

run_model_table(models$flat, I = 500, model.name = "Flat")


## Log model I = 500, sample splitting planning sample 0.10

run_model_table(models$log, I = 500, model.name = "Log")






### The model with square and I = 500, and sample splitting planning
### sample 0.2, this batch will recover Table S.4, the columns of 
### (Wilcox, D-Wilcox, Adaptive, Split (0.20); U, D-U, Adaptive, Split (0.20))

run_model_table(models$square, I = 500, model.name = "Square", pilot.prop = 0.20)


### the model with kink and I = 500, and sample splitting planning sample 0.2

run_model_table(models$kink, I = 500, model.name = "Kink", pilot.prop = 0.20)


### the model with linear and I = 500, and sample splitting planning sample 0.20

run_model_table(models$linear, I = 500, model.name = "Linear", pilot.prop = 0.20)


## Squared Root model I = 500, sample splitting planning sample 0.20

run_model_table(models$sqrt, I = 500, model.name = "Squared Root", pilot.prop = 0.20)


## Flat model I= 500, sample splitting planning sample 0.20

run_model_table(models$flat, I = 500, model.name = "Flat", pilot.prop = 0.20)


## Log model I = 500, sample splitting planning sample 0.20

run_model_table(models$log, I = 500, model.name = "Log", pilot.prop = 0.20)





#########################################################################################
#  Table S.5, I = 1000, adaptive test with sample split (0.10) and sample split (0.20)
#
########################################################################################



######### I = 1000, planning sample 0.10

run_model_table(models$square, I = 1000, model.name = "Square")



run_model_table(models$kink, I = 1000, model.name = "Kink")


run_model_table(models$linear, I = 1000, model.name = "Linear")

run_model_table(models$sqrt, I = 1000, model.name = "Squared Root")

run_model_table(models$flat, I = 1000, model.name = "Flat")

run_model_table(models$log, I = 1000, model.name = "Log")



######### I = 1000, planning sample 0.20

run_model_table(models$square, I = 1000, model.name = "Square", pilot.prop = 0.20)


run_model_table(models$kink, I = 1000, model.name = "Kink", pilot.prop = 0.20)


run_model_table(models$linear, I = 1000, model.name = "Linear", pilot.prop = 0.20)

run_model_table(models$sqrt, I = 1000, model.name = "Squared Root", pilot.prop = 0.20)

run_model_table(models$flat, I = 1000, model.name = "Flat", pilot.prop = 0.20)

run_model_table(models$log, I = 1000, model.name = "Log", pilot.prop = 0.20)






#########################################################################################
#  Table S.6, I = 5000, adaptive test with sample split (0.10) and sample split (0.20)
#
########################################################################################


######### I = 5000, planning sample 0.10

run_model_table(models$square, I = 5000, model.name = "Square")

run_model_table(models$kink, I = 5000, model.name = "Kink")

run_model_table(models$linear, I = 5000, model.name="Linear")

run_model_table(models$sqrt, I = 5000, model.name="Squared Root")

run_model_table(models$flat, I = 5000, model.name="Flat")

run_model_table(model = models$log, I = 5000, model.name = "Log")




######### I = 5000, planning sample 0.20

run_model_table(models$square, I = 5000, model.name = "Square", pilot.prop = 0.20)

run_model_table(models$kink, I = 5000, model.name = "Kink", pilot.prop = 0.20)

run_model_table(models$linear, I = 5000, model.name="Linear", pilot.prop = 0.20)

run_model_table(models$sqrt, I = 5000, model.name="Squared Root", pilot.prop = 0.20)

run_model_table(models$flat, I = 5000, model.name="Flat", pilot.prop = 0.20)

run_model_table(model = models$log, I = 5000, model.name = "Log", pilot.prop = 0.20)







