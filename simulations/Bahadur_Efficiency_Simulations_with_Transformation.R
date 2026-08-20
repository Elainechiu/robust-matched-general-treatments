#############################################################################################################################################
#   The Bahadur Rosenbaum Relative efficiency simulations for Table S.9, Table S.10, Table S.11 and Table S.12                              #
#   considering different dose transformation for test design in Section C.3 in the supplementary material                                  #
#   Another batch of Bahadur rosenbaum relative efficiency simulations is in "Continuous_Outcome_Generalized_Design_Bahadur_Efficiency.R"   #
#   which is equivalent to setting \kappa(z) = z based on the dose-transformation theorem in Section C.3                                    #
#   This file will outputs in csv format due to a higher computational time compared to other simulations                                   #
#   See "bahadur_relative_efficiency_verification_long.csv" for the reproduction of Table S.9 in the supplementary material,                #
#   or run the following "RUN_PART_A" coding part and call "releff_long" to see the results on screen;                                      #
#                                                                                                                                           #
#   See "bahadur_slopes_by_kappa_long.csv" for the reproduction of Table S.10, S.11, and S.12 in the supplementary material,                #
#   or run the following "RUN_PART_B" coding part and call "slope_long" to see the results on screen.                                       #
#                                                                                                                                           #
#############################################################################################################################################


source("matched_sensitivity_analysis_functions.R")

## ============================================================================
## Configuration.
## ============================================================================
ALPHA      <- 0.01          # significance level 
POWER      <- 0.95          # target power for the minimal-sample-size search
MC_ITER    <- 1000          # Monte Carlo replications per candidate I
GAMMA_GRID <- c(1.00, 1.10, 1.30, 1.50)
NPOP       <- 5e5           # population size for the analytic slope and
# design-sensitivity formulas, this setup is the same as other simulations

## U-statistic parameters (m, m_lower, m_upper) = (8, 7, 8), as in the main text.
U_m  <- 8; U_m1 <- 7; U_m2 <- 8

CURVES_TO_RUN <- c("Square", "Kink", "Linear", "Square Root", "Flat", "Log")
KAPPAS_TO_RUN <- c("z", "log(z)", "z^2", "sqrt(z)")

RUN_PART_B <- TRUE          # analytic slopes
RUN_PART_A <- TRUE          # simulated verification

## ============================================================================
## Dose-response curve registry.
##   pop : a single large population (I = NPOP) used for the analytic slope and
##         design-sensitivity formulas.

##   gen : generator gen(I, seed) -> list(dose, response) used by the Monte
##         Carlo minimal-sample-size search, i.e., this part gives I number of 
##         matched pairs under a particular population model (polynomial,kink,...log)
##         and will return a data frame with dose and response, which is the data 
##         structure that the sendose, designsensitivity function rely on

##   ub  : initial upper bound on the number of matched pairs for the size
##         search with bisection method
##         ; extended automatically when the target power is not yet
##         bracketed (see min_pairs_fun).
## ============================================================================
curves <- list(
  "Square"  = list(pop = generate_polynomial_effect(I = NPOP, power = 2),
                       gen = function(I, seed) generate_polynomial_effect(I = I, power = 2, seed = seed),
                       ub  = 1000),
  "Kink"        = list(pop = generate_kink_effect(I = NPOP, linear_coefficient = 1.5, kink = 0.8),
                       gen = function(I, seed) generate_kink_effect(I = I, linear_coefficient = 1.5, kink = 0.8, seed = seed),
                       ub  = 1500),
  "Linear"      = list(pop = generate_linear_effect(I = NPOP, linear_coefficient = 1),
                       gen = function(I, seed) generate_linear_effect(I = I, linear_coefficient = 1, seed = seed),
                       ub  = 1000),
  "Square Root" = list(pop = generate_sqrt_effect(I = NPOP),
                       gen = function(I, seed) generate_sqrt_effect(I = I, seed = seed),
                       ub  = 1000),
  "Flat"        = list(pop = generate_flat_effect(I = NPOP),
                       gen = function(I, seed) generate_flat_effect(I = I, seed = seed),
                       ub  = 1000),
  "Log"         = list(pop = generate_log_effect(I = NPOP),
                       gen = function(I, seed) generate_log_effect(I = I, seed = seed),
                       ub  = 2000)
)

## Dose transformation for the test design, see Section C.3 for definition of 
## \kappa(z) and \phi(z)
##
## The Doses lie in (0.1, 2) from either one of the dose-response curves, 
## so log and sqrt are well defined.
kappas <- list(
  "z"       = function(z) z,
  "log(z)"  = function(z) log(z),
  "z^2"     = function(z) z^2,
  "sqrt(z)" = function(z) sqrt(z)
)

## Test registry: internal method name, display label, and whether the test
## uses the dose (and therefore depends on kappa).
tests <- list(
  list(method = "wilcoxon",        label = "Wilcox",   dose = FALSE),
  list(method = "dose.weighted",   label = "D-Wilcox", dose = TRUE),
  list(method = "U",               label = "U",        dose = FALSE),
  list(method = "dose.weighted.U", label = "D-U",      dose = TRUE)
)

## ============================================================================
## Core helpers.
## the "design_sensitivity_fun" and "Bahadur_slope_fun" are from 
## "matched_sensitivity_analysis_functions.R"
## ============================================================================

## Generalized design sensitivity Gamma_bar_* for a test under transformation
## kappa.
design_sens <- function(pop, method, kappa) {
  design_sensitivity_fun(x = pop$dose, y = pop$response, method = method,
                         m = U_m, m1 = U_m1, m2 = U_m2, kappa = kappa)
}

## Generalized Bahadur-Rosenbaum exact slope Upsilon, gated on
## Gamma_ave < Gamma_bar_*. Returns NA when the slope is not well defined.
slope_guarded <- function(pop, method, Gamma_ave, kappa, Gstar = NULL) {
  if (is.null(Gstar)) Gstar <- design_sens(pop, method, kappa)
  if (!is.finite(Gstar) || Gamma_ave >= Gstar) return(NA_real_)
  tryCatch(
    Bahadur_slope_fun(x = pop$dose, y = pop$response, Gamma_ave = Gamma_ave,
                      method = method, m = U_m, m1 = U_m1, m2 = U_m2, kappa = kappa),
    error = function(e) NA_real_, warning = function(w) NA_real_)
}

## Minimal number of matched pairs to attain POWER at ALPHA, by Monte Carlo.
## A uniroot(bisection) search over I, mirroring the existing sample_size_*_fun routines.
## The upper bracket is extended geometrically (up to max.upper) until the
## target power is bracketed; this is invoked only for Gamma_ave below the
## design sensitivity, where the power tends to one and a finite root exists.
## Returns NA if the target remains unreachable within max.upper.
min_pairs_fun <- function(gen, method, Gamma_ave, kappa = function(z) z,
                          upper.bound = 2000, lower = 10, max.upper = 8000,
                          rank.order = 2, m = U_m, m1 = U_m1, m2 = U_m2,
                          mc.iter = MC_ITER, tol = 1, verbose = FALSE,
                          tag = "") {
  f <- function(I) {
    I <- max(1L, round(I))
    pvals <- numeric(mc.iter)
    for (mc in seq_len(mc.iter)) {
      d <- gen(I = I, seed = mc)
      pvals[mc] <- sendose(x = d$dose, y = d$response, method = method,
                           alternative = "greater than", Gamma_ave = Gamma_ave,
                           rank.order = rank.order, m = m, m1 = m1, m2 = m2,
                           kappa = kappa)$pval
    }
    pw <- mean(pvals < ALPHA)
    if (verbose) {
      cat(sprintf("        %s try I = %5d  ->  power = %.3f (target %.2f)\n",
                  tag, I, pw, POWER)); flush.console()
    }
    pw - POWER
  }
  
  ub      <- upper.bound
  f.lower <- f(lower)
  if (f.lower >= 0) return(ceiling(lower))
  f.upper <- f(ub)
  while (f.upper < 0 && ub < max.upper) {
    ub      <- min(max.upper, ub * 2)
    f.upper <- f(ub)
  }
  if (f.upper < 0) return(NA_integer_)
  
  root <- tryCatch(
    uniroot(f, lower = lower, upper = ub, f.lower = f.lower, f.upper = f.upper,
            tol = tol)$root,
    error = function(e) NA_real_)
  if (is.na(root)) NA_integer_ else ceiling(root)
}

## ============================================================================
## Design sensitivities (gating thresholds), keyed
## DS[[curve]][[kappa]][[method]] = Gamma_bar_*.
## ============================================================================
message("Computing design sensitivities (gating thresholds) ...")
DS <- list()
for (cv in CURVES_TO_RUN) {
  DS[[cv]] <- list()
  pop <- curves[[cv]]$pop
  for (kn in KAPPAS_TO_RUN) {
    kf <- kappas[[kn]]
    DS[[cv]][[kn]] <- list()
    for (tt in tests) {
      DS[[cv]][[kn]][[tt$method]] <- design_sens(pop, tt$method, kf)
    }
  }
}


## The user can check the design sensitivities for all of the settings 
DS

## or can check the design sensitivities per dose-response curve
DS$Kink

## or can check the design sensitivity of a test give a \kappa(z) test design
## transformation under a dose-response curve

DS$Kink$`log(z)`$dose.weighted


## ============================================================================
## Bahadur-Rosenbaum slope Upsilon for each test, curve, kappa, reproduced 
## Table S.10, Table S.11, Table S.12 for the Bahadur Rosenbaum slope,
## And by theorem, the Bahadur Rosenbaum slope for the adaptive test is the 
## highest one of the component tests, we do not compute here.
## Each row of the slope_long follows the format of
##
##      curve   kappa Gamma_bar     test Gamma_star     Upsilon highest_Upsilon_test
## Square       z       1.0   Wilcox   2.320736 0.051959079             D-Wilcox
##       ...      ...       ...      ...        ...         ...                  ...
##       ...      ...       ...      ...        ...         ...                  ...
##
## where the curve represents the dose-response curve, the kappa represents the 
## \kappa(z) taken, the Gamma_bar is the generalized sensitivity parameter for testing
## the Gamma_star is the design sensitivity for this particular test, for this \kappa 
## choice under this dose-response curve, therefore, we make sure that the Gamma_bar 
## is smaller than the \Gamma_star,
## the Upsilon is the Bahadur Rosenbaum slope, and 
## highest_Upsilon_test is the test with the highest Upsilon fixing the dose-response curve,
## the Gamma_bar and the kappa(z) choice among Wilcox, D-Wilcox, U, and D-U tests.
## The first row of "slope_long" then recovers the cell of Dose-response Square, 
## \kappa(z) = z, \Gamma_bar = 1.00, wilcox 
## ============================================================================
if (RUN_PART_B) {
  message("PART B: computing Bahadur slopes ...")
  rows <- list(); k <- 0L
  for (cv in CURVES_TO_RUN) {
    pop <- curves[[cv]]$pop
    for (kn in KAPPAS_TO_RUN) {
      kf <- kappas[[kn]]
      for (g in GAMMA_GRID) {
        slopes <- setNames(rep(NA_real_, length(tests)),
                           vapply(tests, function(t) t$label, ""))
        for (tt in tests) {
          Gstar <- DS[[cv]][[kn]][[tt$method]]
          slopes[tt$label] <- slope_guarded(pop, tt$method, g, kf, Gstar = Gstar)
        }
        highest <- if (all(is.na(slopes))) NA_character_ else names(which.max(slopes))
        for (tt in tests) {
          k <- k + 1L
          rows[[k]] <- data.frame(
            curve   = cv, kappa = kn, Gamma_bar = g,
            test    = tt$label,
            Gamma_star = DS[[cv]][[kn]][[tt$method]],
            Upsilon = slopes[tt$label],
            highest_Upsilon_test = highest,
            stringsAsFactors = FALSE, row.names = NULL)
        }
      }
    }
  }
  slope_long <- do.call(rbind, rows)
  
  ## Wide table per curve: rows = (kappa, Gamma_bar), columns = the four slopes.
  make_slope_table <- function(curve_name, digits = 3) {
    sub <- slope_long[slope_long$curve == curve_name, ]
    wide <- reshape(sub[, c("kappa","Gamma_bar","test","Upsilon")],
                    idvar = c("kappa","Gamma_bar"), timevar = "test",
                    direction = "wide")
    names(wide) <- sub("Upsilon\\.", "", names(wide))
    wide <- wide[order(match(wide$kappa, KAPPAS_TO_RUN), wide$Gamma_bar), ]
    num <- vapply(wide, is.numeric, logical(1)); num["Gamma_bar"] <- FALSE
    wide[num] <- lapply(wide[num], round, digits)
    rownames(wide) <- NULL
    wide
  }
  
  cat("\n================ PART B: Bahadur slopes (Upsilon) ================\n")
  for (cv in CURVES_TO_RUN) {
    cat("\n--- Dose-response:", cv, "--- (NA: Gamma_bar >= design sensitivity)\n")
    print(make_slope_table(cv))
  }
  write.csv(slope_long, "bahadur_slopes_by_kappa_long.csv", row.names = FALSE)
}

## ============================================================================
## PART A  -- Verification: predicted relative efficiency vs simulated size ratio.
##  This part gives Table S.9 in the supplementary material 
##   Two comparisons per (curve, kappa)
##     Pair 1: Wilcox vs D-Wilcox(kappa)
##     Pair 2: U      vs D-U(kappa)
##   predicted = Upsilon_nondose / Upsilon_dose
##   simulated = I_dose / I_nondose        
##   Both gated on Gamma_bar < min(Gamma_*_nondose, Gamma_*_dose).
##
## ============================================================================
if (RUN_PART_A) {
  message("PART A: verification (simulated minimal sample sizes) ...")
  
  pairs <- list(
    list(name = "Wilcox vs D-Wilcox", nondose = "wilcoxon",
         dose = "dose.weighted",   nondose_lab = "Wilcox", dose_lab = "D-Wilcox"),
    list(name = "U vs D-U",          nondose = "U",
         dose = "dose.weighted.U", nondose_lab = "U",      dose_lab = "D-U")
  )
  
  total   <- length(CURVES_TO_RUN) * length(KAPPAS_TO_RUN) *
    length(GAMMA_GRID) * length(pairs)
  out_csv <- "bahadur_relative_efficiency_verification_long.csv"
  fmtI    <- function(x) if (is.na(x)) "  NA" else sprintf("%4d", as.integer(x))
  wrote_header  <- FALSE
  VERBOSE_SEARCH <- TRUE
  
  cat(sprintf("\nPART A: %d cells to compute.\n", total))
  cat("predicted = Upsilon_nondose/Upsilon_dose ;  simulated = I_dose/I_nondose\n")
  cat("(results appended to ", out_csv, " as they complete)\n\n", sep = "")
  flush.console()
  
  rows <- list(); k <- 0L
  for (cv in CURVES_TO_RUN) {
    pop <- curves[[cv]]$pop; ub <- curves[[cv]]$ub
    for (kn in KAPPAS_TO_RUN) {
      kf <- kappas[[kn]]
      for (g in GAMMA_GRID) {
        for (pr in pairs) {
          k <- k + 1L
          t0 <- Sys.time()
          cat(sprintf("[%3d/%3d] %-11s kappa=%-7s Gamma=%.2f  %-18s\n",
                      k, total, cv, kn, g, pr$name)); flush.console()
          
          Gstar_nd <- DS[[cv]][[kn]][[pr$nondose]]
          Gstar_d  <- DS[[cv]][[kn]][[pr$dose]]
          valid    <- is.finite(Gstar_nd) && is.finite(Gstar_d) &&
            (g < Gstar_nd) && (g < Gstar_d)
          
          predicted <- simulated <- I_nd <- I_d <- NA_real_
          if (valid) {
            Ups_nd <- slope_guarded(pop, pr$nondose, g, kf, Gstar = Gstar_nd)
            Ups_d  <- slope_guarded(pop, pr$dose,    g, kf, Gstar = Gstar_d)
            predicted <- Ups_nd / Ups_d
            
            I_nd <- min_pairs_fun(curves[[cv]]$gen, method = pr$nondose,
                                  Gamma_ave = g, kappa = kf, upper.bound = ub,
                                  verbose = VERBOSE_SEARCH,
                                  tag = sprintf("[%s]", pr$nondose_lab))
            I_d  <- min_pairs_fun(curves[[cv]]$gen, method = pr$dose,
                                  Gamma_ave = g, kappa = kf, upper.bound = ub,
                                  verbose = VERBOSE_SEARCH,
                                  tag = sprintf("[%s]", pr$dose_lab))
            simulated <- I_d / I_nd
          }
          secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          
          if (valid) {
            cat(sprintf("    -> pred = %6.3f   sim = %6.3f   (I_%s = %s, I_%s = %s)   [%5.1fs]\n",
                        predicted, simulated, pr$nondose_lab, fmtI(I_nd),
                        pr$dose_lab, fmtI(I_d), secs))
          } else {
            cat(sprintf("    -> NA: Gamma_bar >= design sensitivity   [%5.1fs]\n", secs))
          }
          flush.console()
          
          row <- data.frame(
            curve = cv, kappa = kn, Gamma_bar = g, comparison = pr$name,
            Gamma_star_nondose = Gstar_nd, Gamma_star_dose = Gstar_d,
            predicted_rel_eff = predicted,
            I_nondose = I_nd, I_dose = I_d,
            simulated_size_ratio = simulated,
            stringsAsFactors = FALSE, row.names = NULL)
          rows[[k]] <- row
          
          write.table(row, out_csv, sep = ",", row.names = FALSE,
                      col.names = !wrote_header, append = wrote_header,
                      qmethod = "double")
          wrote_header <- TRUE
        }
      }
    }
  }
  releff_long <- do.call(rbind, rows)
  
  fmt_cell <- function(p, s) {
    if (is.na(p) && is.na(s)) return("NA")
    sprintf("%.2f(%.2f)", p, s)
  }
  make_releff_table <- function(curve_name) {
    sub <- releff_long[releff_long$curve == curve_name, ]
    sub$cell <- mapply(fmt_cell, sub$predicted_rel_eff, sub$simulated_size_ratio)
    wide <- reshape(sub[, c("kappa","Gamma_bar","comparison","cell")],
                    idvar = c("kappa","Gamma_bar"), timevar = "comparison",
                    direction = "wide")
    names(wide) <- sub("cell\\.", "", names(wide))
    wide <- wide[order(match(wide$kappa, KAPPAS_TO_RUN), wide$Gamma_bar), ]
    rownames(wide) <- NULL
    wide
  }
  
  cat("\n=========== PART A: VERIFICATION  predicted(simulated) ===========\n")
  cat("predicted = Upsilon_nondose / Upsilon_dose ;  simulated = I_dose / I_nondose\n")
  for (cv in CURVES_TO_RUN) {
    cat("\n--- Dose-response:", cv, "--- (NA: Gamma_bar >= a design sensitivity)\n")
    print(make_releff_table(cv))
  }
  write.csv(releff_long, "bahadur_relative_efficiency_verification_long.csv",
            row.names = FALSE)
}

message("Done. CSV outputs (if enabled): ",
        "bahadur_slopes_by_kappa_long.csv, ",
        "bahadur_relative_efficiency_verification_long.csv")

## can call releff_long to show the results on the screen 
releff_long