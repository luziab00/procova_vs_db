source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")
source("analysis/procova_analysis.R")
source("analysis/comm_prior.R")
source("analysis/SAM_prior.R")
source("analysis/adaptive_lasso.R")
source("analysis/bayesian_procova.R")


covar_distr_hist <- list(
  mu  = c(X1=0, X2=0, X3=0, X4=0),
  sd  = c(X1=1, X2=1, X3=1, X4=1),
  pX5 = 0.5
)

covar_distr_current <- list(
  mu  = c(X1=0, X2=0, X3=0, X4=0),
  sd  = c(X1=1, X2=1, X3=1, X4=1),
  pX5 = 0.5
)


coeff_hist <-  list(
  beta0  = 0,
  beta1  = 0.6, beta2  = 0.3, beta3  = 0.2, beta4  = 0.4,
  gamma1 = 0.5, gamma2 = 0.4, gamma3 = 0.6, gamma4 = 0.3, gamma5=0.5
)

coeff_current <-  list(
  beta0  = 0,
  beta1  = 0.6, beta2  = 0.3, beta3  = 0.2, beta4  = 0.4,
  gamma1 = 0.5, gamma2 = 0.4, gamma3 = 0.6, gamma4 = 0.3, gamma5=0.5
)

treatment_effect <- 0


nsim <- 1
N_hist_default    <- 200
N_current_default <- 200

# Choose your prognostic strength via pi and sigma:
sigma_total <- 1
R2_target <- 0.25
pi_value <- sqrt(R2_target) * sigma_total # since R^2 = pi^2 / sigma^2

seed<-1

# Compute population-level standardization ONCE for this setting:
score_std_setting <- compute_score_standardization(
  covars_distr = covar_distr_bl,
  coeff = coeff_hist,
  M = 200000,
  seed = seed + 9999
)




df_historical_trial <- simulate_single_trial_standardized(
  N = N_hist_default,
  covars_distr = covar_distr_hist,
  coeff = coeff_hist,
  sigma = sigma_total,
  pi = pi_value,
  score_std = score_std_setting,
  trial_name = "historical",
  treatment_arm = FALSE,
  tau = 0,
  randomization_ratio = NULL
)

df_current_trial <- simulate_single_trial_standardized(
  N = N_current_default,
  covars_distr = covar_distr_current,
  coeff = coeff_current,
  sigma = sigma_total,
  pi = pi_value,
  score_std = score_std_setting,
  trial_name = "current",
  treatment_arm = TRUE,
  tau = treatment_effect,
  randomization_ratio = 0.5
)

bayesian_procova_fit <- bayesian_procova_analysis(
  df_historical = df_historical_trial,
  df_current = df_current_trial,
  model_fun = prognostic_model_superlearn,
  seed = seed + 300000
)





tic()
fit_jags_superlearner   <- bayesian_procova_jags(df_historical_trial, df_current_trial, oracle=FALSE, K0H=1/sqrt(200),  model_fun = prognostic_model_superlearn)
toc()


tic()
fit_custom_superlearner <- bayesian_procova_analysis(df_historical_trial, df_current_trial, oracle=FALSE, K0H=1/sqrt(200), model_fun = prognostic_model_superlearn)
toc()



tic()
fit_jags_oracle   <- bayesian_procova_jags(df_historical_trial, df_current_trial, oracle=TRUE, K0H=1/sqrt(200))
toc()


tic()
fit_custom_oracle <- bayesian_procova_analysis(df_historical_trial, df_current_trial, oracle=TRUE, K0H=1/sqrt(200))
toc()



summ_treat_bayesian_procova(fit_jags)

summ_treat_bayesian_procova(fit_custom)



# 1. Compare posterior omega distributions
plot(density(fit_custom_oracle$draws$omega), main = "omega posterior")
lines(density(fit_jags_oracle$draws$omega), col = "red")

# 2. Compare Z frequencies
mean(fit_custom_oracle$draws$Z)   # fraction of draws from informative component
mean(fit_jags_oracle$draws$Z)


mean(fit_custom_oracle$draws$sigma2)


# 1. Compare posterior omega distributions
plot(density(fit_custom_superlearner$draws$omega), main = "omega posterior")
lines(density(fit_jags_superlearner$draws$omega), col = "red")

# 2. Compare Z frequencies
mean(fit_custom_superlearner$draws$Z)   # fraction of draws from informative component
mean(fit_jags_superlearner$draws$Z)


mean(fit_custom_oracle$draws$sigma2)



par(mfrow = c(1, 2))

# Oracle
plot(density(fit_custom_oracle$draws$omega), 
     main = "Oracle: omega posterior", col = "black")
lines(density(fit_jags_oracle$draws$omega), col = "red")
legend("topright", c("custom", "jags"), col = c("black","red"), lty = 1)

# Superlearner  
plot(density(fit_custom_superlearner$draws$omega),
     main = "Superlearner: omega posterior", col = "black")
lines(density(fit_jags_superlearner$draws$omega), col = "red")
legend("topright", c("custom", "jags"), col = c("black","red"), lty = 1)



par(mfrow = c(1, 2))

# Oracle
plot(density(fit_custom_oracle$draws$beta[,"beta1"]),
     main = "Oracle: beta1 posterior", col = "black")
lines(density(fit_jags_oracle$draws$beta[,"beta1"]), col = "red")

# Superlearner
plot(density(fit_custom_superlearner$draws$beta[,"beta1"]),
     main = "Superlearner: beta1 posterior", col = "black")
lines(density(fit_jags_superlearner$draws$beta[,"beta1"]), col = "red")

# Under drift, Z should drop toward 0 for both implementations
# The drop should be more pronounced for oracle (more sensitive to drift
# because it borrows more aggressively under no drift)
cat("Z mean under no drift, oracle custom:     ", mean(fit_custom_oracle$draws$Z), "\n")
cat("Z mean under no drift, superlearner custom:", mean(fit_custom_superlearner$draws$Z), "\n")

# because it borrows more aggressively under no drift)
cat("Z mean under no drift, oracle jags:     ", mean(fit_jags_oracle$draws$Z), "\n")
cat("Z mean under no drift, superlearner jags:", mean(fit_jags_superlearner$draws$Z), "\n")









# Decompose q_I into contributions from each source
r_I_raw <- y_c - V %*% beta_H

# How much of q_I comes from the slope mismatch vs intercept vs noise?
# Approximate decomposition:
cat("=== q_I DECOMPOSITION ===\n")
cat("K2H current:", K2H, "\n")
cat("beta_H[3]:", beta_H[3], "\n")

# RCT slope
rct_fit <- lm(y_c ~ I(m - mean(m)))
rct_slope <- coef(rct_fit)[2]
cat("RCT slope:", rct_slope, "\n")
cat("Slope discrepancy:", beta_H[3] - rct_slope, "\n")

# How much does the slope discrepancy contribute to q_I?
# Replace beta_H[3] with RCT slope and recompute q_I
beta_H_corrected <- beta_H
beta_H_corrected[3] <- rct_slope
r_I_corrected <- y_c - V %*% beta_H_corrected
q_I_corrected <- sum(r_I_corrected * (M_I_inv %*% r_I_corrected))

cat("q_I original:  ", q_I, "\n")
cat("q_I corrected slope:", q_I_corrected, "\n")
cat("q_I increase from slope mismatch:", q_I - q_I_corrected, "\n")

# What log_I would be with corrected slope?
log_I_corrected <- lgamma((N + N_H - 2)/2) - lgamma((N_H-2)/2) -
  N/2 * log((N_H-2)*pi*s2_H) - 0.5*log_det_M_I -
  (N + N_H - 2)/2 * log(1 + q_I_corrected/((N_H-2)*s2_H))
cat("log_I original:  ", log_I, "\n")
cat("log_I corrected: ", log_I_corrected, "\n")
cat("log_I - log_F corrected:", log_I_corrected - log_F, "\n")


# The prior variance on beta_H[3] needs to cover the slope discrepancy
slope_discrepancy <- beta_H[3] - rct_slope
required_K2H <- slope_discrepancy^2 / s2_H
cat("Required K2H to tolerate slope discrepancy:", required_K2H, "\n")

# Compare to current value
cat("Current K2H:", K2H, "\n")
cat("Ratio:", required_K2H / K2H, "\n")





# Try a range that covers the required magnitude
k2h_grid <- c(
  K2H,                    # current attempt
  required_K2H * 0.5,          # just below required
  required_K2H,                # exact match
  required_K2H * 2,            # above required
  K1H                          # match the treatment effect prior — fully agnostic
)

results <- lapply(k2h_grid, function(k2h) {
  fit <- bayesian_procova_analysis(
    df_historical_trial, df_current_trial,
    model_fun = prognostic_model_superlearn,
    oracle = FALSE, K2H = k2h
  )
  data.frame(
    K2H        = k2h,
    Z          = mean(fit$draws$Z),
    beta1_mean = mean(fit$draws$beta[,"beta1"]),
    beta1_sd   = sd(fit$draws$beta[,"beta1"])
  )
})

do.call(rbind, results)


cat("=== FULL q_I DECOMPOSITION ===\n")

# 1. What is q_I if we use OLS fit to RCT data as beta_H?
rct_fit_full <- lm(y_c ~ w + I(m - mean(m)))
beta_rct <- c(coef(rct_fit_full)[1], 
              coef(rct_fit_full)[2], 
              coef(rct_fit_full)[3])
r_I_rct <- y_c - V %*% beta_rct
q_I_rct <- sum(r_I_rct * (M_I_inv %*% r_I_rct))
cat("q_I with RCT beta (best case):", q_I_rct, "\n")

# 2. Intercept contribution
beta_H_int_fix <- beta_H
beta_H_int_fix[1] <- coef(lm(y_c ~ 1))[1]  # RCT intercept
r_I_int <- y_c - V %*% beta_H_int_fix
q_I_int <- sum(r_I_int * (M_I_inv %*% r_I_int))
cat("q_I with corrected intercept:", q_I_int, "\n")
cat("q_I increase from intercept mismatch:", q_I - q_I_int, "\n")  

# 3. Pure noise floor — what q_I would be with perfect beta
cat("Pure noise q_I (RCT OLS residuals):", q_I_rct, "\n")
cat("Expected q_I under correct model (≈N):", N, "\n")

# 4. log_I with best-case beta
log_I_bestcase <- lgamma((N + N_H - 2)/2) - lgamma((N_H-2)/2) -
  N/2 * log((N_H-2)*pi*s2_H) - 0.5*log_det_M_I -
  (N + N_H - 2)/2 * log(1 + q_I_rct/((N_H-2)*s2_H))
cat("log_I best case:", log_I_bestcase, "\n")
cat("log_I - log_F best case:", log_I_bestcase - log_F, "\n")



# Set sigma2_0 to RCT residual variance estimate
rct_sigma2_est <- var(residuals(lm(y_c ~ w + I(m - mean(m)))))
cat("RCT sigma2 estimate:", rct_sigma2_est, "\n")

# Or more conservatively, use a larger nu0 with matched scale
fit_sl_fixed <- bayesian_procova_analysis(
  df_historical_trial, df_current_trial,
  model_fun = prognostic_model_superlearn,
  oracle = FALSE,
  K2H = 0.13,
  nu0 = 3,
  sigma2_0 = rct_sigma2_est  # match actual variance scale
)
cat("Z mean:", mean(fit_sl_fixed$draws$Z), "\n")


# Use historical residual variance as sigma2_0
sigma2_0_principled <- s2_H
cat("s2_H as sigma2_0:", sigma2_0_principled, "\n")

fit_sl_principled <- bayesian_procova_analysis(
  df_historical_trial, df_current_trial,
  model_fun = prognostic_model_superlearn,
  oracle = FALSE,
  nu0 = 3,
  sigma2_0 = s2_H   # principled: set from historical data
)
cat("Z mean:", mean(fit_sl_principled$draws$Z), "\n")



cat("\n=== log_I components ===\n")
cat("lgamma term:    ", lgamma((N+N_H-2)/2) - lgamma((N_H-2)/2), "\n")
cat("log scale term: ", -N/2 * log((N_H-2)*pi*s2_H), "\n")
cat("log det term:   ", -0.5*log_det_M_I, "\n")
cat("quadratic term: ", -(N+N_H-2)/2 * log(1 + q_I/((N_H-2)*s2_H)), "\n")
cat("log_I total:    ", log_I, "\n")

cat("\n=== log_F components ===\n")
cat("lgamma term:    ", lgamma((N+nu0)/2) - lgamma(nu0/2), "\n")
cat("log scale term: ", -N/2 * log(nu0*pi*sigma2_0), "\n")
cat("log det term:   ", -0.5*log_det_M_F, "\n")
cat("quadratic term: ", -(N+nu0)/2 * log(1 + q_F/(nu0*sigma2_0)), "\n")
cat("log_F total:    ", log_F, "\n")

cat("\n=== key quantities ===\n")
cat("N:", N, "\n")
cat("N_H:", N_H, "\n")
cat("nu0:", nu0, "\n")
cat("sigma2_0:", sigma2_0, "\n")
cat("s2_H:", s2_H, "\n")
cat("k:", k, "\n")
cat("q_F:", q_F, "\n")
cat("q_F/(nu0*sigma2_0):", q_F/(nu0*sigma2_0), "\n")
cat("log_det_M_I:", log_det_M_I, "\n")
cat("log_det_M_F:", log_det_M_F, "\n")



fit_sl_scaled <- bayesian_procova_analysis(
  df_historical_trial, df_current_trial,
  model_fun = prognostic_model_superlearn,
  oracle = FALSE,
  nu0 = 1,
  sigma2_0 = (N_H - 2) * s2_H   # ≈ 74
)
cat("Z mean:", mean(fit_sl_scaled$draws$Z), "\n")

# Also try nu0=3 version
fit_sl_scaled_nu3 <- bayesian_procova_analysis(
  df_historical_trial, df_current_trial,
  model_fun = prognostic_model_superlearn,
  oracle = TRUE,
  nu0 = 3,
  sigma2_0 = (N_H - 2) * s2_H / 3   # ≈ 24.7
)
cat("Z mean nu0=3:", mean(fit_sl_scaled_nu3$draws$Z), "\n")



compute_hyperparameters <- function(
  df_historical,
  model_covar_names,
  nu0 = 3,
  N_current,        # expected RCT sample size
  n_boot = 1000,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  y_H <- df_historical$y
  m_H <- df_historical$progn_score
  N_H <- nrow(df_historical)
  
  m_bar_H <- mean(m_H)
  fit_H   <- lm(I(y_H - m_bar_H) ~ I(m_H - m_bar_H))
  s2_H    <- sum(residuals(fit_H)^2) / (N_H - 2)
  beta_H  <- c(coef(fit_H)[1], 0, coef(fit_H)[2])
  
  # -----------------------------------------------------------
  # K0H: bootstrap intercept variability at RCT sample size
  # (Appendix B of paper)
  # -----------------------------------------------------------
  boot_intercepts <- replicate(n_boot, {
    idx    <- sample(N_H, N_current, replace = TRUE)
    m_b    <- m_H[idx]; y_b <- y_H[idx]
    m_bar_b <- mean(m_b)
    coef(lm(I(y_b - m_bar_b) ~ I(m_b - m_bar_b)))[1]
  })
  var_delta <- var(boot_intercepts)
  K0H <- 1/N_H + var_delta/s2_H
  
  # -----------------------------------------------------------
  # K2H: bootstrap slope variability at RCT sample size
  # Same logic as K0H but for the prognostic score slope
  # -----------------------------------------------------------
  boot_slopes <- replicate(n_boot, {
    idx    <- sample(N_H, N_current, replace = TRUE)
    m_b    <- m_H[idx]; y_b <- y_H[idx]
    m_bar_b <- mean(m_b)
    coef(lm(I(y_b - m_bar_b) ~ I(m_b - m_bar_b)))[2]
  })
  var_slope <- var(boot_slopes)
  K2H <- 1/sum((m_H - m_bar_H)^2) + var_slope/s2_H
  
  # -----------------------------------------------------------
  # sigma2_0: match flat component scale to informative component
  # so both components are evaluated on comparable scales
  # -----------------------------------------------------------
  sigma2_0 <- (N_H - 2) * s2_H / nu0
  
  # -----------------------------------------------------------
  # K1H: keep large (weakly informative on treatment effect)
  # -----------------------------------------------------------
  K1H <- 100
  
  list(
    K0H      = K0H,
    K1H      = K1H,
    K2H      = K2H,
    nu0      = nu0,
    sigma2_0 = sigma2_0,
    s2_H     = s2_H,
    beta_H   = beta_H,
    # diagnostics
    var_delta_boot = var_delta,
    var_slope_boot = var_slope
  )
}


hyperparams_oracle <- compute_hyperparameters(
  df_historical = df_historical_with_oracle_scores,
  N_current = N_current_default,
  nu0 = 3,
  seed = seed
)

hyperparams_sl <- compute_hyperparameters(
  df_historical = df_historical_with_sl_scores,
  N_current = N_current_default,
  nu0 = 3,
  seed = seed
)