

covar_distr_bl <- list(
  mu = c(X1 = 0, X2 = 0, X3 = 0, X4 = 0),
  sd = c(X1 = 1, X2 = 1, X3 = 1, X4 = 1),
  pX5 = 0.5
)


coeff <- list(
  beta1 = 1,
  beta2 = -1,
  beta3 = 0,
  beta4 = 2,
  gamma1 = 0.5,
  gamma2 = 0.5,
  gamma3 = 0.5,
  gamma4 = 0.5
)


gen_covs <- function(N, covars) {
  X1 <- rnorm(N, covars$mu["X1"], covars$sd["X1"])
  X2 <- rnorm(N, covars$mu["X2"], covars$sd["X2"])
  X3 <- rnorm(N, covars$mu["X3"], covars$sd["X3"])
  X4 <- rnorm(N, covars$mu["X4"], covars$sd["X4"])
  X5  <- rbinom(N, 1, covars$pX5)
  data.frame(X1, X2, X3, X4, X5)
}


# DGM for the continuous outcome
# df must contain: X1, X2, X3, X4, X5, and a treatment indicator T_var (default "XTreat")

dgm_outcome <- function(
  df,
  coeff,                 # list: beta1..beta4, gamma1..gamma4
  tau   = 0,             # treatment effect 
  residual_error_sd = 0.05,          # residual SD σ
  seed  = NULL           # optional seed
) {
  if (!is.null(seed)) set.seed(seed)

  # extract treatment indicator
  Tvec <- df[["XTreat"]]

  # linear part: β1 X1 + β2 X2 + β3 X3 + β4 X4
  lin <- coeff$beta1 * df$X1 +
         coeff$beta2 * df$X2 +
         coeff$beta3 * df$X3 +
         coeff$beta4 * df$X4

  # nonlinear / interaction part:
  # γ1 X1^2 + γ2 sin(X3) + γ3 I(X4 > 1) + γ4 X1 X5
  nonlinear <- coeff$gamma1 * (df$X1^2) +
               coeff$gamma2 * sin(df$X3) +
               coeff$gamma3 * (df$X4 > 1) +
               coeff$gamma4 * (df$X1 * df$X5)

  
  if (is.null(tau)) tau <- 0
  # treatment contribution: tau
  treat_term <- tau * Tvec

  # linear predictor
  eta <- coeff$beta0 + lin + nonlinear + treat_term

  # add Normal noise
  eps <- rnorm(nrow(df), mean = 0, sd = residual_error_sd)

  y <- eta + eps
  y
}


simulate_single_trial <- function(
  N,                     # numeric/integer > 0, sample size of the dataset
  covars_distr,          # list describing covariate distribution:
                         #   covars_distr$mu  = c(X1 = ..., X2 = ..., X3 = ..., X4 = ...)
                         #   covars_distr$sd  = c(X1 = ..., X2 = ..., X3 = ..., X4 = ...)
                         #   covars_distr$pX4 = probability for binary X5
  coeff,                 # list of regression coefficients:
                         #   coeff$beta_lin
                         #   coeff$beta_quad
                         #   coeff$beta_sin
                         #   coeff$beta_interaction
  residual_error_sd = 0.05,   # numeric >= 0, SD of residual noise
  trial_name       = "historical", # character, label identifying the trial
  treatment_arm    = FALSE,   # boolean, TRUE if treatment arm included
  tau              = NULL,    # numeric, treatment effect (only used if treatment_arm=TRUE)
  randomization_ratio = NULL, # numeric in (0,1), probability of XTreat=1
  seed = NULL                 # optional integer for reproducibility
) {

  # generate the covariates
  df <- gen_covs(N, covars_distr)
  
  if (treatment_arm) {
    n_treated<- round(N * randomization_ratio)
    df$XTreat <- c(rep(1L, n_treated), rep(0L, N - n_treated))
  } else {
    df$XTreat <- 0L
  }

  # generate outcome based on covariates
  df$y <- dgm_outcome(
    df   = df,
    coeff = coeff,
    tau   = tau,
    residual_error_sd = residual_error_sd
  )

  df$trial_name<-trial_name


  df
}


dgm_score_raw <- function(df, coeff) {
  # Expect coeff to contain:
  #   coeff$beta_lin          named c(X1=..., X2=..., X3=..., X4=...)
  #   coeff$beta_quad         scalar for X1^2   (gamma1 in screenshot)
  #   coeff$beta_sin          scalar for sin(X3) (gamma2)
  #   coeff$beta_indicator    scalar for 1(X4 > 1) (gamma3)
  #   coeff$beta_interaction  scalar for X1*X5 (gamma4)

  coeff$beta1 * df$X1 +
  coeff$beta2 * df$X2 +
  coeff$beta3 * df$X3 +
  coeff$beta4 * df$X4
  coeff$gamma1 * (df$X1^2) +
  coeff$gamma2 * sin(df$X3) +
  coeff$gamma3 * (df$X4 > 1) +
  coeff$gamma4 * (df$X1 * df$X5)
}


compute_score_standardization <- function(covars_distr, coeff, M = 200000, seed = 1L) {

  if (!is.null(seed)) set.seed(seed)
  df_pop <- gen_covs(M, covars_distr)  # uses your existing covariate generator
  s_raw  <- dgm_score_raw(df_pop, coeff)

  list(
    mean = mean(s_raw),
    sd   = sd(s_raw)
  )
}


dgm_outcome_standardized <- function(
  df,
  s,           # standardized prognostic score
  coeff,
  tau = 0,
  pi,
  sigma
) {
  eps <- rnorm(nrow(df), mean = 0, sd = 1)

  y <- coeff$beta0 +
    pi * s +
    tau * df$XTreat +
    sqrt(sigma^2 - pi^2) * eps

  y
}


simulate_single_trial_standardized <- function(
  N,
  covars_distr,
  coeff,
  sigma = 2,
  pi = 1,
  score_std,
  trial_name = "historical",
  treatment_arm = FALSE,
  tau = 0,
  randomization_ratio = NULL,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # generate covariates
  df <- gen_covs(N, covars_distr)

  # treatment assignment
  if (treatment_arm) {
    n_treated <- round(N * randomization_ratio)
    df$XTreat <- c(rep(1L, n_treated), rep(0L, N - n_treated))
    df$XTreat <- sample(df$XTreat)
  } else {
    df$XTreat <- 0L
  }

  # compute standardized prognostic score ONCE
  s_raw <- dgm_score_raw(df, coeff)
  df$S <- (s_raw - score_std$mean) / score_std$sd

  # outcome
  df$y <- dgm_outcome_standardized(
    df = df,
    s = df$S,
    coeff = coeff,
    tau = tau,
    pi = pi,
    sigma = sigma
  )

  df$trial_name <- trial_name
  df
}
