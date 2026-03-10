covar_distr_hist <- list(
  mu = c(X1 = 0, X2 = 0, X3 = 0, X4 = 0),
  sd = c(X1 = 1, X2 = 1, X3 = 1, X4 = 1),
  pX5 = 0.5
)

covar_distr_current <- list(
  mu = c(X1 = 0, X2 = 0, X3 = 0, X4 = 0),
  sd = c(X1 = 1, X2 = 1, X3 = 1, X4 = 1),
  pX5 = 0.5
)



coeff_hist <- list(
  beta0 = 0,
  beta1 = 0.6,
  beta2 = 0.3,
  beta3 = 0.2,
  beta4 = 0.4,
  gamma1 = 0.5,
  gamma2 = 0.4,
  gamma3 = 0.6,
  gamma4 = 0.3
  )


coeff_current <- list(
  beta0 = 0,
  beta1 = 0.6,
  beta2 = 0.3,
  beta3 = 0.2,
  beta4 = 0.4,
  gamma1 = 0.5,
  gamma2 = 0.4,
  gamma3 = 0.6,
  gamma4 = 0.3
)

treatment_effect <- 0


nsim <- 1
seed <- 123


# Choose your prognostic strength via pi and sigma:
sigma_total <- 1
R2_target   <- 0.25
pi_value    <- sqrt(R2_target) * sigma_total   # since R^2 = pi^2 / sigma^2