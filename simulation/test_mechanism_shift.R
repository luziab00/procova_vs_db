source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")
source("analysis/procova_analysis.R")
source("analysis/db_analysis.R")
source("analysis/bayesian_procova.R")


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
  beta1 = 0.5,
  beta2 = -0.5,
  beta3 = 0.5,
  beta4 = 1,
  gamma1 = 0.5,
  gamma2 = -0.5,
  gamma3 = 0.5,
  gamma4 = 1
)


coeff_current <- list(
  beta0 = 0,
  beta1 = 0.5,
  beta2 = -0.5,
  beta3 = -0.5,
  beta4 = 0.5,
  gamma1 = -0.5,
  gamma2 = 0.5,
  gamma3 = 1,
  gamma4 = 1
)

treatment_effect <- 0


nsim <- 1
seed <- 123


# Choose your prognostic strength via pi and sigma:
sigma_total <- 1
R2_target <- 0.25
pi_value <- sqrt(R2_target) * sigma_total # since R^2 = pi^2 / sigma^2


df_historical_trial <- simulate_single_trial_standardized(
  N = 10000,
  covars_distr = covar_distr_hist,
  coeff = coeff_hist,
  sigma = sigma_total,
  pi = pi_value,
  score_std = score_std_setting,
  trial_name = "historical",
  treatment_arm = FALSE,
  tau = 0,
  randomization_ratio = NULL,
  seed = seed + 100000 + i
)

df_current_trial <- simulate_single_trial_standardized(
  N = 10000,
  covars_distr = covar_distr_current,
  coeff = coeff_current,
  sigma = sigma_total,
  pi = pi_value,
  score_std = score_std_setting,
  trial_name = "current",
  treatment_arm = TRUE,
  tau = treatment_effect,
  randomization_ratio = 0.5,
  seed = seed + 200000 + i
)




df_current_trial$S_raw_old<-dgm_score_raw(df_current_trial, coeff=coeff_hist)

ggplot(df_current_trial, aes(x=S_raw, y=S_raw_old))+
  geom_point()


cor(df_current_trial$S_raw, df_current_trial$S_raw_old)

