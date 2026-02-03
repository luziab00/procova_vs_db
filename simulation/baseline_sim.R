source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")
source("analysis/procova_analysis.R")



covar_distr_bl <- list(
  mu = c(X1 = 0, X2 = 0, X3 = 0, X4 = 0),
  sd = c(X1 = 1, X2 = 1, X3 = 1, X4 = 1),
  pX5 = 0.5
)


coeff_hist <- list(
  beta0 = 0,
  beta1 = 0.5,
  beta2 = -1,
  beta3 = -0.5,
  beta4 = 1,
  gamma1 = 0.5,
  gamma2 = 0.5,
  gamma3 = -0.5,
  gamma4 = 1
)


coeff_current <- list(
  beta0 = 0,
  beta1 = 0.5,
  beta2 = -1,
  beta3 = -0.5,
  beta4 = 1,
  gamma1 = 0.5,
  gamma2 = 0.5,
  gamma3 = -0.5,
  gamma4 = 1
)

treatment_effect <- 0.5

df_historical_trial <- simulate_single_trial(
  N = 200,
  covars_distr=covar_distr_bl ,
  coeff=coeff_hist, 
  residual_error_sd = 0.25,
  trial_name       = "historical",
  treatment_arm    = FALSE,   # boolean, TRUE if treatment arm included
  tau              = NULL,    # numeric, treatment effect (only used if treatment_arm=TRUE)
  randomization_ratio = NULL, # numeric in (0,1), probability of XTreat=1
  seed = NULL                 # optional integer for reproducibility
)

df_current_trial <- simulate_single_trial(
  N=200,
  covars_distr=covar_distr_bl ,
  coeff=coeff_current, 
  residual_error_sd = 0.25,
  trial_name       = "current",
  treatment_arm    = TRUE,   # boolean, TRUE if treatment arm included
  tau              = treatment_effect,    # numeric, treatment effect (only used if treatment_arm=TRUE)
  randomization_ratio = 0.5, # numeric in (0,1), probability of XTreat=1
  seed = NULL                 # optional integer for reproducibility
)




procova_glm<-procova_analysis_glm(
  df_historical=df_historical_trial, 
  df_current=df_current_trial, 
  model=prognostic_model_oracle_lm, 
  model_covar_names = c("X1","X2","X3","X4", "X5"), 
  seed = NULL)

procova_glm$fit_anova %>% summ_treat_glm()
procova_glm$fit_ancova %>% summ_treat_glm()
procova_glm$fit_procova %>% summ_treat_glm()


procova_stan<-procova_analysis_stan(
  df_historical=df_historical_trial, 
  df_current=df_current_trial, 
  model=prognostic_model_superlearn, 
  model_covar_names = c("X1","X2","X3","X4", "X5"), 
  seed = NULL)

procova_stan$fit_anova %>% summ_treat_stan_brms()
procova_stan$fit_ancova %>% summ_treat_stan_brms()
procova_stan$fit_procova %>% summ_treat_stan_brms()

summary_df<- rbind(
      cbind(method="glm",  model="ANOVA", t(summ_treat_glm(procova_glm$fit_anova))),
      cbind(method="glm",  model="ANCOVA",  t(summ_treat_glm(procova_glm$fit_ancova))),
      cbind(method="glm",  model="PROCOVA", t(summ_treat_glm(procova_glm$fit_procova))),
      cbind(method="stan",  model="ANOVA",   t(summ_treat_stan_brms(procova_stan$fit_anova))),
      cbind(method="stan",  model="ANCOVA",  t(summ_treat_stan_brms(procova_stan$fit_ancova))),
      cbind(method="stan",  model="PROCOVA", t(summ_treat_stan_brms(procova_stan$fit_procova)))
    )


 # combine the dfs
df_all_trial<-rbind(cbind(df_current_trial, ext=0), cbind(df_historical_trial, ext=1))


ggplot(df_all_trial, aes(x=ext, y=y, color=XTreat))+geom_point()

comm_db_psborrow2 <-   comm_prior_db_analysis_psborrow2(
    df_historical_trial, 
    df_current_trial, 
    seed = seed, 
    tau_rate = 0.01
  )

summ_treat_comm_psborrow2(comm_db_psborrow2$fit_comm)


#draws_df <- posterior::as_draws_df(comm_db_psborrow2$fit$draws())
#draws <- as.numeric(draws_df[["beta_trt"]])
#hist(draws, breaks=50, xlim=c(1, 2))


adaptive_lasso_db <- adaptive_lasso_db_analysis(
  df_current_trial, 
  df_historical_trial, 
  model_covar_names = c("X1","X2","X3","X4","X5"), 
  lambda = 2,
  gamma = 1.5)

summ_treat_adaptive_lasso(adaptive_lasso_db$fit)


samprior_fit<-sam_rbest_fit(
  df_historical=df_historical_trial,
  df_current=df_current_trial,
  y_col = "y",
  trt_col = "XTreat",
  delta_mult = 0.5)

sam_rbest_summary(samprior_fit)

comm_prior_db_fit_RBest(
  df_historical=df_historical_trial,
  df_current=df_current_trial,
  seed = NULL,
  tau_rate = 0.01,
  tau_grid_m = 200,
  huge_sd = 1e6
) 



# Pure oracle score: returns a function m(x) = E[Y|X, T=0] using known coefficients
prognostic_score_oracle_true <- function(coeff) {
  force(coeff)
  function(newdata) {
    with(newdata,
      coeff$beta0 +
        coeff$beta1 * X1 + coeff$beta2 * X2 + coeff$beta3 * X3 + coeff$beta4 * X4 +
        coeff$gamma1 * (X1^2) + coeff$gamma2 * sin(X3) + coeff$gamma3 * (X4 > 1) +
        coeff$gamma4 * (X1 * X5)
    )
  }
}


# 1) Create the oracle score function (true m(x))
m_oracle <- prognostic_score_oracle_true(coeff_hist)

# 2) Monte Carlo estimate Var{m(X)} under your covariate distribution
set.seed(1)
M <- 200000
X <- gen_covs(M, covar_distr_bl)
mX <- m_oracle(X)

var_m <- var(mX)
sigma <- 2  # your residual_error_sd
R2_prog <- var_m / (var_m + sigma^2)

R2_prog
