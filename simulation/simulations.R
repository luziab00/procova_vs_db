source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")
source("analysis/procova_analysis.R")
source("analysis/db_analysis.R")


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


nsim <- 50
seed <- 123
set.seed(seed)

sims <- vector("list", nsim)


# Choose your prognostic strength via pi and sigma:
sigma_total <- 1
R2_target   <- 0.25
pi_value    <- sqrt(R2_target) * sigma_total   # since R^2 = pi^2 / sigma^2

# Compute population-level standardization ONCE for this setting:
score_std_setting <- compute_score_standardization(
  covars_distr = covar_distr_bl,
  coeff = coeff_hist,
  M = 200000,
  seed = seed + 9999
)


for (i in seq_len(nsim)) {

  df_historical_trial <- simulate_single_trial_standardized(
    N = 200,
    covars_distr = covar_distr_bl,
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
    N = 200,
    covars_distr = covar_distr_bl,
    coeff = coeff_hist,
    sigma = sigma_total,
    pi = pi_value,
    score_std = score_std_setting,
    trial_name = "current",
    treatment_arm = TRUE,
    tau = treatment_effect,
    randomization_ratio = 0.5,
    seed = seed + 200000 + i
  )

  sims[[i]] <- list(
    sim = i,
    df_historical = df_historical_trial,
    df_current = df_current_trial,
    true_tau = treatment_effect
  )
}

r_sqared_list <- vector("numeric", nsim)

for (i in seq_len(nsim)) {

  oracle_model<-lm(y ~ S, data = sims[[i]]$df_historical)

  r_sqared_list[i] <- summary(oracle_model)$r.squared
}

model_covar_names <- c("X1","X2","X3","X4","X5")

res_list <- vector("list", nsim)

for (i in seq_len(nsim)) {
  if (i %% 10 == 0) message("analyzing ", i, "/", nsim)

  df_h <- sims[[i]]$df_historical
  df_c <- sims[[i]]$df_current
  tau  <- sims[[i]]$true_tau

  # ----- Frequentist PROCOVA -----
  procova_glm <- procova_analysis_glm(
    df_historical = df_h,
    df_current = df_c,
    model = prognostic_model_superlearn,
    model_covar_names = model_covar_names,
    seed = seed + 300000 + i
  )

r_glm <- rbind(
  data.frame(
    method = "glm",
    model  = "ANOVA",
    t(summ_treat_glm(procova_glm$fit_anova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "glm",
    model  = "ANCOVA",
    t(summ_treat_glm(procova_glm$fit_ancova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "glm",
    model  = "PROCOVA - superlearner",
    t(summ_treat_glm(procova_glm$fit_procova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
)


  # ----- rstanarm PROCOVA - superlearner -----
  procova_stan <- procova_analysis_stan(
    df_historical = df_h,
    df_current = df_c,
    model = prognostic_model_superlearn,
    model_covar_names = model_covar_names,
    seed = seed + 300000 + i
  )

r_stan <- rbind(
  data.frame(
    method = "stan",
    model  = "ANOVA",
    t(summ_treat_stan_brms(procova_stan$fit_anova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "stan",
    model  = "ANCOVA",
    t(summ_treat_stan_brms(procova_stan$fit_ancova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "stan",
    model  = "PROCOVA - superlearner",
    t(summ_treat_stan_brms(procova_stan$fit_procova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
)

  
    # ----- rstanarm PROCOVA - oracle -----
  procova_stan <- procova_analysis_stan(
    df_historical = df_h,
    df_current = df_c,
    model = "oracle",
    model_covar_names = model_covar_names,
    seed = seed + 300000 + i
  )

r_stan_oracle <- data.frame(
    method = "stan",
    model  = "PROCOVA - oracle",
    t(summ_treat_stan_brms(procova_stan$fit_procova)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )



  # ----- comm prior using psborrow2 -----

  comm_db_psborrow2 <-   comm_prior_db_analysis_psborrow2(
    df_historical = df_h,
    df_current = df_c,
    tau_rate = 0.01,
    seed = seed + 400000 + i
  )

r_comm <- rbind(
  data.frame(
    method = "psborrow2",
    model  = "COMM_DB",
    t(summ_treat_comm_psborrow2(comm_db_psborrow2$fit_comm)),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "psborrow2",
    model  = "NO_BORROWING",
    t(summ_treat_comm_psborrow2(comm_db_psborrow2$fit_none)),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "psborrow2",
    model  = "POOLING",
    t(summ_treat_comm_psborrow2(comm_db_psborrow2$fit_full)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
)

   # ----- adaptive lasso -----
  
  
  adaptive_lasso_db <- adaptive_lasso_db_analysis(
    df_historical = df_h,
    df_current = df_c,
    lambda = 2,
    gamma = 1.5)
  
  r_alasso <- data.frame(method="glmnet", model="A_LASSO",   t(summ_treat_adaptive_lasso(adaptive_lasso_db$fit)),       
      row.names = NULL,
      stringsAsFactors = FALSE)
  
  adaptive_lasso_adjusted_db <- adaptive_lasso_db_analysis(
    df_historical = df_h,
    df_current = df_c, 
    model_covar_names = c("X1","X2","X3","X4","X5"), 
    lambda = 2,
    gamma = 1.5)
  
  r_alasso_adj <- data.frame(method="glmnet", model="A_LASSO_ADJ",   t(summ_treat_adaptive_lasso(adaptive_lasso_adjusted_db$fit)),       
      row.names = NULL,
      stringsAsFactors = FALSE)

   # ----- comm prior using psborrow2 -----

  
  samprior_fit<-sam_analysis_rbest(
    df_historical = df_h,
    df_current = df_c,
    y_col = "y",
    trt_col = "XTreat",
    delta_mult = 0.5)


  r_sam <- rbind(
    data.frame(
    method = "RBest",
    model  = "SAM_PRIOR",
    t(summ_treat_sam_rbest(samprior_fit, model="SAM")),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "RBest",
    model  = "NO_BORROWING",
    t(summ_treat_sam_rbest(samprior_fit, model="none")),
    row.names = NULL,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "RBest",
    model  = "FULL_BORROWING",
    t(summ_treat_sam_rbest(samprior_fit, model="full")),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  )


  # Combine into one data.frame and add performance components
  tmp <- bind_rows(r_glm, r_stan, r_stan_oracle, r_comm, r_alasso,r_alasso_adj, r_sam)

  tmp$sim <- i

  num_cols <- setdiff(names(tmp), c("method","model","sim"))
  tmp[num_cols] <- lapply(tmp[num_cols], as.numeric)

  tmp$true_tau <- tau
  tmp$bias     <- tmp$est_mean - tau
  tmp$se       <- (tmp$est_mean - tau)^2

  res_list[[i]] <- tmp
}

res <- do.call(rbind, res_list)


# simulate normal data



plot_df <- res %>%
  group_by(method, model) %>%
  summarise(
    bias   = mean(bias),
    MSE    = mean(se),
    reject = mean(reject),
    emp_var = var(est_mean)
  )%>%ungroup()

plot_df <- plot_df %>%
  mutate(
    method = unlist(method),
    model  = factor(unlist(model), levels=c("NO_BORROWING", "ANOVA", "ANCOVA", "PROCOVA - superlearner", "PROCOVA - oracle", "A_LASSO", "A_LASSO_ADJ", "COMM_DB", "SAM_PRIOR", "FULL_BORROWING", "POOLING")),
  )

ggplot(plot_df,
       aes(x = model, y = bias,
           colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Bias", x = NULL)

ggplot(plot_df,
       aes(x = model, y = MSE,
           colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "MSE", x = NULL)

ggplot(plot_df,
       aes(x = model, y = reject,
           colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Rej", x = NULL)

ggplot(plot_df,
       aes(x = model, y = emp_var,
           colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Empirical Variance", x = NULL)

p_bias <- ggplot(plot_df,
                 aes(x = model, y = bias,
                     colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Bias", x = NULL)

p_mse <- ggplot(plot_df,
                aes(x = model, y = MSE,
                    colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "MSE", x = NULL)

p_rej <- ggplot(plot_df,
                aes(x = model, y = reject,
                    colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Rej", x = NULL)

p_empvar <- ggplot(plot_df,
       aes(x = model, y = emp_var,
           colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Empirical Variance", x = NULL)

ggsave("results/test_sim/bias_plot.png", p_bias, width = 12, height = 5, dpi = 300)
ggsave("results/test_sim/mse_plot.png",  p_mse,  width = 12, height = 5, dpi = 300)
ggsave("results/test_sim/rej_plot.png",  p_rej,  width = 12, height = 5, dpi = 300)
ggsave("results/test_sim/empvar_plot.png",  p_empvar,  width = 12, height = 5, dpi = 300)

save(
  res, plot_df,
  p_bias, p_mse, p_rej,
  file = "results/test_sim/simulation_results.RData"
)
