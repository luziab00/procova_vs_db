source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")
source("analysis/procova_analysis.R")
source("analysis/comm_prior.R")
source("analysis/SAM_prior.R")
source("analysis/adaptive_lasso.R")
source("analysis/bayesian_procova.R")

# Function to parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("No input parameters file provided")
  }
  return(args[1])
}


# Source input parameters from the specified file
#input_file <- parse_args()
input_file <- "config/input_baseline.R"
source(input_file)



# Define output directory
output_dir <- paste0("results/", str_remove(str_remove(input_file, "config/input_"), "\\.R"))


# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
#
file.copy(input_file,
          paste0("results/", str_remove(
            str_remove(input_file, "config/input_"), "\\.R"
          ), "/input_file.R"),
          overwrite = TRUE)


cat("Results are saved in:", output_dir, "\n")

set.seed(seed)

sims <- vector("list", nsim)



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
    N = 200,
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

print(paste0("Oracle R-squared across simulations: mean=", round(mean(r_sqared_list),3),
             ", sd=", round(sd(r_sqared_list),3)))



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
    oracle = TRUE,
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
    lambda = 0.005,
    gamma = 1.5)
  
  r_alasso <- data.frame(method="glmnet", model="A_LASSO",   t(summ_treat_adaptive_lasso(adaptive_lasso_db)),       
      row.names = NULL,
      stringsAsFactors = FALSE)
  
  adaptive_lasso_adjusted_db <- adaptive_lasso_db_analysis(
    df_historical = df_h,
    df_current = df_c, 
    model_covar_names = c("X1","X2","X3","X4","X5"), 
    lambda = 0.005,
    gamma = 1.5)
  
  r_alasso_adj <- data.frame(method="glmnet", model="A_LASSO_ADJ",   t(summ_treat_adaptive_lasso(adaptive_lasso_adjusted_db)),       
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

      
  # ----- Bayesian PROCOVA - superlearner -----
  bayesian_procova_fit <- bayesian_procova_analysis(
    df_historical = df_h,
    df_current = df_c,
    model = prognostic_model_superlearn,
    seed = seed + 300000 + i
  )

  r_bayesian_procova <- data.frame(
    method = "gibbs",
    model  = "Bayesian_PROCOVA",
    t(summ_treat_bayesian_procova(bayesian_procova_fit)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )



  # Combine into one data.frame and add performance components
  tmp <- bind_rows(r_glm, r_stan, r_stan_oracle, r_comm, r_alasso,r_alasso_adj, r_sam, r_bayesian_procova)

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
    emp_var = var(est_mean),
    mc_sd = sd(est_mean)
  )%>%ungroup()

plot_df <- plot_df %>%
  mutate(
    method = unlist(method),
    model  = factor(unlist(model), levels=c("NO_BORROWING", "ANOVA", "ANCOVA", "PROCOVA - superlearner", "PROCOVA - oracle", "A_LASSO", "A_LASSO_ADJ", "COMM_DB", "SAM_PRIOR", "FULL_BORROWING", "POOLING", "Bayesian_PROCOVA")),
  )

ggplot(plot_df,
       aes(x = model, y = bias, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Bias", x = NULL)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(plot_df,
       aes(x = model, y = MSE,
         shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "MSE", x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(plot_df,
       aes(x = model, y = reject, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Rej", x = NULL) +
  geom_hline(yintercept = 0.05, linewidth=0.5, color="red", linetype="dashed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(plot_df,
       aes(x = model, y = emp_var, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Empirical Variance", x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_bias <- ggplot(plot_df,
                 aes(x = model, y = bias,
                     colour = model, shape = method)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(y = "Bias", x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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

ggsave(paste(output_dir, "bias_plot.png", sep = "/"), p_bias, width = 12, height = 5, dpi = 300)
ggsave(paste(output_dir, "mse_plot.png", sep = "/"),  p_mse,  width = 12, height = 5, dpi = 300)
ggsave(paste(output_dir, "rej_plot.png", sep = "/"),  p_rej,  width = 12, height = 5, dpi = 300)
ggsave(paste(output_dir, "empvar_plot.png", sep = "/"),  p_empvar,  width = 12, height = 5, dpi = 300)

save(
  res, plot_df,
  file = paste(output_dir, "simulation_results.RData", sep = "/")
)
