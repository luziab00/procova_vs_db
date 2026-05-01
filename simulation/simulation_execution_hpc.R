#!/usr/bin/env Rscript
# ==============================================================================
# run_simulation_hpc.R
#
# Follows the same structure as the HPC example:
#   - SLURM array task ID = one iteration of one scenario
#   - Each task runs run_one_sim() exactly once and saves one .rds file
#
# Array size = nrow(scenario_grid) * nsim
# Task ID decoding (example: 11 scenarios, 100 iterations each):
#   task   1-100  -> scenario 1, iter 1-100
#   task 101-200  -> scenario 2, iter 1-100
#   ...
#   task 1001-1100 -> scenario 11, iter 1-100
#
# Usage: Rscript run_simulation_hpc.R config/input_baseline.R
# ==============================================================================

library(tidyverse)

source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")
source("analysis/procova_analysis.R")
source("analysis/comm_prior.R")
source("analysis/SAM_prior.R")
source("analysis/adaptive_lasso.R")
source("analysis/bayesian_procova.R")

# ---- Read config -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_simulation_hpc.R <input_file>")
input_file <- args[1]

source(input_file)
# Defines: scenario_grid (with $score_std list column), make_scenario_inputs(),
#          nsim, seed

stopifnot(
  "scenario_grid must be defined"        = exists("scenario_grid"),
  "make_scenario_inputs must be defined" = exists("make_scenario_inputs"),
  "nsim must be defined"                 = exists("nsim"),
  "seed must be defined"                 = exists("seed")
)

# ---- Decode task ID into (scenario, iteration) -------------------------------
task_id     <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
n_scenarios <- nrow(scenario_grid)
total_tasks <- n_scenarios * nsim

if (is.na(task_id) || task_id < 1 || task_id > total_tasks) {
  stop(sprintf("SLURM_ARRAY_TASK_ID=%d out of range [1, %d]", task_id, total_tasks))
}

scenario_idx <- ceiling(task_id / nsim)
iter         <- task_id - (scenario_idx - 1L) * nsim

scenario_row <- scenario_grid[scenario_idx, , drop = FALSE]
inputs       <- make_scenario_inputs(scenario_row)
score_std    <- scenario_row$score_std[[1]]   # pre-computed in config

cat(sprintf("[task %d] scenario %d/%d | iter %d/%d\n",
            task_id, scenario_idx, n_scenarios, iter, nsim))

# ---- Output directory --------------------------------------------------------
config_name <- str_remove(str_remove(basename(input_file), "^input_"), "\\.R$")
output_dir  <- file.path("results", config_name,
                         sprintf("scenario_%03d", scenario_row$scenario_id))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Log file setup ----------------------------------------------------------
log_dir  <- file.path("results", config_name, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, sprintf("task_%05d_scenario_%03d_iter_%04d.log",
                                       task_id, scenario_row$scenario_id, iter))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n")
  cat(msg)                          # also print to SLURM stdout
  cat(msg, file = log_file, append = TRUE)
}

log_msg("START task %d | scenario %d/%d | iter %d/%d",
        task_id, scenario_idx, n_scenarios, iter, nsim)
log_msg("Host: %s | R version: %s", Sys.info()[["nodename"]], R.version$version.string)

# Save scenario metadata once (overwrite races are fine)
file.copy(input_file, file.path(dirname(output_dir), "input_file.R"), overwrite = TRUE)

# ---- Run one simulation (with error catching) --------------------------------
sim_seed_offset <- 10000L * scenario_idx + iter

result <- tryCatch({

  # -- Data generation ---------------------------------------------------------
  log_msg("Generating historical trial data")
  df_h <- simulate_single_trial_standardized(
    N             = inputs$N_hist,
    covars_distr  = inputs$covar_distr_hist,
    coeff         = inputs$coeff_hist,
    sigma         = inputs$sigma_total,
    pi            = inputs$pi_value,
    score_std     = score_std,
    trial_name    = "historical",
    treatment_arm = FALSE,
    tau           = 0,
    randomization_ratio = NULL,
    seed          = seed + 100000L + sim_seed_offset
  )

  log_msg("Generating current trial data")
  df_c <- simulate_single_trial_standardized(
    N             = inputs$N_current,
    covars_distr  = inputs$covar_distr_current,
    coeff         = inputs$coeff_current,
    sigma         = inputs$sigma_total,
    pi            = inputs$pi_value,
    score_std     = score_std,
    trial_name    = "current",
    treatment_arm = TRUE,
    tau           = inputs$treatment_effect,
    randomization_ratio = 0.5,
    seed          = seed + 200000L + sim_seed_offset
  )

  # ==============================================================================
  # ANALYSES — add / remove methods here; nothing else needs to change
  # ==============================================================================

  # -- Frequentist PROCOVA -----------------------------------------------------
  log_msg("Running frequentist PROCOVA")
  procova_glm <- procova_analysis_glm(
    df_historical     = df_h, df_current = df_c,
    model             = prognostic_model_superlearn,
    model_covar_names = inputs$model_covar_names,
    seed              = seed + 300000L + sim_seed_offset
  )
  r_glm <- rbind(
    data.frame(method="glm", model="ANOVA",
               t(summ_treat_glm(procova_glm$fit_anova)),   row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="glm", model="ANCOVA",
               t(summ_treat_glm(procova_glm$fit_ancova)),  row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="glm", model="PROCOVA - superlearner",
               t(summ_treat_glm(procova_glm$fit_procova)), row.names=NULL, stringsAsFactors=FALSE)
  )

  # -- Bayesian PROCOVA (stan) – superlearner ----------------------------------
  log_msg("Running Bayesian PROCOVA (stan) - superlearner")
  procova_stan <- procova_analysis_stan(
    df_historical     = df_h, df_current = df_c,
    model             = prognostic_model_superlearn,
    model_covar_names = inputs$model_covar_names,
    seed              = seed + 300000L + sim_seed_offset
  )
  r_stan <- rbind(
    data.frame(method="stan", model="ANOVA",
               t(summ_treat_stan_brms(procova_stan$fit_anova)),   row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="stan", model="ANCOVA",
               t(summ_treat_stan_brms(procova_stan$fit_ancova)),  row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="stan", model="PROCOVA - superlearner",
               t(summ_treat_stan_brms(procova_stan$fit_procova)), row.names=NULL, stringsAsFactors=FALSE)
  )

  # -- Bayesian PROCOVA (stan) – oracle ----------------------------------------
  log_msg("Running Bayesian PROCOVA (stan) - oracle")
  procova_stan_oracle <- procova_analysis_stan(
    df_historical     = df_h, df_current = df_c,
    oracle            = TRUE,
    model_covar_names = inputs$model_covar_names,
    seed              = seed + 300000L + sim_seed_offset
  )
  r_stan_oracle <- data.frame(
    method="stan", model="PROCOVA - oracle",
    t(summ_treat_stan_brms(procova_stan_oracle$fit_procova)),
    row.names=NULL, stringsAsFactors=FALSE
  )

  # -- Commensurate prior (psborrow2) ------------------------------------------
  log_msg("Running commensurate prior (psborrow2)")
  comm_db <- comm_prior_db_analysis_psborrow2(
    df_historical = df_h, df_current = df_c,
    tau_rate = inputs$tau_rate,
    seed     = seed + 400000L + sim_seed_offset
  )
  r_comm <- rbind(
    data.frame(method="psborrow2", model="COMM_DB",
               t(summ_treat_comm_psborrow2(comm_db$fit_comm)), row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="psborrow2", model="POOLING",
               t(summ_treat_comm_psborrow2(comm_db$fit_full)), row.names=NULL, stringsAsFactors=FALSE)
  )

  # -- Adaptive lasso ----------------------------------------------------------
  log_msg("Running adaptive lasso")
  r_alasso <- data.frame(
    method="glmnet", model="A_LASSO",
    t(summ_treat_adaptive_lasso(
      adaptive_lasso_db_analysis(df_h, df_c,
                                 lambda = inputs$lambda_alasso,
                                 gamma  = inputs$gamma_alasso)
    )), row.names=NULL, stringsAsFactors=FALSE
  )
  r_alasso_adj <- data.frame(
    method="glmnet", model="A_LASSO_ADJ",
    t(summ_treat_adaptive_lasso(
      adaptive_lasso_db_analysis(df_h, df_c,
                                 model_covar_names = inputs$model_covar_names,
                                 lambda = inputs$lambda_alasso,
                                 gamma  = inputs$gamma_alasso)
    )), row.names=NULL, stringsAsFactors=FALSE
  )

  # -- SAM prior (RBest) -------------------------------------------------------
  log_msg("Running SAM prior (RBest)")
  samprior_fit <- sam_analysis_rbest(
    df_historical = df_h, df_current = df_c,
    y_col="y", trt_col="XTreat",
    delta_mult = inputs$delta_mult_sam
  )
  r_sam <- rbind(
    data.frame(method="RBest", model="SAM_PRIOR",
               t(summ_treat_sam_rbest(samprior_fit, model="SAM")),  row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="RBest", model="NO_BORROWING",
               t(summ_treat_sam_rbest(samprior_fit, model="none")), row.names=NULL, stringsAsFactors=FALSE),
    data.frame(method="RBest", model="FULL_BORROWING",
               t(summ_treat_sam_rbest(samprior_fit, model="full")), row.names=NULL, stringsAsFactors=FALSE)
  )

  # -- Bayesian PROCOVA (Gibbs) ------------------------------------------------
  log_msg("Running Bayesian PROCOVA (Gibbs)")
  bp_fit <- bayesian_procova_analysis(
    df_historical = df_h, df_current = df_c,
    model = prognostic_model_superlearn,
    seed  = seed + 300000L + sim_seed_offset
  )
  r_bp <- data.frame(
    method="gibbs", model="Bayesian_PROCOVA",
    t(summ_treat_bayesian_procova(bp_fit)),
    row.names=NULL, stringsAsFactors=FALSE
  )

  # -- Assemble ----------------------------------------------------------------
  res <- dplyr::bind_rows(r_glm, r_stan, r_stan_oracle, r_comm,
                          r_alasso, r_alasso_adj, r_sam, r_bp)

  num_cols      <- setdiff(names(res), c("method", "model"))
  res[num_cols] <- lapply(res[num_cols], function(x) suppressWarnings(as.numeric(x)))

  res$sim      <- iter
  res$true_tau <- inputs$treatment_effect
  res$bias     <- res$est_mean - inputs$treatment_effect
  res$sq_error <- (res$est_mean - inputs$treatment_effect)^2

  scenario_cols <- setdiff(names(scenario_row), "score_std")
  res <- dplyr::bind_cols(res, scenario_row[rep(1, nrow(res)), scenario_cols, drop=FALSE])

  res  # return value of tryCatch on success

}, warning = function(w) {
  log_msg("WARNING: %s", conditionMessage(w))
  invokeRestart("muffleWarning")   # log but don't abort

}, error = function(e) {
  log_msg("ERROR: %s", conditionMessage(e))
  log_msg("Traceback:\n%s", paste(capture.output(traceback()), collapse = "\n"))
  log_msg("Task FAILED — no .rds saved for task %d", task_id)
  NULL  # return NULL so the script exits cleanly without crashing the whole array
})

# ---- Save (only if simulation succeeded) -------------------------------------
if (!is.null(result)) {
  out_file <- file.path(output_dir, sprintf("iter_%04d.rds", iter))
  saveRDS(result, out_file)
  log_msg("SUCCESS — saved -> %s", out_file)
} else {
  log_msg("SKIPPED saving .rds due to earlier error")
}

log_msg("END task %d", task_id)