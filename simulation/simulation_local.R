# ==============================================================================
# run_simulation_local.R
# Runs the simulation in parallel on your local machine.
# Mirrors what SLURM does on the HPC but uses all available CPU cores.
#
# Usage: Rscript run_simulation_local.R config/input_baseline.R
#    or just source it interactively in RStudio
# ==============================================================================

library(parallel)

# Disable renv sandbox for parallel workers (avoids ~12s startup delay per worker)
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")

# ---- Config ------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "config/input_mean_shift_variation.R"

source(input_file)
# Defines: scenario_grid, make_scenario_inputs(), nsim, seed

n_scenarios <- nrow(scenario_grid)
total_tasks <- n_scenarios * nsim

# Number of cores to use — leave 1 free so your Mac stays responsive
n_cores <- max(1, detectCores() - 1)

# How often to print a progress update to the terminal
print_every <- 100

# Capture working directory now — workers won't inherit it automatically
main_wd <- getwd()

# ---- Shared progress file & log directory ------------------------------------
config_name   <- sub("\\.R$", "", sub("^input_", "", basename(input_file)))
log_dir       <- file.path("results", config_name, "logs")
progress_file <- file.path(log_dir, ".progress")   # one line appended per completed task
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(character(0), progress_file)             # reset/create on each run

cat(sprintf("[%s] Starting %d tasks on %d cores (reporting every %d completions)\n",
            format(Sys.time(), "%H:%M:%S"), total_tasks, n_cores, print_every))

# ---- Worker function ---------------------------------------------------------
run_task <- function(task_id, input_file, main_wd) {

  # Workers don't inherit the working directory — set it explicitly first
  # so all relative paths (source(), file.path(), etc.) resolve correctly
  setwd(main_wd)

  # Disable renv sandbox in workers too
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")

  library(tidyverse)
  source("utils/__init__.R")
  source("utils/helper_functions.R")
  source("data_generating/trial_generation.R")
  source("analysis/procova_analysis.R")
  source("analysis/comm_prior.R")
  source("analysis/SAM_prior.R")
  source("analysis/adaptive_lasso.R")
  source("analysis/bayesian_procova.R")

  source(input_file)

  scenario_idx <- ceiling(task_id / nsim)
  iter         <- task_id - (scenario_idx - 1L) * nsim

  scenario_row <- scenario_grid[scenario_idx, , drop = FALSE]
  inputs       <- make_scenario_inputs(scenario_row)
  score_std    <- scenario_row$score_std[[1]]
  sim_seed_offset <- 10000L * scenario_idx + iter

  # ---- Log file & progress helpers -------------------------------------------
  config_name   <- sub("\\.R$", "", sub("^input_", "", basename(input_file)))
  log_dir       <- file.path("results", config_name, "logs")
  log_file      <- file.path(log_dir, sprintf("task_%05d_scenario_%03d_iter_%04d.log",
                                              task_id, scenario_row$scenario_id, iter))
  progress_file <- file.path(log_dir, ".progress")

  log_msg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n")
    cat(msg, file = log_file, append = TRUE)
  }

  # Append one line to the shared progress file when a task finishes.
  # Uses a simple lock file to avoid garbled writes from simultaneous workers.
  report_progress <- function(status) {
    line <- sprintf("%d\t%d\t%d\t%s\t%s\n",
                    task_id, scenario_idx, iter, status,
                    format(Sys.time(), "%H:%M:%S"))
    lock <- paste0(progress_file, ".lock")
    while (file.exists(lock)) Sys.sleep(0.01)
    file.create(lock, showWarnings = FALSE)
    cat(line, file = progress_file, append = TRUE)
    file.remove(lock)
  }

  log_msg("START task %d | scenario %d | iter %d", task_id, scenario_idx, iter)

  # ---- Safe per-method wrapper -----------------------------------------------
  # Runs expr, logs warnings without aborting, catches errors and returns NULL
  run_method <- function(method_name, expr) {
    withCallingHandlers(
      tryCatch({
        log_msg("Running %s", method_name)
        expr
      }, error = function(e) {
        log_msg("ERROR in %s: %s", method_name, conditionMessage(e))
        NULL
      }),
      warning = function(w) {
        log_msg("WARNING in %s: %s", method_name, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }

  # ---- Data generation (always saved, even if analyses fail) -----------------
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

  # Save datasets immediately before any analysis that could fail
  output_dir <- file.path("results", config_name,
                          sprintf("scenario_%03d", scenario_row$scenario_id))
  data_dir   <- file.path(output_dir, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(df_h, file.path(data_dir, sprintf("df_h_iter_%04d.rds", iter)))
  saveRDS(df_c, file.path(data_dir, sprintf("df_c_iter_%04d.rds", iter)))
  log_msg("Datasets saved to %s", data_dir)

  # ==========================================================================
  # ANALYSES — add / remove methods here; nothing else needs to change
  # ==========================================================================

  # -- Frequentist PROCOVA -----------------------------------------------------
  r_glm <- run_method("frequentist PROCOVA", {
    procova_glm <- procova_analysis_glm(
      df_historical     = df_h, df_current = df_c,
      model             = prognostic_model_superlearn,
      model_covar_names = inputs$model_covar_names,
      seed              = seed + 300000L + sim_seed_offset
    )
    rbind(
      data.frame(method="glm", model="ANOVA",
                 t(summ_treat_glm(procova_glm$fit_anova)),   row.names=NULL, stringsAsFactors=FALSE),
      data.frame(method="glm", model="ANCOVA",
                 t(summ_treat_glm(procova_glm$fit_ancova)),  row.names=NULL, stringsAsFactors=FALSE),
      data.frame(method="glm", model="PROCOVA - superlearner",
                 t(summ_treat_glm(procova_glm$fit_procova)), row.names=NULL, stringsAsFactors=FALSE)
    )
  })

  # -- Bayesian PROCOVA (stan) – superlearner ----------------------------------
  #r_stan <- run_method("Bayesian PROCOVA (stan) - superlearner", {
  #  procova_stan <- procova_analysis_stan(
  #    df_historical     = df_h, df_current = df_c,
  #    model             = prognostic_model_superlearn,
  #    model_covar_names = inputs$model_covar_names,
  #    seed              = seed + 300000L + sim_seed_offset
  #  )
  #  rbind(
  #    data.frame(method="stan", model="ANOVA",
  #               t(summ_treat_stan_brms(procova_stan$fit_anova)),   row.names=NULL, stringsAsFactors=FALSE),
  #    data.frame(method="stan", model="ANCOVA",
  #               t(summ_treat_stan_brms(procova_stan$fit_ancova)),  row.names=NULL, stringsAsFactors=FALSE),
  #    data.frame(method="stan", model="PROCOVA - superlearner",
  #               t(summ_treat_stan_brms(procova_stan$fit_procova)), row.names=NULL, stringsAsFactors=FALSE)
  #  )
  #})

  # -- Bayesian PROCOVA (stan) – oracle ----------------------------------------
  #r_stan_oracle <- run_method("Bayesian PROCOVA (stan) - oracle", {
  #  procova_stan_oracle <- procova_analysis_stan(
  #    df_historical     = df_h, df_current = df_c,
  #    oracle            = TRUE,
  #    model_covar_names = inputs$model_covar_names,
  #    seed              = seed + 300000L + sim_seed_offset
  #  )
  #  data.frame(
  #    method="stan", model="PROCOVA - oracle",
  #    t(summ_treat_stan_brms(procova_stan_oracle$fit_procova)),
  #    row.names=NULL, stringsAsFactors=FALSE
  #  )
  #})

  # -- Commensurate prior (psborrow2) ------------------------------------------
  #r_comm <- run_method("commensurate prior (psborrow2)", {
  #  comm_db <- comm_prior_db_analysis_psborrow2(
  #    df_historical = df_h, df_current = df_c,
  #    tau_rate = inputs$tau_rate,
  #    seed     = seed + 400000L + sim_seed_offset
  #  )
  #  rbind(
  #    data.frame(method="psborrow2", model="COMM_DB",
  #               t(summ_treat_comm_psborrow2(comm_db$fit_comm)), row.names=NULL, stringsAsFactors=FALSE),
  #    data.frame(method="psborrow2", model="POOLING",
  #               t(summ_treat_comm_psborrow2(comm_db$fit_full)), row.names=NULL, stringsAsFactors=FALSE)
  #  )
  #})

  # -- Adaptive lasso ----------------------------------------------------------
  r_alasso <- run_method("adaptive lasso", {
    rbind(
      data.frame(method="glmnet", model="A_LASSO",
                 t(summ_treat_adaptive_lasso(
                   adaptive_lasso_db_analysis(df_h, df_c,
                                              lambda = inputs$lambda_alasso,
                                              gamma  = inputs$gamma_alasso)
                 )), row.names=NULL, stringsAsFactors=FALSE),
      data.frame(method="glmnet", model="A_LASSO_ADJ",
                 t(summ_treat_adaptive_lasso(
                   adaptive_lasso_db_analysis(df_h, df_c,
                                              model_covar_names = inputs$model_covar_names,
                                              lambda = inputs$lambda_alasso,
                                              gamma  = inputs$gamma_alasso)
                 )), row.names=NULL, stringsAsFactors=FALSE)
    )
  })

  # -- SAM prior (RBest) -------------------------------------------------------
  #r_sam <- run_method("SAM prior (RBest)", {
  #  samprior_fit <- sam_analysis_rbest(
  #    df_historical = df_h, df_current = df_c,
  #    y_col="y", trt_col="XTreat",
  #    delta_mult = inputs$delta_mult_sam
  #  )
  #  rbind(
  #    data.frame(method="RBest", model="SAM_PRIOR",
  #               t(summ_treat_sam_rbest(samprior_fit, model="SAM")),  row.names=NULL, stringsAsFactors=FALSE),
  #    data.frame(method="RBest", model="NO_BORROWING",
  #               t(summ_treat_sam_rbest(samprior_fit, model="none")), row.names=NULL, stringsAsFactors=FALSE),
 #     data.frame(method="RBest", model="FULL_BORROWING",
 #                t(summ_treat_sam_rbest(samprior_fit, model="full")), row.names=NULL, stringsAsFactors=FALSE)
 #   )
 # })

  # -- Bayesian PROCOVA (Gibbs) ------------------------------------------------
  r_bp <- run_method("Bayesian PROCOVA (Gibbs)", {
    bp_fit <- bayesian_procova_analysis(
      df_historical = df_h, df_current = df_c,
      model = prognostic_model_superlearn,
      seed  = seed + 300000L + sim_seed_offset
    )
    data.frame(
      method="gibbs", model="Bayesian_PROCOVA",
      t(summ_treat_bayesian_procova(bp_fit)),
      row.names=NULL, stringsAsFactors=FALSE
    )
  })

    # -- Bayesian PROCOVA (Gibbs) ------------------------------------------------
  r_bp_conservative <- run_method("Bayesian PROCOVA (Gibbs)", {
    bp_fit <- bayesian_procova_analysis(
      df_historical = df_h, df_current = df_c,
      model = prognostic_model_superlearn,
      seed  = seed + 300000L + sim_seed_offset,
      K0H=1/sqrt(200)
    )
    data.frame(
      method="gibbs", model="Bayesian_PROCOVA_conservative",
      t(summ_treat_bayesian_procova(bp_fit)),
      row.names=NULL, stringsAsFactors=FALSE
    )
  })


  # ---- Assemble (drop NULLs from any failed methods) -------------------------
  results_list <- list(r_glm, r_alasso, r_bp, r_bp_conservative) #list(r_glm, r_stan, r_stan_oracle, r_comm, r_alasso, r_sam, r_bp, r_bp_conservative)

  n_failed     <- sum(sapply(results_list, is.null))
  n_ok         <- length(results_list) - n_failed

  if (n_failed > 0) log_msg("%d/%d method(s) failed and will be absent from results",
                             n_failed, length(results_list))

  res <- dplyr::bind_rows(Filter(Negate(is.null), results_list))

  # ---- Save & report ---------------------------------------------------------
  if (nrow(res) == 0) {
    log_msg("All methods failed — no results to save for task %d", task_id)
    report_progress("FAIL")
    cat(sprintf("[task %d] scenario %d | iter %d FAILED — see log\n",
                task_id, scenario_idx, iter))
  } else {
    num_cols      <- setdiff(names(res), c("method", "model"))
    res[num_cols] <- lapply(res[num_cols], function(x) suppressWarnings(as.numeric(x)))

    res$sim      <- iter
    res$true_tau <- inputs$treatment_effect
    res$bias     <- res$est_mean - inputs$treatment_effect
    res$sq_error <- (res$est_mean - inputs$treatment_effect)^2

    scenario_cols <- setdiff(names(scenario_row), "score_std")
    res <- dplyr::bind_cols(res, scenario_row[rep(1, nrow(res)), scenario_cols, drop=FALSE])

    saveRDS(res, file.path(output_dir, sprintf("iter_%04d.rds", iter)))
    log_msg("SUCCESS — saved %d/%d methods -> iter_%04d.rds", n_ok, length(results_list), iter)
    report_progress(if (n_failed > 0) "PARTIAL" else "OK")
    cat(sprintf("[task %d] scenario %d | iter %d done (%d/%d methods OK)\n",
                task_id, scenario_idx, iter, n_ok, length(results_list)))
  }

  log_msg("END task %d", task_id)
  invisible(NULL)
}

# ---- Run in parallel, printing progress every `print_every` tasks ------------
cl         <- makeCluster(n_cores)
start_time <- Sys.time()

clusterExport(cl, varlist = c("input_file", "main_wd"), envir = environment())

# Build non-overlapping chunks of size print_every
task_chunks <- lapply(
  seq(1, total_tasks, by = print_every),
  function(start) start:min(start + print_every - 1, total_tasks)
)

tryCatch({
  for (chunk in task_chunks) {
    parLapply(cl, chunk, run_task, input_file = input_file, main_wd = main_wd)

    # Read shared progress file and print update
    progress_lines <- readLines(progress_file)
    n_done    <- length(progress_lines)
    n_failed  <- sum(grepl("\tFAIL\t",    progress_lines))
    n_partial <- sum(grepl("\tPARTIAL\t", progress_lines))
    n_ok      <- n_done - n_failed - n_partial
    elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate      <- if (elapsed > 0) n_done / elapsed else NA
    eta_min   <- if (!is.na(rate) && rate > 0) (total_tasks - n_done) / rate else NA
    pct       <- round(100 * n_done / total_tasks)

    cat(sprintf("[%s] %d/%d (%d%%) | %d OK | %d partial | %d failed | %.1f min elapsed | ETA: %s\n",
                format(Sys.time(), "%H:%M:%S"),
                n_done, total_tasks, pct, n_ok, n_partial, n_failed,
                elapsed,
                if (!is.na(eta_min)) sprintf("%.1f min", eta_min) else "?"))
  }
}, finally = {
  stopCluster(cl)
})

# ---- Final summary -----------------------------------------------------------
elapsed_total  <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
progress_lines <- readLines(progress_file)
n_failed       <- sum(grepl("\tFAIL\t",    progress_lines))
n_partial      <- sum(grepl("\tPARTIAL\t", progress_lines))

cat(sprintf("\nFinished %d tasks in %.1f minutes.\n", total_tasks, elapsed_total))

if (n_failed == 0 && n_partial == 0) {
  cat("All tasks completed successfully.\n")
} else {
  if (n_failed > 0) {
    cat(sprintf("%d task(s) FAILED completely. Check these logs:\n", n_failed))
    failed_logs <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)
    failed_logs <- failed_logs[sapply(failed_logs, function(f) any(grepl("FAILED", readLines(f))))]
    cat(paste(" -", failed_logs, collapse = "\n"), "\n")
  }
  if (n_partial > 0) {
    cat(sprintf("%d task(s) completed with some methods FAILING. Check these logs:\n", n_partial))
    partial_logs <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)
    partial_logs <- partial_logs[sapply(partial_logs, function(f) any(grepl("ERROR in", readLines(f))))]
    cat(paste(" -", partial_logs, collapse = "\n"), "\n")
  }
}

cat("Run aggregate.R to collect results.\n")