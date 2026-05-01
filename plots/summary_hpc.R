# ==============================================================================
# aggregate.R  -  run after all SLURM jobs complete
#
# Usage: Rscript aggregate.R config/input_baseline.R
# ==============================================================================
library(tidyverse)

args        <- commandArgs(trailingOnly = TRUE)
input_file  <- if (length(args) >= 1) args[1] else "config/input_baseline.R"

# Source config to get scenario_grid (needed for varying col detection)
source(input_file)

config_name <- str_remove(str_remove(basename(input_file), "^input_"), "\\.R$")
output_dir  <- file.path("results", config_name)

# ---- Collect per-iteration RDS files ----------------------------------------
# Structure on disk: results/<config>/scenario_001/iter_0001.rds
rds_files <- list.files(output_dir, pattern = "^iter_.*\\.rds$",
                        full.names = TRUE, recursive = TRUE)

expected <- nrow(scenario_grid) * nsim
cat(sprintf("Found %d / %d iteration files\n", length(rds_files), expected))

if (length(rds_files) == 0) {
  stop("No iteration files found under: ", output_dir,
       "\nExpected path pattern: results/<config>/scenario_XXX/iter_XXXX.rds")
}

res <- dplyr::bind_rows(lapply(rds_files, readRDS))
cat(sprintf("Total rows: %d\n", nrow(res)))

# ---- Detect varying columns --------------------------------------------------
# Any scenario_grid column (excluding scenario_id, score_std, treatment_effect)
# with >1 unique value. treatment_effect is excluded here because it is handled
# explicitly via plot_df_t1e / plot_df_power splits rather than faceting.
meta_cols    <- setdiff(names(scenario_grid), c("score_std", "scenario_id"))
varying_cols <- setdiff(
  meta_cols[sapply(meta_cols, function(v) dplyr::n_distinct(res[[v]]) > 1)],
  "treatment_effect"
)
cat("Varying parameters (excl. treatment_effect):", paste(varying_cols, collapse=", "), "\n")

# ---- Summarise ---------------------------------------------------------------
model_levels <- c("NO_BORROWING","ANOVA","ANCOVA",
                  "PROCOVA - superlearner","PROCOVA - oracle",
                  "A_LASSO","A_LASSO_ADJ",
                  "COMM_DB","SAM_PRIOR","FULL_BORROWING","POOLING",
                  "Bayesian_PROCOVA")

group_vars <- c("scenario_id", meta_cols, "method", "model")

plot_df <- res %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    bias         = mean(bias),
    MSE          = mean(sq_error),
    reject       = mean(reject),
    mc_se_reject = sqrt(mean(reject) * (1 - mean(reject)) / n()),
    emp_var      = var(est_mean),
    mc_sd        = sd(est_mean),
    n_sims       = n(),
    .groups = "drop"
  ) %>%
  mutate(model = factor(model, levels = model_levels))

# Split by treatment effect for separate T1E / power plots
plot_df_t1e   <- filter(plot_df, treatment_effect == 0)
plot_df_power <- filter(plot_df, treatment_effect == 0.35)

summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

save(res, plot_df, file = file.path(summary_dir, "simulation_results.RData"))
write_csv(plot_df,  file.path(summary_dir, "simulation_summary.csv"))

# ---- Plot helpers ------------------------------------------------------------
theme_set(theme_bw())
x_theme <- theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

facet_dot_plot <- function(df, y_var, y_lab, extra=NULL) {
  p <- ggplot(df, aes(x=model, y=.data[[y_var]], colour=model, shape=method)) +
    geom_point(size=2) + labs(y=y_lab, x=NULL) + x_theme
  if (y_var == "reject") {
    p <- p + geom_errorbar(
      aes(ymin = reject - 1.96 * mc_se_reject,
          ymax = reject + 1.96 * mc_se_reject),
      width = 0.3
    )
  }
  if (length(varying_cols) > 0) {
    rhs <- varying_cols[length(varying_cols)]
    lhs <- if (length(varying_cols) > 1)
              paste(varying_cols[-length(varying_cols)], collapse=" + ") else "."
    p <- p + facet_grid(as.formula(paste(lhs, "~", rhs)), labeller=label_both)
  }
  if (!is.null(extra)) p <- p + extra
  p
}

# Save a T1E and power version of a dot plot
save_pair <- function(y_var, y_lab, file_stem, extra=NULL, width=16, height=8) {
  ggsave(file.path(summary_dir, paste0(file_stem, "_t1e.png")),
         facet_dot_plot(plot_df_t1e,   y_var, paste(y_lab, "(T1E, trt=0)")),
         width=width, height=height, dpi=300)
  ggsave(file.path(summary_dir, paste0(file_stem, "_power.png")),
         facet_dot_plot(plot_df_power, y_var, paste(y_lab, "(Power, trt=0.35)")),
         width=width, height=height, dpi=300)
}

# ---- Dot plots ---------------------------------------------------------------
save_pair("bias",    "Bias",               "bias_plot")
save_pair("MSE",     "MSE",                "mse_plot")
save_pair("reject",  "Rejection rate",     "rej_plot",
          geom_hline(yintercept=0.05, linewidth=0.5, colour="red", linetype="dashed"))
save_pair("emp_var", "Empirical Variance", "empvar_plot")

# ---- Line plots per varying parameter ----------------------------------------
for (param in varying_cols) {

  make_lp <- function(df, y_var, y_lab, extra=NULL) {
    p <- ggplot(df,
      aes(x=.data[[param]], y=.data[[y_var]],
          colour=model, shape=method,
          group=interaction(model, method))) +
      geom_line() + geom_point(size=2) +
      labs(y=y_lab, x=param) + theme_bw()
    if (y_var == "reject") {
      p <- p + geom_errorbar(
        aes(ymin = reject - 1.96 * mc_se_reject,
            ymax = reject + 1.96 * mc_se_reject),
        width = 0.02
      )
    }
    if (!is.null(extra)) p <- p + extra
    p
  }

  rej_extra <- geom_hline(yintercept=0.05, linewidth=0.5,
                           colour="red", linetype="dashed")

  for (split in list(list(df=plot_df_t1e,   suffix="t1e",   label="T1E, trt=0"),
                     list(df=plot_df_power, suffix="power", label="Power, trt=0.35"))) {
    ggsave(file.path(summary_dir, paste0("bias_vs_",  param, "_", split$suffix, ".png")),
           make_lp(split$df, "bias",   paste("Bias -",            split$label)),
           width=10, height=5, dpi=300)
    ggsave(file.path(summary_dir, paste0("mse_vs_",   param, "_", split$suffix, ".png")),
           make_lp(split$df, "MSE",    paste("MSE -",             split$label)),
           width=10, height=5, dpi=300)
    ggsave(file.path(summary_dir, paste0("rej_vs_",   param, "_", split$suffix, ".png")),
           make_lp(split$df, "reject", paste("Rejection rate -",  split$label), rej_extra),
           width=10, height=5, dpi=300)
  }
}

cat("Done. Summary saved in:", summary_dir, "\n")