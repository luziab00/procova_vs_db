# ==============================================================================
# config/input_baseline.R
#
# This file defines EVERYTHING needed for the simulation study:
#   1. Fixed settings
#   2. The scenario grid (one row per scenario)
#   3. make_scenario_inputs() — converts a grid row into full parameter lists
#   4. Pre-computed score standardisation per scenario (stored in the grid)
#
# Sourced once by run_simulation_hpc.R; the script stays fully generic.
# ==============================================================================

# ---- Dependencies (needed for grid construction and score_std) ---------------
library(tidyverse)
source("utils/__init__.R")
source("utils/helper_functions.R")
source("data_generating/trial_generation.R")

# ==============================================================================
# 1. Fixed settings
# ==============================================================================
nsim     <- 20
seed     <- 123

N_hist_default    <- 200
N_current_default <- 200
sigma_total       <- 1
#treatment_effect  <- 0
R2_target         <- 0.1
pi_value          <- sqrt(R2_target) * sigma_total

tau_rate      <- 0.1
lambda_alasso <- 0.001
gamma_alasso  <- 1.5
delta_mult_sam <- 0.25

model_covar_names <- c("X1", "X2", "X3", "X4", "X5")

# Base parameter objects — make_scenario_inputs() modifies copies of these
coeff_hist_base <- list(
  beta0  = 0,
  beta1  = 0.6, beta2  = 0.3, beta3  = 0.2, beta4  = 0.4,
  gamma1 = 0.5, gamma2 = 0.4, gamma3 = 0.6, gamma4 = 0.3, gamma5 = 0.5
)

coeff_current_base <- list(
  beta0  = 0,
  beta1  = 0.6, beta2  = 0.3, beta3  = 0.2, beta4  = 0.4,
  gamma1 = 0.5, gamma2 = 0.4, gamma3 = 0.6, gamma4 = 0.3, gamma5 = 0.5
)

covar_distr_hist_base <- list(
  mu  = c(X1=0, X2=0, X3=0, X4=0),
  sd  = c(X1=1, X2=1, X3=1, X4=1),
  pX5 = 0.5
)

covar_distr_current_base <- list(
  mu  = c(X1=0, X2=0, X3=0, X4=0),
  sd  = c(X1=1, X2=1, X3=1, X4=1),
  pX5 = 0.5
)

# ==============================================================================
# 2. Scenario grid
#    One row per scenario. Add crossing() dimensions here to expand the design.
# ==============================================================================
scenario_grid <- tidyr::crossing(
  hist_mean_shift  = seq(-0.5, 0.5, 0.25),
  treatment_effect = c(0, 0.35)
  # add more axes here, e.g.:

  # R2_target       = c(0.1, 0.25),
  # N_hist          = c(100, 200)
) %>%
  dplyr::mutate(scenario_id = dplyr::row_number())

# ==============================================================================
# 3. make_scenario_inputs()
#    Converts one scenario_grid row into a fully-specified named list.
#    The simulation script calls this and unpacks the result — it never
#    touches scenario_grid directly.
# ==============================================================================
make_scenario_inputs <- function(scenario_row) {

  # Start from base objects and apply scenario-specific overrides
  coeff_hist          <- coeff_hist_base
  coeff_current       <- coeff_current_base
  covar_distr_hist    <- covar_distr_hist_base
  covar_distr_current <- covar_distr_current_base

  # --- Apply overrides from the scenario row ----------------------------------
  coeff_hist$beta0 <- coeff_hist_base$beta0 + scenario_row$hist_mean_shift

  # Add more overrides here as your grid grows, e.g.:
  # coeff_current$beta0   <- scenario_row$current_mean_shift
  # covar_distr_hist$mu[] <- scenario_row$hist_mu
  # N_hist                <- scenario_row$N_hist

  list(
    coeff_hist          = coeff_hist,
    coeff_current       = coeff_current,
    covar_distr_hist    = covar_distr_hist,
    covar_distr_current = covar_distr_current,
    N_hist              = N_hist_default,
    N_current           = N_current_default,
    pi_value            = pi_value,
    sigma_total         = sigma_total,
    treatment_effect    = scenario_row$treatment_effect,
    tau_rate            = tau_rate,
    lambda_alasso       = lambda_alasso,
    gamma_alasso        = gamma_alasso,
    delta_mult_sam      = delta_mult_sam,
    model_covar_names   = model_covar_names
  )
}

# ==============================================================================
# 4. Pre-compute score standardisation for every scenario
#    Stored as a list column in scenario_grid so each task can look it up
#    directly without recomputing.
#
#    Note: score_std depends on covar_distr_hist and coeff_hist, both of which
#    can vary per scenario, so we compute it per row. The seed is made
#    scenario-specific but deterministic.
# ==============================================================================
cat("Pre-computing score standardisation for", nrow(scenario_grid), "scenarios...\n")

scenario_grid$score_std <- lapply(seq_len(nrow(scenario_grid)), function(i) {
  inputs <- make_scenario_inputs(scenario_grid[i, , drop = FALSE])
  compute_score_standardization(
    covars_distr = inputs$covar_distr_hist,
    coeff        = inputs$coeff_hist,
    M            = 200000,
    seed         = seed + 9999 + i    # unique per scenario, deterministic
  )
})

cat("Done.\n")