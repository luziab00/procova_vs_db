#### load all necessary packages
# Package names
packages <- c(
"rlang",
  "ggplot2",
  "survival",
  "nph",
  "gtsummary",
  "survminer",
  "ggsurvfit",
  "parallel",    
  "MASS",
  "Matrix",
  "pbmcapply",
  "ggpubr",
  "tidyverse",
  "glmnet",
  "rpact",
  "purrr",
  "future",
  "furrr",
  "tictoc",
  "numbers",
  "future.apply",
  "data.table",
  "patchwork",
  "progressr", 
  "rstanarm",
  "SAMprior",
  "SuperLearner",
  "ranger",
  "xgboost",
  "psborrow2", 
  "cmdstanr",
  "RBesT"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

check_cmdstan_toolchain()