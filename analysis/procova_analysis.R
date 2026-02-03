
# generate prognostic models

prognostic_model_superlearn<-function(
  df_historical,
  covar_names   = c("X1","X2","X3","X4", "X5"),  
  SL.library    = c("SL.glm",
                    "SL.glmnet",
                    "SL.ranger",
                    "SL.xgboost"),
  seed          = NULL
){
  if (!is.null(seed)) set.seed(seed)
  
  y <- df_historical[["y"]]
  X <- df_historical[, covar_names, drop = FALSE]
  
  fit_SL <- SuperLearner(
    Y          = y,
    X          = X,
    SL.library = SL.library,
  family        = gaussian(),      
  method        = "method.NNLS"     # non-negative least squares ensemble
  )

  fit_SL
}


prognostic_model_lm <- function(
  df_historical,
  covar_names = c("X1","X2","X3","X4"),
  seed        = NULL
) {
  if (!is.null(seed)) set.seed(seed) 

  fml <- as.formula(
    paste("y ~", paste(covar_names, collapse = " + "))
  )
  
  fit_lm <- lm(fml, data = df_historical)

  fit_lm
}

prognostic_model_oracle_lm <- function(
  df_historical,
  covar_names  = NULL,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Oracle feature construction matching your DGM:
  # linear terms: X1..X4
  # nonlinear terms: X1^2, sin(X3), I(X4>1), X1*X5
  fml <- y ~ X1 + X2 + X3 + X4 +
    I(X1^2) + I(sin(X3)) + I(X4 > 1) + I(X1 * X5)

  fit_oracle <- lm(fml, data = df_historical)

  fit_oracle
}


# procova analysis frequentist
procova_analysis_glm <- function(
  df_historical, 
  df_current, 
  model, 
  model_covar_names = NULL, 
  seed = NULL
){


  if (model=="oracle"){
    df_current$progn_score <- df_current$S
    prog_fit <- lm(y ~ S, data = df_historical)

  } else {

    prog_fit <- model(
      df_historical = df_historical,
      covar_names   = model_covar_names,
      seed          = seed
    )
  
    if (inherits(prog_fit, "SuperLearner")) {
      X_new <- df_current[, model_covar_names, drop = FALSE]
      progn_score <- as.numeric(predict(prog_fit, newdata = X_new)$pred)
    } else {
      progn_score <- as.numeric(predict(prog_fit, newdata = df_current))
    }
  
    df_current$progn_score <- progn_score
  }



  ## --- Build formulas ---
  # ANOVA: treatment only
  f_anova <- as.formula(
    paste0("y", " ~ ", "XTreat")
  )
  
  # ANCOVA: treatment + raw baseline covariates
  f_ancova <- as.formula(
    paste0(
      "y", " ~ ", "XTreat", " + ",
      paste(model_covar_names, collapse = " + ")
    )
  )
  
  # PROCOVA: ANCOVA + prognostic score
  f_procova <- as.formula(
    paste0(
      "y", " ~ ", "XTreat", " + progn_score + ",
      paste(model_covar_names, collapse = " + ")
    )
  )
  
  ## --- Fit models ---
  fit_anova <- lm(f_anova, data = df_current) 
  fit_ancova <- lm(f_ancova, data = df_current) 
  fit_procova <- lm(f_procova, data = df_current)
  
  
  list(
    prognostic_model   = prog_fit,
    current_with_score = df_current,
    fit_anova          = fit_anova,
    fit_ancova         = fit_ancova,
    fit_procova        = fit_procova
  )
}



# procova analysis using stand via brms
procova_analysis_stan <- function(
  df_historical, 
  df_current, 
  model, 
  model_covar_names = NULL, 
  seed = NULL
){


  if (model=="oracle"){
    df_current$progn_score <- df_current$S
    prog_fit <- lm(y ~ S, data = df_historical)

  } else {

    prog_fit <- model(
      df_historical = df_historical,
      covar_names   = model_covar_names,
      seed          = seed
    )
  
    if (inherits(prog_fit, "SuperLearner")) {
      X_new <- df_current[, model_covar_names, drop = FALSE]
      progn_score <- as.numeric(predict(prog_fit, newdata = X_new)$pred)
    } else {
      progn_score <- as.numeric(predict(prog_fit, newdata = df_current))
    }
  
    df_current$progn_score <- progn_score
  }




  ## --- Build formulas ---
  # ANOVA: treatment only
  f_anova <- as.formula(
    paste0("y", " ~ ", "XTreat")
  )
  
  # ANCOVA: treatment + raw baseline covariates
  f_ancova <- as.formula(
    paste0(
      "y", " ~ ", "XTreat", " + ",
      paste(model_covar_names, collapse = " + ")
    )
  )
  
  # PROCOVA: ANCOVA + prognostic score
  f_procova <- as.formula(
    paste0(
      "y", " ~ ", "XTreat", " + progn_score + ",
      paste(model_covar_names, collapse = " + ")
    )
  )
  
  ## --- Fit models ---
  fit_anova  <- stan_glm(
  f_anova,
  data   = df_current,
  family = gaussian(),

  # default flat priors on regression coefficients
  prior           = NULL,
  prior_intercept = NULL,

  # default variance prior 
  prior_aux = NULL,

  # using this so we get as close to the typical ml estimates
  algorithm = "optimizing"
)

  fit_ancova  <-  stan_glm(
  f_ancova,
  data   = df_current,
  family = gaussian(),

  # default flat priors on regression coefficients
  prior           = NULL,
  prior_intercept = NULL,

  # default variance prior 
  prior_aux = NULL,

  # using this so we get as close to the typical ml estimates
  algorithm = "optimizing"

)
  
  fit_procova <- stan_glm(
  f_procova,
  data   = df_current,
  family = gaussian(),

  # default flat priors on regression coefficients
  prior           = NULL,
  prior_intercept = NULL,

  # default variance prior 
  prior_aux = NULL,

  # using this so we get as close to the typical ml estimates
  algorithm = "optimizing"
)
  
  
  list(
    prognostic_model   = prog_fit,
    current_with_score = df_current,
    fit_anova          = fit_anova,
    fit_ancova         = fit_ancova,
    fit_procova        = fit_procova
  )
}



summ_treat_stan_brms <- function(fit, treat_name = "XTreat", alpha = 0.05) {
  draws <- as.numeric(as.matrix(fit, pars = treat_name))

  est_mean <- mean(draws)
  est_sd   <- sd(draws)
  ci    <- quantile(draws, probs = c(alpha/2, 1 - alpha/2), names = FALSE)

  # two-sided "reject H0: beta=0" if 0 not in (1-alpha) CrI
  reject <- (ci[1] > 0) || (ci[2] < 0)

  # optional: posterior probability of benefit
  post_prob_gt0 <- mean(draws > 0)

  c(
    est_mean = est_mean,
    est_sd   = est_sd,
    ci_low   = ci[1],
    ci_high  = ci[2],
    reject   = as.numeric(reject),
    post_prob_gt0 = post_prob_gt0
  )
}


summ_treat_glm <- function(fit, treat_name = "XTreat", alpha = 0.05) {
  cf <- summary(fit)$coef[treat_name, ]

  est_mean <- cf["Estimate"]
  est_sd   <- cf["Std. Error"]

  tcrit <- qt(1 - alpha/2, df = fit$df.residual)
  ci    <- est_mean + c(-1, 1) * tcrit * est_sd

  # two-sided reject H0: beta = 0 if 0 not in CI
  reject <- (ci[1] > 0) || (ci[2] < 0)

  # frequentist analogue of P(beta > 0)
  prob_gt0 <- pt(est_mean / est_sd, df = fit$df.residual)

  c(
    est_mean = as.numeric(est_mean),
    est_sd   = as.numeric(est_sd),
    ci_low   = ci[1],
    ci_high  = ci[2],
    reject   = as.numeric(reject),
    post_prob_gt0 = as.numeric(prob_gt0)
  )
}
