# adaptive lasso dynamic borrowing
adaptive_lasso_db_analysis<- function(
  df_historical, 
  df_current, 
  lambda = 0.001,
  gamma = 0.5,
  model_covar_names = NULL, 
  seed = NULL
){


  # combine the dfs
  df_all<-rbind(cbind(df_current, ext=0), cbind(df_historical, ext=1))

  # build formula
  f_covars <- as.formula(
    paste0("y", " ~ ", paste(c("XTreat", "ext", model_covar_names), collapse = " + "))
  )

  # fit inital model to get weights
  fit_inital <- lm(f_covars, data = df_all, x = TRUE) 
  predictor_matrix <- fit_inital$x

  # caluclate weights
  inital_coeff_ext <- as.numeric(coef(fit_inital)["ext"])
  weight_ext <- (abs(inital_coeff_ext))^-gamma

  # set all other weights to zero (since we dont want to penalize them)
  adaptive_weights <- rep(0, ncol(predictor_matrix))  
  adaptive_weights[3] <- weight_ext

  # since glmnet does automatically rescale the weights of the penalization to summarize to the number of penealized parameters
  # we standardize the weights here and adjust the general penalization strength (lambda) accordingly
  scaling_factor <- sum(adaptive_weights) / length(adaptive_weights)
  adjusted_adaptive_weights <- adaptive_weights / scaling_factor

  lambda_rescaled <- lambda * scaling_factor


  # fit second model with adaptive lasso 
  fit_penalized <- glmnet(
    x = predictor_matrix, y = df_all$y,
    family = "gaussian",
    alpha = 1,                    # lasso
    lambda = lambda_rescaled,             # single value allowed
    penalty.factor = adjusted_adaptive_weights
  )



  # extract ext_hat to decide borrowing
  b_pen <- as.matrix(coef.glmnet(fit_penalized, s = lambda_rescaled))
  ext_hat <- as.numeric(b_pen["ext", 1])

  # if coefficient was fully shrinked to 0
  full_borrowing <- abs(ext_hat) <  1e-8 

  # refit unpenalized model for SE/CI (drop ext if it was shrunk to 0)
  if (full_borrowing) {
    f_final <- as.formula(
      paste0("y", " ~ ", paste(c("XTreat", model_covar_names), collapse = " + "))
    )
  } else {
    f_final <- f_covars
  }
  fit_final <- lm(f_final, data = df_all)



  list(
    fit_penalized = fit_penalized,
    fit_final = fit_final,
    ext_hat = ext_hat,
    lambda_rescaled = lambda_rescaled
  )
}


summ_treat_adaptive_lasso <- function(fit, treat_name = "XTreat", alpha = 0.05) {
  est_mean <- coef(fit$fit_penalized)[treat_name,]


  #SE and DF from the Refitted Model (Oracle property)
  refit_sum <- summary(fit$fit_final)
  est_sd    <- refit_sum$coefficients[treat_name, "Std. Error"]
  df_res    <- fit$fit_final$df.residual


  # Critical value and Confidence Interval
  tcrit <- qt(1 - alpha/2, df = df_res)
  ci    <- est_mean + c(-1, 1) * tcrit * est_sd

  # two-sided reject H0: beta = 0 if 0 not in CI
  reject <- (ci[1] > 0) || (ci[2] < 0)

  # frequentist analogue of P(beta > 0)
  prob_gt0 <- pt(est_mean / est_sd, df = df_res)

  c(
    est_mean = est_mean,
    est_sd   = est_sd,
    ci_low   = ci[1],
    ci_high  = ci[2],
    reject   = as.numeric(reject),
    prob_gt0 = as.numeric(prob_gt0)
  )
}

