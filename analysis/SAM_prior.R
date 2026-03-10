

# SAMprior dynamic borrowing
sam_analysis_rbest <- function(
  df_historical,
  df_current,
  y_col = "y",
  trt_col = "XTreat",
  delta_mult = 0.5
){
  y_hist <- df_historical[[y_col]]
  y_cc   <- df_current[df_current[[trt_col]] == 0, y_col, drop = TRUE]
  y_trt  <- df_current[df_current[[trt_col]] == 1, y_col, drop = TRUE]

  # historical prior component
  mu_h <- mean(y_hist)
  sd_h <- sd(y_hist)
  se_h <- sd_h / sqrt(length(y_hist))

  # current control summaries (for nf prior + sigma)
  mu_cc <- mean(y_cc) 
  sd_cc <- sd(y_cc)
  se_cc <- sd_cc / sqrt(length(y_cc))


  # current-only scale (no borrowing)
  sigma_unit <- sd(c(y_cc, y_trt))
  
  
  if.prior <- mixnorm(c(1, mu_h, se_h), sigma = sigma_unit)
  nf.prior <- mixnorm(c(1, 0, sigma_unit), sigma = sigma_unit)

  # --- 1) SAM dynamic borrowing (between if.prior and nf.prior) ---
  wSAM <- SAM_weight(if.prior = if.prior,
                     delta    = delta_mult * sd_cc,
                     data     = y_cc)

  SAM.prior <- SAM_prior(if.prior = if.prior, nf.prior = nf.prior,
                         weight = wSAM, sigma = sd_cc)

  post_cc_SAM  <- postmix(prior = SAM.prior, data = y_cc, sigma = sd_cc)

  # --- 2) Full borrowing (use historical prior directly) ---
  post_cc_full <- postmix(prior = if.prior, data = y_cc, sigma = sd_cc)

  # --- 3) No borrowing (use vague prior only) ---
  post_cc_none <- postmix(prior = nf.prior, data = y_cc, sigma = sd_cc)



  # treatment posterior (mirror your code: use historical prior)
  trt.prior <- mixnorm(c(1, mu_h, sd_h))
  post_trt  <- postmix(prior = trt.prior, data = y_trt, sigma = sd(y_trt))

  list(
    post_cc_SAM  = post_cc_SAM,
    post_cc_full = post_cc_full,
    post_cc_none = post_cc_none,
    post_trt     = post_trt,
    wSAM         = wSAM
  )
}



summ_treat_sam_rbest <- function(
  fit,
  n_draws = 20000,
  model = c("SAM", "full", "none"),
  seed = NULL,
  alpha = 0.05
){
  if (!is.null(seed)) set.seed(seed)
  
  # pick the right control posterior from the fit object
  post_cc <- switch(
    model,
    "SAM"  = fit$post_cc_SAM,
    "full" = fit$post_cc_full,
    "none" = fit$post_cc_none
  )

  delta <- rmix(n = n_draws, mix = fit$post_trt) - rmix(n = n_draws, mix = post_cc)

  est_mean <- mean(delta)
  est_sd   <- sd(delta)
  ci       <- unname(quantile(delta, c(alpha/2, 1 - alpha/2)))
  reject   <- (ci[1] > 0) || (ci[2] < 0)
  post_prob_gt0 <- mean(delta > 0)

  c(
    #sam_weight       = fit$wSAM,
    est_mean          = est_mean,
    est_sd            = est_sd,
    ci_low           = ci[1],
    ci_high          = ci[2],
    reject    = as.numeric(reject),
    post_prob_gt0 = post_prob_gt0
  )
}

