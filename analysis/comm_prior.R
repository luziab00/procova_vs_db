# commensurate prior dynamic borrowing 


comm_prior_db_analysis_psborrow2 <- function(
  df_historical, 
  df_current, 
  seed = NULL,
  tau_rate = 0.01,
  iter_warmup = 1000, 
  iter_sampling = 10000
){

  # combine the dfs
  df_all<-rbind(cbind(df_current, ext=0), cbind(df_historical, ext=1))

  # create matrix object
  df_all_matrix <- create_data_matrix(
  df_all,
  outcome = c("y"),
  trt_flag_col = "XTreat",
  ext_flag_col = "ext"
)
  
  # create outcome 
  outcome <- outcome_cont_normal(
    continuous_var = "y",
    baseline_prior = prior_normal(0, 1000),
    std_dev_prior = prior_half_cauchy(0, 50)
  )

treatment <- treatment_details(
    trt_flag_col = "XTreat",
    trt_prior    = prior_normal(0, 1000)
  )



  # FULL BORROWING

  borrowing_full_obj <- borrowing_full(ext_flag_col = "ext")

  anls_full <- create_analysis_obj(
    data_matrix = df_all_matrix,
    outcome     = outcome,
    borrowing   = borrowing_full_obj,
    treatment   = treatment,
    quiet       = TRUE,
  )

  fit_full <- suppressMessages(
  suppressWarnings(
    quietly(
      mcmc_sample(
    anls_full,
    iter_warmup = 1000, 
    iter_sampling = 10000, 
    chains = 4,
    seed          = if (is.null(seed)) NULL else (seed + 2L),
    verbose       = FALSE
   )
    )
  )
)

  # COMMENSURATE
  borrowing_comm_obj <- borrowing_hierarchical_commensurate(
    ext_flag_col = "ext",
    tau_prior    = prior_gamma(alpha = 1, beta = tau_rate)
  )

  anls_comm <- create_analysis_obj(
    data_matrix = df_all_matrix,
    outcome     = outcome,
    borrowing   = borrowing_comm_obj,
    treatment   = treatment,
    quiet       = TRUE
  )

  fit_comm <- suppressMessages(
  suppressWarnings(
    quietly(
      mcmc_sample(
    anls_comm,
    iter_warmup = iter_warmup, 
    iter_sampling = iter_sampling, 
    chains = 4,
    seed          = if (is.null(seed)) NULL else (seed + 3L),
    verbose       = FALSE
   )
    )
  )
)

  list(
    fit_full = fit_full,
    fit_comm = fit_comm
  )
}



summ_treat_comm_psborrow2 <- function(
  fit,
  treat_name = "XTreat",
  alpha = 0.05
) {
  # extract posterior draws as a data.frame
  draws_df <- posterior::as_draws_df(fit$draws())

  draws <- as.numeric(draws_df[["beta_trt"]])

  est_mean <- mean(draws)
  est_sd   <- sd(draws)

  ci    <-  quantile(draws, probs = c(alpha/2, 1 - alpha/2), names = FALSE)

  # two-sided "reject H0: beta=0" if 0 not in (1-alpha) CrI
  reject <- (ci[1] > 0) || (ci[2] < 0)

  # optional: posterior probability of benefit
  post_prob_gt0 <- mean(draws > 0)

  c(
    est_mean        = est_mean,
    est_sd          = est_sd,
    ci_low          = ci[1],
    ci_high         = ci[2],
    reject          = as.numeric(reject),
    post_prob_gt0   = post_prob_gt0
  )
}




## commensurate prior dynamic borrowing using RBest

comm_prior_db_fit_RBest <- function(
  df_historical,
  df_current,
  seed = NULL,
  tau_rate = 0.01,
  tau_grid_m = 200,
  huge_sd = 1e6
) {
  if (!is.null(seed)) set.seed(seed)

  # split arms
  y_h_c <- df_historical$y[df_historical$XTreat == 0]
  y_c_c <- df_current$y[df_current$XTreat == 0]
  y_c_t <- df_current$y[df_current$XTreat == 1]

  mu_h <- mean(y_h_c)
  se_h <- sd(y_h_c) / sqrt(length(y_h_c))  # SE of historical control mean

  # plug-in common sigma if not provided
  sigma <- sd(c(y_h_c, y_c_c, y_c_t))


  # =========================
  # 1) NO BORROWING (vague prior on current control mean)
  # =========================
  prior_none_ctrl <- mixnorm(c(1, 0, huge_sd), sigma = sigma)
  post_none_ctrl  <- postmix(prior = prior_none_ctrl, data = y_c_c, sigma = sigma)

  # =========================
  # 2) FULL BORROWING (prior N(mu_h, se_h^2) for current control mean)
  # =========================
  prior_full_ctrl <- mixnorm(c(1, mu_h, se_h), sigma = sigma)
  post_full_ctrl  <- postmix(prior = prior_full_ctrl, data = y_c_c, sigma = sigma)

  # =========================
  # 3) COMMENSURATE borrowing via mixture over tau
  # tau ~ Gamma(shape=1, rate=tau_rate)
  # mu_c | tau ~ N(mu_h, se_h^2 + 1/tau)
  # =========================
  shape <- 1
  probs <- seq(0.001, 0.999, length.out = tau_grid_m)
  tau   <- qgamma(probs, shape = shape, rate = tau_rate)

  sd_grid <- sqrt(se_h^2 + 1 / tau)

  dens <- dgamma(tau, shape = shape, rate = tau_rate)
  dtau <- c(diff(tau), tail(diff(tau), 1))
  w    <- dens * dtau
  w    <- w / sum(w)

  # pack (w_i, mu_h, sd_i) triples
  comp_comm <- as.vector(rbind(w, rep(mu_h, tau_grid_m), sd_grid))
  prior_comm_ctrl <- mixnorm(comp_comm, sigma = sigma)
  post_comm_ctrl  <- postmix(prior = prior_comm_ctrl, data = y_c_c, sigma = sigma)

  # Treatment arm posterior with vague prior
  prior_trt <- mixnorm(c(1, 0, huge_sd), sigma = sigma)
  post_trt  <- postmix(prior = prior_trt, data = y_c_t, sigma = sigma)

  # Treatment effect mixes: trt - ctrl
  te_none <- pmixdiff(post_trt, post_none_ctrl)
  te_full <- pmixdiff(post_trt, post_full_ctrl)
  te_comm <- pmixdiff(post_trt, post_comm_ctrl)

  list(
    call = match.call(),
    sigma = sigma,
    inputs = list(mu_h = mu_h, se_h = se_h, tau_rate = tau_rate, tau_grid_m = tau_grid_m),
    post = list(
      ctrl_none = post_none_ctrl,
      ctrl_full = post_full_ctrl,
      ctrl_comm = post_comm_ctrl,
      trt       = post_trt
    ),
    te = list(
      none = te_none,
      full = te_full,
      comm = te_comm
    )
  )
}


comm_prior_db_summ_RBest <- function(
  fit,
  alpha = 0.05,
  n_draws = 200000,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)


  summ_one <- function(te_mix) {
    draws <- rmix(n_draws, te_mix)

    est_mean <- mean(draws)
    est_sd   <- sd(draws)
    ci       <- quantile(draws, probs = c(alpha/2, 1 - alpha/2), names = FALSE)

    reject <- (ci[1] > 0) || (ci[2] < 0)
    post_prob_gt0 <- mean(draws > 0)

    list(
      est_mean      = est_mean,
      est_sd        = est_sd,
      ci_low        = ci[1],
      ci_high       = ci[2],
      reject        = as.numeric(reject),
      post_prob_gt0 = post_prob_gt0
    )
  }

  list(
    sigma = fit$sigma,
    none  = summ_one(fit$te$none),
    full  = summ_one(fit$te$full),
    comm  = summ_one(fit$te$comm)
  )
}


