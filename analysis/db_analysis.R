# dynamic borrowing analyses

#install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
library(cmdstanr)

# Install the external CmdStan program
check_cmdstan_toolchain()
#install_cmdstan(cores = 2)


## commensurate prior dynamic borrowing using psborrow2


quietly <- function(expr) {
  val <- NULL
  capture.output(
    val <- force(expr),
    file = NULL
  )
  val
}



comm_prior_db_analysis_psborrow2 <- function(
  df_historical, 
  df_current, 
  seed = NULL,
  tau_rate = 0.01
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

  # =========================
  # 1) NO BORROWING
  # =========================
  borrowing_none_obj <- borrowing_none(ext_flag_col = "ext")

  anls_none <- create_analysis_obj(
    data_matrix = df_all_matrix,
    outcome     = outcome,
    borrowing   = borrowing_none_obj,
    treatment   = treatment,
    quiet = TRUE
  )

  fit_none <- suppressMessages(
  suppressWarnings(
    quietly(
      mcmc_sample(
        anls_none,
        iter_warmup   = 2000,
        iter_sampling = 50000,
        chains        = 4,
        seed          = if (is.null(seed)) NULL else (seed + 1L),
        verbose       = FALSE
   )
    )
  )
)

  # =========================
  # 2) FULL BORROWING
  # =========================
  borrowing_full_obj <- borrowing_full(ext_flag_col = "ext")

  anls_full <- create_analysis_obj(
    data_matrix = df_all_matrix,
    outcome     = outcome,
    borrowing   = borrowing_full_obj,
    treatment   = treatment,
    quiet       = TRUE
  )

  fit_full <- suppressMessages(
  suppressWarnings(
    quietly(
      mcmc_sample(
    anls_full,
    iter_warmup = 2000, 
    iter_sampling = 50000, 
    chains = 4,
    seed          = if (is.null(seed)) NULL else (seed + 2L),
    verbose       = FALSE
   )
    )
  )
)

  # =========================
  # 3) HIERARCHICAL COMMENSURATE
  # =========================
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
    iter_warmup = 2000, 
    iter_sampling = 50000, 
    chains = 4,
    seed          = if (is.null(seed)) NULL else (seed + 3L),
    verbose       = FALSE
   )
    )
  )
)

  list(
    fit_none = fit_none,
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

  list(
    fit_penalized = fit_penalized
  )
}


summ_treat_adaptive_lasso <- function(fit, treat_name = "XTreat", alpha = 0.05) {
  cf <- coef(fit)[treat_name,]

  est_mean <- cf
  #est_sd   <- cf["Std. Error"]

  #tcrit <- qt(1 - alpha/2, df = fit$df.residual)
  #ci    <- est_mean + c(-1, 1) * tcrit * est_sd

  # two-sided reject H0: beta = 0 if 0 not in CI
  #reject <- (ci[1] > 0) || (ci[2] < 0)

  # frequentist analogue of P(beta > 0)
  #prob_gt0 <- pt(est_mean / est_sd, df = fit$df.residual)

  c(
    est_mean = est_mean#,
    #est_sd   = est_sd,
    #ci_low   = ci[1],
    #ci_high  = ci[2],
    #reject   = as.numeric(reject),
    #prob_gt0 = as.numeric(prob_gt0)
  )
}





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

