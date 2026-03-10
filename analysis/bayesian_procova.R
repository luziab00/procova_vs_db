# Bayesian PROCOVA 


bayesian_procova_analysis <- function(
  df_historical,
  df_current,
  model_fun = NULL,
  oracle = FALSE,
  K0H = NULL, # default: 1/N_H
  K1H = 100, # large = weakly informative
  K2H = NULL, # default: from data
  k = 100, # flat component variance
  nu0 = 1, # flat component df
  sigma2_0 = 1, # flat component scale
  a_omega = 1, # beta prior for the weight
  b_omega = 1,
  iter_warmup = 1000,
  iter_sampling = 3000,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (oracle) {
    df_current$progn_score <- df_current$S
    prog_fit <- lm(y ~ S, data = df_historical)

    df_historical$progn_score <- df_historical$S
  } else {
    prog_fit <- model_fun(
      df_historical = df_historical,
      seed = seed
    )

    if (inherits(prog_fit, "SuperLearner")) {
      X_new_current <- df_current[, model_covar_names, drop = FALSE]
      df_current$progn_score <- as.numeric(
        predict(prog_fit, newdata = X_new_current)$pred
      )

      X_new_historical <- df_current[, model_covar_names, drop = FALSE]
      df_historical$progn_score <- as.numeric(
        predict(prog_fit, newdata = X_new_historical)$pred
      )
    } else {
      df_current$progn_score <- as.numeric(predict(
        prog_fit,
        newdata = df_current
      ))
      df_historical$progn_score <- as.numeric(predict(
        prog_fit,
        newdata = df_historical
      ))
    }
  }

  # --- Extract data ---
  y <- df_current$y
  w <- df_current$XTreat
  m <- df_current$progn_score
  y_H <- df_historical$y
  m_H <- df_historical$progn_score

  N <- length(y)
  N_H <- length(y_H)

  # --- Center and create design matrix ---
  m_bar <- mean(m)
  y_c <- y - m_bar
  V <- cbind(1, w, m - m_bar)

  # --- Historical estimates (OLS) ---
  m_bar_H <- mean(m_H)
  fit_H <- lm(I(y_H - m_bar_H) ~ I(m_H - m_bar_H))
  beta_H <- c(coef(fit_H)[1], 0, coef(fit_H)[2])
  s2_H <- sum(residuals(fit_H)^2) / (N_H - 2)

  # --- Set hyperparameters ---
  if (is.null(K0H)) {
    K0H <- 1 / N_H
  }
  if (is.null(K2H)) {
    K2H <- 1 / sum((m_H - m_bar_H)^2)
  }
  K <- diag(c(K0H, K1H, K2H))



  # --- Precompute everything that is loop-invariant ---

  # Identity matrix (never changes)
  I_N <- diag(N)

  # K diagonal (used in posterior updates and omega draw)
  K_d <- diag(K)

  # Informative component: M_I and its inverse
  M_I     <- I_N + V %*% K %*% t(V)   # N x N, fixed
  M_I_inv <- solve(M_I)                # O(N^3) -- done ONCE

  # Flat component: M_F and its inverse
  M_F     <- I_N + k * V %*% t(V)     # N x N, fixed
  M_F_inv <- solve(M_F)                # O(N^3) -- done ONCE

  # Log-determinants (fixed)
  log_det_M_I <- determinant(M_I, log = TRUE)$modulus
  log_det_M_F <- determinant(M_F, log = TRUE)$modulus


  # Informative integral
  r_I  <- y_c - V %*% beta_H
  q_I  <- sum(r_I * (M_I_inv %*% r_I))

  log_I <- lgamma((N + N_H - 2) / 2) -
    lgamma((N_H - 2) / 2) -
    N / 2 * log((N_H - 2) * pi * s2_H) -
    0.5 * log_det_M_I -
    (N + N_H - 2) / 2 * log(1 + q_I / ((N_H - 2) * s2_H))

  # Flat integral
  q_F  <- sum(y_c * (M_F_inv %*% y_c))
  log_F <- lgamma((N + nu0) / 2) -
    lgamma(nu0 / 2) -
    N / 2 * log(nu0 * pi * sigma2_0) -
    0.5 * log_det_M_F -
    (N + nu0) / 2 * log(1 + q_F / (nu0 * sigma2_0))
  
  
  # Informative posterior mean/cov building blocks (fixed)
  KVt      <- K %*% t(V)              # 3 x N
  KVtMinvI <- KVt %*% M_I_inv         # 3 x N
  Sigma_s_I <- K - KVtMinvI %*% V %*% K  # 3 x 3 posterior cov (informative)

  # Flat posterior cov building block (fixed)
  kVt       <- k * t(V)               # 3 x N
  kVtMinvF  <- kVt %*% M_F_inv        # 3 x N
  S_s_F     <- k * diag(3) - k * kVtMinvF %*% V  # 3 x 3 posterior cov (flat)
  # note: k^2 * t(V) %*% M_F_inv %*% V = k * (k * t(V) %*% M_F_inv) %*% V

  # OPT 2: Precompute omega grid and its logs (vectorized omega draw)
  om_grid     <- seq(0.001, 0.999, length.out = 500)
  log_om      <- log(om_grid)
  log_1m_om   <- log(1 - om_grid)
  log_prior_om <- dbeta(om_grid, a_omega, b_omega, log = TRUE)



  # --- Storage ---
  n_total <- iter_warmup + iter_sampling
  beta_out <- matrix(NA, iter_sampling, 3)
  sigma2_out <- omega_out <- Z_out <- numeric(iter_sampling)

  # --- Initialize ---
  omega <- 0.5

  # --- Gibbs sampler ---
  for (it in 1:n_total) {
    # 1. Compute omega_* (weight for informative component)

    # omega_* = omega * exp(log_I) / [omega*exp(log_I) + (1-omega)*exp(log_F)]
    max_log <- max(log(omega) + log_I, log(1 - omega) + log_F)
    omega_star <- exp(log(omega) + log_I - max_log) /
      (exp(log(omega) + log_I - max_log) +
        exp(log(1 - omega) + log_F - max_log))

    # 2. Draw component indicator
    Z <- rbinom(1, 1, omega_star)

    # 3. Draw (beta, sigma2) from selected component
    if (Z == 1) {
      # Informative
      beta_s <- beta_H + KVtMinvI %*% (y_c - V %*% beta_H)

      r <- y_c - V %*% beta_s
      s2_s <- (sum(r^2) +
        (N_H - 2) * s2_H +
        (beta_s[1] - beta_H[1])^2 / K_d[1] +
        beta_s[2]^2 / K_d[2] +
        (beta_s[3] - beta_H[3])^2 / K_d[3]) /
        (N + N_H - 2)

      sigma2 <- (N + N_H - 2) * s2_s / rchisq(1, N + N_H - 2)
      beta <- mvrnorm(1, beta_s, sigma2 * Sigma_s_I)
    } else {
      # Flat
      b_s <- kVtMinvF %*% y_c

      r <- y_c - V %*% b_s
      s2_s <- (sum(r^2) + nu0 * sigma2_0 + sum(b_s^2) / k) / (N + nu0)

      sigma2 <- (N + nu0) * s2_s / rchisq(1, N + nu0)
      beta <- mvrnorm(1, b_s, sigma2 * S_s_F)
    }

    # 4. Draw omega | beta, sigma2 (numerical CDF inversion)

    ## define unnormalized density
    ## calculate the constant terms
    log_pI <- -((N_H + 3) / 2) *
      log(sigma2) -
      ((N_H - 2) *
        s2_H +
        (beta[1] - beta_H[1])^2 / K_d[1] +
        beta[2]^2 / K_d[2] +
        (beta[3] - beta_H[3])^2 / K_d[3]) /
        (2 * sigma2)

    log_pF <- -((nu0 + 4) / 2) *
      log(sigma2) -
      (nu0 * sigma2_0 + sum(beta^2) / k) / (2 * sigma2)

    # Vectorized over om_grid
    max_l    <- pmax(log_om + log_pI, log_1m_om + log_pF)
    log_lik  <- max_l + log(exp(log_om + log_pI - max_l) +
                              exp(log_1m_om + log_pF - max_l))
    log_post <- log_prior_om + log_lik
    post     <- exp(log_post - max(log_post))

    post <- exp(log_post - max(log_post))

    ## caluclate cdf (normalized)
    cdf <- cumsum(post) / sum(post)

    ## prob integral transform via the inverse of the cdf
    omega <- om_grid[which(cdf >= runif(1))[1]]

    # 5. Store post-warmup
    if (it > iter_warmup) {
      idx <- it - iter_warmup
      beta_out[idx, ] <- beta
      sigma2_out[idx] <- sigma2
      omega_out[idx] <- omega
      Z_out[idx] <- Z
    }
  }

  colnames(beta_out) <- c("beta0", "beta1", "beta2")

  list(
    draws = list(
      beta = beta_out,
      sigma2 = sigma2_out,
      omega = omega_out,
      Z = Z_out
    ),
    setup = list(N = N, N_H = N_H, beta_H = beta_H, s2_H = s2_H)
  )
}


summ_treat_bayesian_procova <- function(fit, alpha = 0.05) {
  draws <- fit$draws$beta[, "beta1"]
  est_mean <- mean(draws)
  est_sd <- sd(draws)
  ci <- quantile(draws, c(alpha / 2, 1 - alpha / 2))

  c(
    est_mean = est_mean,
    est_sd = est_sd,
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[2]),
    reject = as.numeric(ci[1] > 0 | ci[2] < 0),
    post_prob_gt0 = mean(draws > 0)
  )
}
