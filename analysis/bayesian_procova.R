# Bayesian PROCOVA 


bayesian_procova_analysis <- function(
  df_historical,
  df_current,
  model_fun = NULL,
  oracle = FALSE,
  model_covar_names = c("X1", "X2", "X3", "X4", "X5"), 
  K0H = NULL, # default: 1/N_H
  K1H = 100, # large = weakly informative
  K2H = NULL, # default: from data
  k = 100, # flat component variance
  nu0 = 1, # flat component df
  sigma2_0 = NULL, # flat component scale
  a_omega = 1, # beta prior for the weight
  b_omega = 1,
  iter_warmup = 1000,
  iter_sampling = 3000,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }


  # ---------------------------------------------------------------------------
  # STEP 1: Obtain prognostic scores  m_i = E[Y_i(0) | x_i]
  # ---------------------------------------------------------------------------

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

      X_new_historical <- df_historical[, model_covar_names, drop = FALSE]
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

  # ---------------------------------------------------------------------------
  # STEP 2: Extract data vectors
  # ---------------------------------------------------------------------------
 
  y <- df_current$y
  w <- df_current$XTreat
  m <- df_current$progn_score
  y_H <- df_historical$y
  m_H <- df_historical$progn_score

  N <- length(y)
  N_H <- length(y_H)

  # ---------------------------------------------------------------------------
  # STEP 3: Centre outcomes and build the design matrix V 
  # ---------------------------------------------------------------------------
  # From Eq. (7): centre by m_bar so that beta_0 is interpretable as the
  # average bias of the prognostic score for control participants.

  m_bar <- mean(m)
  y_c <- y - m_bar
  V <- cbind(1, w, m - m_bar)

 # ---------------------------------------------------------------------------
  # STEP 4: Estimate informative prior parameters from historical control data (Section 3.2).
  # ---------------------------------------------------------------------------


  m_bar_H <- mean(m_H)
  fit_H <- lm(I(y_H - m_bar_H) ~ I(m_H - m_bar_H))
  beta_H <- c(coef(fit_H)[1], 0, coef(fit_H)[2])
  s2_H <- sum(residuals(fit_H)^2) / (N_H - 2)

  # ---------------------------------------------------------------------------
  # STEP 5: Set hyperparameters for K = diag(K0H, K1H, K2H)
  # ---------------------------------------------------------------------------

  if (is.null(K0H)) {
    K0H <- 1 / N_H
  }
  if (is.null(K2H)) {
    K2H <- 1 / sum((m_H - m_bar_H)^2)
  }
  K <- diag(c(K0H, K1H, K2H))

  # --- sigma2_0: prior point estimate of sigma^2 (paper Section 3.2) ---
  # Paper uses sigma2_0=1 because their sigma^2=1
  # Paper-faithful for general sigma^2: use s2_H
  if (is.null(sigma2_0)) {
    sigma2_0 <- s2_H
  }


  # ===========================================================================
  # PRECOMPUTE LOOP-INVARIANT QUANTITIES
  # ===========================================================================

  # Identity matrix (never changes)
  I_N <- diag(N)

  # K diagonal (used in posterior updates and omega draw)
  K_d <- diag(K)

  # ---------------------------------------------------------------------------
  # Informative component: M_I = I_N + V K V^T  (N x N)
  # ---------------------------------------------------------------------------
  # The marginal distribution of y^(c) under the informative prior is
  # (Section 3.2):
  #   y^(c) | omega=1 ~ t_{N_H - 2}( V * beta_H,  s2_H * (I_N + V K V^T) )
  #
  # M_I and its inverse appear in:
  #   - The marginal log-likelihood log_I (integrated over beta and sigma^2)
  #   - The posterior mean and covariance of beta (Section 3.4, Eq. for beta_*)
 
  M_I     <- I_N + V %*% K %*% t(V)   # N x N, fixed
  M_I_inv <- solve(M_I)                # O(N^3) -- done ONCE

  # ---------------------------------------------------------------------------
  # Flat component: M_F = I_N + k V V^T  (N x N)
  # ---------------------------------------------------------------------------
  # The marginal distribution of y^(c) under the flat prior is (Section 3.2):
  #   y^(c) | omega=0 ~ t_{nu0}( m_bar * 1_N,  sigma2_0 * (I_N + k V V^T) )
  #
  # M_F and its inverse appear analogously for the flat component.
 
  M_F     <- I_N + k * V %*% t(V)     # N x N, fixed
  M_F_inv <- solve(M_F)                # O(N^3) -- done ONCE

 
  # Log-determinants of M_I and M_F (needed in log_I and log_F below).
  log_det_M_I <- determinant(M_I, log = TRUE)$modulus
  log_det_M_F <- determinant(M_F, log = TRUE)$modulus


  # ---------------------------------------------------------------------------
  # Marginal log-likelihood under the informative component: log_I
  # ---------------------------------------------------------------------------
  # log_I = log integral L(beta, sigma^2 | data) * p_I(beta, sigma^2) d(beta) d(sigma^2)
  #
  # Closed-form (Appendix, C^{-1} formula, informative term):
  #
  #   log_I = log Gamma((N + N_H - 2)/2)  -  log Gamma((N_H - 2)/2)
  #           - N/2 * log((N_H - 2) * pi * s2_H)
  #           - 1/2 * log|M_I|
  #           - (N + N_H - 2)/2 * log(1 + q_I / ((N_H - 2) * s2_H))
  #
  # where:
  #   r_I = y^(c) - V * beta_H          (residual from prior mean)
  #   q_I = r_I^T * M_I^{-1} * r_I     (scaled quadratic form)
 
  r_I  <- y_c - V %*% beta_H
  q_I  <- sum(r_I * (M_I_inv %*% r_I))

  log_I <- lgamma((N + N_H - 2) / 2) -
    lgamma((N_H - 2) / 2) -
    N / 2 * log((N_H - 2) * pi * s2_H) -
    0.5 * log_det_M_I -
    (N + N_H - 2) / 2 * log(1 + q_I / ((N_H - 2) * s2_H))

  # ---------------------------------------------------------------------------
  # Marginal log-likelihood under the flat component: log_F
  # ---------------------------------------------------------------------------
  # Analogous closed-form for the flat component:
  #
  #   log_F = log Gamma((N + nu0)/2)  -  log Gamma(nu0/2)
  #           - N/2 * log(nu0 * pi * sigma2_0)
  #           - 1/2 * log|M_F|
  #           - (N + nu0)/2 * log(1 + q_F / (nu0 * sigma2_0))
  #
  # where:
  #   q_F = y^(c)^T * M_F^{-1} * y^(c)   (prior mean for flat component is 0)
 
  q_F  <- sum(y_c * (M_F_inv %*% y_c))
  log_F <- lgamma((N + nu0) / 2) -
    lgamma(nu0 / 2) -
    N / 2 * log(nu0 * pi * sigma2_0) -
    0.5 * log_det_M_F -
    (N + nu0) / 2 * log(1 + q_F / (nu0 * sigma2_0))
  
  
  # ---------------------------------------------------------------------------
  # Posterior mean and covariance building blocks — informative component
  # ---------------------------------------------------------------------------
  # From Section 3.4, the posterior mean of beta under the informative component
  # (conditional on sigma^2) is (Eq. for beta_*):
  #
  #   beta_* = beta_H + K V^T (I_N + V K V^T)^{-1} (y^(c) - V beta_H)
  #          = beta_H + KVtMinvI * (y^(c) - V beta_H)
  #
  # where KVtMinvI = K V^T M_I^{-1}  (3 x N precomputed matrix).
  #
  # The posterior covariance (independent of sigma^2 up to the sigma^2 factor) is:
  #   Sigma_* = sigma^2 * [K - K V^T M_I^{-1} V K]


  KVt      <- K %*% t(V)              # 3 x N
  KVtMinvI <- KVt %*% M_I_inv         # 3 x N
  Sigma_s_I <- K - KVtMinvI %*% V %*% K  # 3 x 3 posterior cov (informative)

  # ---------------------------------------------------------------------------
  # Posterior mean and covariance building blocks — flat component
  # ---------------------------------------------------------------------------
  # From Section 3.4, the posterior mean of beta under the flat component is:
  #
  #   b_* = k V^T (I_N + k V V^T)^{-1} y^(c)
  #       = kVtMinvF * y^(c)
  #
  # where kVtMinvF = k V^T M_F^{-1}  (3 x N precomputed matrix).
  #
  # The posterior covariance scale is:
  #   S_* = sigma^2 * [k I_{3x3} - k^2 V^T M_F^{-1} V]

  kVt       <- k * t(V)               # 3 x N
  kVtMinvF  <- kVt %*% M_F_inv        # 3 x N
  S_s_F     <- k * diag(3) - k * kVtMinvF %*% V  # 3 x 3 posterior cov (flat)
  # note: k^2 * t(V) %*% M_F_inv %*% V = k * (k * t(V) %*% M_F_inv) %*% V


  # ---------------------------------------------------------------------------
  # Precompute omega grid for the conditional omega draw (Step 4 of Gibbs)
  # ---------------------------------------------------------------------------
  # p(omega | beta, sigma^2, data) ∝ p(omega) * [omega * f_I + (1-omega) * f_F]
  om_grid     <- seq(0.001, 0.999, length.out = 500)
  log_om      <- log(om_grid)
  log_1m_om   <- log(1 - om_grid)
  log_prior_om <- dbeta(om_grid, a_omega, b_omega, log = TRUE)



  # ===========================================================================
  # STORAGE AND INITIALISATION
  # ===========================================================================
 
  n_total <- iter_warmup + iter_sampling
  beta_out <- matrix(NA, iter_sampling, 3)
  sigma2_out <- omega_out <- Z_out <- numeric(iter_sampling)

  omega <- 0.5

  # ===========================================================================
  # GIBBS SAMPLER  (Section 3.5)
  # ===========================================================================
  for (it in 1:n_total) {
    # -------------------------------------------------------------------------
    # GIBBS STEP 1: Compute the posterior mixture weight omega_*
    # -------------------------------------------------------------------------
    # Given the current omega, the posterior weight on the informative component
    # (Section 3.4) is:
    #
    #   omega_* = [omega * C_I^{-1}] / [omega * C_I^{-1} + (1-omega) * C_F^{-1}]
    #
    # where C_I^{-1} = exp(log_I) and C_F^{-1} = exp(log_F) are the marginal
    # likelihoods integrated over (beta, sigma^2) for each component.
    #
    # Numerically stable computation via log-sum-exp:

    max_log <- max(log(omega) + log_I, log(1 - omega) + log_F)
    omega_star <- exp(log(omega) + log_I - max_log) /
      (exp(log(omega) + log_I - max_log) +
        exp(log(1 - omega) + log_F - max_log))

    # -------------------------------------------------------------------------
    # GIBBS STEP 2: Draw component indicator Z ~ Bernoulli(omega_*)
    # -------------------------------------------------------------------------
    # Z = 1 -> draw from informative posterior component
    # Z = 0 -> draw from flat posterior component
 
    Z <- rbinom(1, 1, omega_star)

    # -------------------------------------------------------------------------
    # GIBBS STEP 3a (informative component, Z = 1):
    #   Draw sigma^2 then beta from their conjugate posteriors
    # -------------------------------------------------------------------------
    if (Z == 1) {

      # Posterior mean of beta under the informative component (Section 3.4):
      #   beta_* = beta_H + K V^T M_I^{-1} (y^(c) - V beta_H)
      beta_s <- beta_H + KVtMinvI %*% (y_c - V %*% beta_H)

      # Posterior scale for sigma^2 under the informative component (Section 3.4):
      #   s^2_* = 1/(N + N_H - 2) * [
      #             (y^(c) - V beta_*)^T (y^(c) - V beta_*)     <- RCT residual SS
      #             + (N_H - 2) * s2_H                           <- historical SS
      #             + (beta_*_0 - beta_hat_0H)^2 / K0H           <- prior penalty on beta_0
      #             + beta_*_1^2 / K1H                           <- prior penalty on beta_1
      #             + (beta_*_2 - beta_hat_2H)^2 / K2H           <- prior penalty on beta_2
      #           ]

      r <- y_c - V %*% beta_s   # (y^(c) - V beta_*)
      s2_s <- (sum(r^2) +
        (N_H - 2) * s2_H +
        (beta_s[1] - beta_H[1])^2 / K_d[1] +
        beta_s[2]^2 / K_d[2] +
        (beta_s[3] - beta_H[3])^2 / K_d[3]) /
        (N + N_H - 2)

      # Draw sigma^2 from its marginal posterior (scaled inverse chi-squared):
      sigma2 <- (N + N_H - 2) * s2_s / rchisq(1, N + N_H - 2)

      # Draw beta from its conditional posterior (multivariate normal):
      beta <- mvrnorm(1, beta_s, sigma2 * Sigma_s_I)

    } else {
      # -----------------------------------------------------------------------
      # GIBBS STEP 3b (flat component, Z = 0):
      #   Draw sigma^2 then beta from their conjugate posteriors
      # -----------------------------------------------------------------------
 
      # Posterior mean of beta under the flat component (Section 3.4):
      #   b_* = k V^T M_F^{-1} y^(c)   (kVtMinvF precomputed)
      
      b_s <- kVtMinvF %*% y_c

      # Posterior scale for sigma^2 under the flat component (Section 3.4):
      #   sigma^2_{0,*} = 1/(N + nu0) * [
      #                     (y^(c) - V b_*)^T (y^(c) - V b_*)    <- RCT residual SS
      #                     + nu0 * sigma2_0                       <- prior SS
      #                     + (b_0^2 + b_1^2 + b_2^2) / k         <- prior penalty on beta
      #                   ]
      r <- y_c - V %*% b_s
      s2_s <- (sum(r^2) + nu0 * sigma2_0 + sum(b_s^2) / k) / (N + nu0)

      # Draw sigma^2:
      #   sigma^2 | y, w, X ~ (N + nu0) * sigma^2_{0,*} / chi^2_{N + nu0}
      sigma2 <- (N + nu0) * s2_s / rchisq(1, N + nu0)

      # Draw beta | sigma^2:
      #   beta | sigma^2, y, w, X ~ N(b_*, sigma^2 * S_*)
      beta <- mvrnorm(1, b_s, sigma2 * S_s_F)
    }

    # -------------------------------------------------------------------------
    # GIBBS STEP 4: Draw omega from p(omega | beta, sigma^2, data)
    # -------------------------------------------------------------------------
    # The conditional posterior of omega (Section 3.4) is:
    #
    #   p(omega | beta, sigma^2, data) ∝ p(omega) * [omega * f_I(beta, sigma^2)
    #                                               + (1-omega) * f_F(beta, sigma^2)]
    #
    # where the component densities evaluated at the current (beta, sigma^2) are:
    #
    #   f_I(beta, sigma^2):
    #     log f_I ∝ -((N_H + 3)/2) * log(sigma^2)
    #               - 1/(2*sigma^2) * [(N_H-2)*s2_H
    #                                   + (beta_0 - beta_hat_0H)^2/K0H
    #                                   + beta_1^2/K1H
    #                                   + (beta_2 - beta_hat_2H)^2/K2H]
    #
    #   f_F(beta, sigma^2):
    #     log f_F ∝ -((nu0 + 3)/2 + 1) * log(sigma^2)
    #               - 1/(2*sigma^2) * [nu0*sigma2_0 + (beta_0^2+beta_1^2+beta_2^2)/k]
    #
    # The normalising constant is computed by numerical integration over the
    # grid om_grid, then omega is drawn via the probability integral transform
    # (inverse CDF method).
 
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

    log_pF <- -((nu0 + 5) / 2) *
      log(sigma2) -
      (nu0 * sigma2_0 + sum(beta^2) / k) / (2 * sigma2)

    # Vectorised evaluation over om_grid using log-sum-exp for numerical stability:
    #   log p(omega_j | beta, sigma^2) ∝ log p(omega_j)
    #       + log[ omega_j * exp(log_pI) + (1 - omega_j) * exp(log_pF) ]
 
    max_l    <- pmax(log_om + log_pI, log_1m_om + log_pF)
    log_lik  <- max_l + log(exp(log_om + log_pI - max_l) +
                              exp(log_1m_om + log_pF - max_l))
    log_post <- log_prior_om + log_lik
    post     <- exp(log_post - max(log_post))

    post <- exp(log_post - max(log_post))

    # Compute CDF and draw via probability integral transform (Step 5, Section 3.5):
    cdf <- cumsum(post) / sum(post)

    ## prob integral transform via the inverse of the cdf
    omega <- om_grid[which(cdf >= runif(1))[1]]

    # -------------------------------------------------------------------------
    # Store post-warmup draws
    # -------------------------------------------------------------------------
    if (it > iter_warmup) {
      idx <- it - iter_warmup
      beta_out[idx, ] <- beta
      sigma2_out[idx] <- sigma2
      omega_out[idx] <- omega
      Z_out[idx] <- Z
    }
  }

  colnames(beta_out) <- c("beta0", "beta1", "beta2")

  # Return posterior draws and setup metadata.
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



library(rjags)
library(coda)

bayesian_procova_jags <- function(
  df_historical,
  df_current,
  model_fun   = NULL,
  model_covar_names = c("X1", "X2", "X3", "X4", "X5"), 
  oracle      = FALSE,
  K0H         = NULL,
  K1H         = 100,
  K2H         = NULL,
  k           = 100,
  nu0         = 1,
  sigma2_0    = 1,
  a_omega     = 1,
  b_omega     = 1,
  iter_warmup = 1000,
  iter_sampling = 3000,
  seed        = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # ── Step 1: prognostic scores (identical to your implementation) ────────────
  if (oracle) {
    df_current$progn_score    <- df_current$S
    df_historical$progn_score <- df_historical$S
  } else {
    prog_fit <- model_fun(df_historical = df_historical, seed = seed)
    if (inherits(prog_fit, "SuperLearner")) {
      df_current$progn_score    <- as.numeric(predict(prog_fit, newdata = df_current[, model_covar_names, drop=FALSE])$pred)
      df_historical$progn_score <- as.numeric(predict(prog_fit, newdata = df_historical[, model_covar_names, drop=FALSE])$pred)
    } else {
      df_current$progn_score    <- as.numeric(predict(prog_fit, newdata = df_current))
      df_historical$progn_score <- as.numeric(predict(prog_fit, newdata = df_historical))
    }
  }

  # ── Step 2: extract vectors ─────────────────────────────────────────────────
  y   <- df_current$y
  w   <- df_current$XTreat
  m   <- df_current$progn_score
  y_H <- df_historical$y
  m_H <- df_historical$progn_score

  N   <- length(y)
  N_H <- length(y_H)

  # ── Step 3: centre and build design matrix ──────────────────────────────────
  m_bar   <- mean(m)
  y_c     <- y - m_bar
  m_c     <- m - m_bar
  m_bar_H <- mean(m_H)
  m_c_H   <- m_H - m_bar_H
  y_c_H   <- y_H - m_bar_H

  # ── Step 4: historical hyperparameters ──────────────────────────────────────
  fit_H   <- lm(y_c_H ~ m_c_H)
  beta0_H <- coef(fit_H)[1]   # beta-hat_0H
  beta2_H <- coef(fit_H)[2]   # beta-hat_2H
  s2_H    <- sum(residuals(fit_H)^2) / (N_H - 2)

  if (is.null(K0H)) K0H <- 1 / N_H
  if (is.null(K2H)) K2H <- 1 / sqrt(sum((m_c_H)^2)) # 1 / sum(m_c_H^2)

  # ── Step 5: package data for JAGS ───────────────────────────────────────────
  # JAGS will sample β = (β0, β1, β2) and σ² within each component,
  # selected by the latent Z.  We pass the informative prior parameters
  # and flat prior parameters as data.

  jags_data <- list(
    # Observed RCT data
    y_c   = y_c,
    w     = w,
    m_c   = m_c,
    N     = N,

    # Informative prior parameters
    beta_H   = c(beta0_H, 0, beta2_H),   # prior mean vector
    K0H      = K0H,
    K1H      = K1H,
    K2H      = K2H,
    N_H      = N_H,
    s2_H     = s2_H,

    # Flat prior parameters
    k        = k,
    nu0      = nu0,
    sigma2_0 = sigma2_0,

    # Mixture weight prior
    a_omega  = a_omega,
    b_omega  = b_omega
  )

  # ── Step 6: JAGS model string ────────────────────────────────────────────────
  # Design: Z ~ Bernoulli(omega) selects the component.
  # Component densities are encoded via the ones trick so JAGS can handle
  # the non-standard joint prior on (beta, sigma2).
  #
  # The likelihood is: y_c[i] ~ N(beta0 + beta1*w[i] + beta2*m_c[i],  sigma2)
  # The prior on (beta, sigma2) is the mixture selected by Z.
  #
  # We use the "ones trick":
  #   ones[i] ~ Bernoulli(p)   with   p = target_density / C
  # where C is a large enough constant so that p < 1 always holds.

  model_string <- "
  model {

    # ── Mixture weight ──────────────────────────────────────────────
    omega ~ dbeta(a_omega, b_omega)

    # ── Component indicator ─────────────────────────────────────────
    Z ~ dbern(omega)

    # ── Variance parameter ──────────────────────────────────────────
    # Informative component: sigma2 ~ scaled-Inv-chi2(N_H-2, s2_H)
    #   Equivalent: (N_H-2)*s2_H / sigma2_I  ~  chi2(N_H-2)
    #   In JAGS:  chi2_I_raw ~ dgamma((N_H-2)/2, 1/2) (shape, rate)
    #             then sigma2_I = (N_H-2)*s2_H / chi2_I_raw
    tau_I ~ dgamma( (N_H-2)/2,  (N_H-2)*s2_H/2 )
    sigma2_I <- 1 / tau_I         # = s2_H * (N_H-2) / chi2 draw

    # Flat component: sigma2 ~ scaled-Inv-chi2(nu0, sigma2_0)
    tau_F    ~ dgamma( nu0/2,  nu0 * sigma2_0 / 2 )
    sigma2_F <- 1 / tau_F

    # Active sigma2 determined by Z
    sigma2 <- Z * sigma2_I + (1 - Z) * sigma2_F

    # ── Regression coefficients ─────────────────────────────────────
    # Informative component: beta | sigma2_I ~ MVN(beta_H, sigma2_I * K)
    beta0_I ~ dnorm(beta_H[1],  1 / (sigma2_I * K0H))
    beta1_I ~ dnorm(beta_H[2],  1 / (sigma2_I * K1H))
    beta2_I ~ dnorm(beta_H[3],  1 / (sigma2_I * K2H))

    # Flat component: beta | sigma2_F ~ MVN(0, sigma2_F * k * I)
    beta0_F ~ dnorm(0,  1 / (sigma2_F * k))
    beta1_F ~ dnorm(0,  1 / (sigma2_F * k))
    beta2_F ~ dnorm(0,  1 / (sigma2_F * k))

    # Active coefficients determined by Z
    beta0 <- Z * beta0_I + (1 - Z) * beta0_F
    beta1 <- Z * beta1_I + (1 - Z) * beta1_F
    beta2 <- Z * beta2_I + (1 - Z) * beta2_F

    # ── Likelihood ──────────────────────────────────────────────────
    for (i in 1:N) {
      mu[i] <- beta0 + beta1 * w[i] + beta2 * m_c[i]
      y_c[i] ~ dnorm(mu[i],  1 / sigma2)
    }
  }
  "

  # ── Step 7: compile and run ──────────────────────────────────────────────────
  jags_model <- jags.model(
    textConnection(model_string),
    data    = jags_data,
    n.chains = 1,
    n.adapt  = iter_warmup,
    quiet    = FALSE
  )

  # burn-in
  update(jags_model, n.iter = iter_warmup)

  # sampling
  jags_samples <- coda.samples(
    jags_model,
    variable.names = c("beta0", "beta1", "beta2", "sigma2", "omega", "Z"),
    n.iter = iter_sampling * 10,  # increase iterations to compensate
    thin   = 10
  )

  # ── Step 8: extract and return in same format as your implementation ─────────
  draws_mat <- as.matrix(jags_samples)

  list(
    draws = list(
      beta   = draws_mat[, c("beta0", "beta1", "beta2")],
      sigma2 = draws_mat[, "sigma2"],
      omega  = draws_mat[, "omega"],
      Z      = draws_mat[, "Z"]
    ),
    setup = list(
      N      = N,
      N_H    = N_H,
      beta_H = c(beta0_H, 0, beta2_H),
      s2_H   = s2_H
    ),
    jags_samples = jags_samples   # keep mcmc object for diagnostics
  )
}
