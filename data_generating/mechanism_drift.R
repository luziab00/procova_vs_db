
source("utils/helper_functions.R")
source("utils/__init__.R")



### configuration 
# start with a base setting - all datasets have the same covariate distribution


covars <- list(
  hist1  = list(mu = c(X1= 0, X2=0, X3=0, X4=0), sd = c(X1=1, X2=1, X3=1, X4=1), pX5 = 0.5),
  hist2  = list(mu = c(X1= 0, X2=0, X3= 0, X4=0), sd = c(X1=1, X2=1, X3=1, X4=1), pX5 = 0.5),
  current = list(mu = c(X1= 0, X2=0, X3= 0, X4=0), sd = c(X1=1, X2=1, X3=1, X4=1), pX5 = 0.5)
)

sample_size <- list(
  hist1 = 200,
  hist2 = 200,
  current = 400)


## params for current trial
treatment_effect <- 2  
randomization_ratio <- 0.5
residual_error_sd <- 0.05


coeff <- list(
  hist1  = list(beta1=1, beta2=-1, beta3=0, beta4=2, gamma1=0.5, gamma2=0.5, gamma3=0.5, gamma4=0.5),
  hist2  = list(beta1=1, beta2=-1, beta3=0, beta4=2, gamma1=0.5, gamma2=0.5, gamma3=0.5, gamma4=0.5),
  current = list(beta1=1, beta2=-1, beta3=0, beta4=2, gamma1=0.5, gamma2=0.5, gamma3=0.5, gamma4=0.5)
)


gen_covs <- function(n, covars) {
  X1 <- rnorm(n, covars$mu["X1"], covars$sd["X1"])
  X2 <- rnorm(n, covars$mu["X2"], covars$sd["X2"])
  X3 <- rnorm(n, covars$mu["X3"], covars$sd["X3"])
  X4 <- rnorm(n, covars$mu["X4"], covars$sd["X4"])
  X5  <- rbinom(n, 1, covars$pX5)
  data.frame(X1, X2, X3, X4, X5)
}

hist1 <- gen_covs(sample_size$hist1, covars$hist1) 
hist2 <- gen_covs(sample_size$hist2, covars$hist2) 
current <- gen_covs(sample_size$current, covars$current) 


hist1$XTreat <- rbinom(sample_size$hist1, 1, 0)
hist2$XTreat <- rbinom(sample_size$hist2, 1, 0)
current$XTreat <- rbinom(sample_size$current, 1, randomization_ratio)


outcome_from_covs_scaled_residuals <- function(
  df,
  coeffs,
  residual_error_sd,
  treatment_effect = 0,
  target_cor = NULL  # if not NULL, enforce this correlation between S and Y (control arm)
) {
  # df must contain: X1, X2, X3, X4, X5, XTreat (0/1)
  n <- nrow(df)

  # ----- Prognostic score S(X) from the new DGM -----
  S <- with(df, {
      coeffs$beta1  * X1 +
      coeffs$beta2  * X2 +
      coeffs$beta3  * X3 +
      coeffs$beta4  * X4 +
      coeffs$gamma1 * (X2^2) +
      coeffs$gamma2 * sin(X3) +
      coeffs$gamma3 * as.numeric(X4 > 1) +
      coeffs$gamma4 * X1 * X5
  })

  # ----- Option 1: standard DGM (no forced correlation) -----
  if (is.null(target_cor)) {
    residual_error <- rnorm(n, mean = 0, sd = residual_error_sd)
    Y <- S + treatment_effect * df$XTreat + residual_error
    return(Y)
  }

  # ----- Option 2: enforce sample correlation between S and control outcome -----
  # We construct Y0 such that cor(S, Y0) == target_cor in the sample.

  # standardize prognostic score
  S_std <- as.numeric(scale(S))

  # draw noise and make it orthogonal to S
  U      <- rnorm(n)
  R      <- residuals(lm(U ~ S_std))         # R ⟂ S_std in sample
  R_std  <- as.numeric(scale(R))

  # core control outcome with exact correlation target_cor
  Y0_core <- target_cor * S_std + sqrt(1 - target_cor^2) * R_std

  # rescale / shift: residual_error_sd controls overall scale,
  # mean(S) sets the approximate level (does not affect correlation)
  Y0 <- mean(S) + residual_error_sd * Y0_core

  # add treatment effect
  Y <- Y0 + treatment_effect * df$XTreat

  return(Y)
}



outcome_from_covs_adjusted_error_sd <- function(
  df,
  coeffs,
  treatment_effect = 0,
  residual_error_sd = NULL,  # set directly, OR
  target_cor = NULL          # set desired population correlation |r| in (0, 1)
) {
  # df must contain: X1, X2, X3, X4, X5, XTreat (0/1)
  n <- nrow(df)

  # ----- Prognostic score S(X) from the new DGM -----
  S <- with(df, {
      coeffs$beta1  * X1 +
      coeffs$beta2  * X2 +
      coeffs$beta3  * X3 +
      coeffs$beta4  * X4 +
      coeffs$gamma1 * (X2^2) +
      coeffs$gamma2 * sin(X3) +
      coeffs$gamma3 * as.numeric(X4 > 1) +
      coeffs$gamma4 * X1 * X5
  })

  if (!is.null(target_cor)) {
    1 <- var(S)
    residual_error_sd <- sqrt(var_S * (1 / target_cor^2 - 1))
  }

  if (is.null(residual_error_sd)) {
    residual_error_sd <- 1  # some default if neither is supplied
  }

  # ----- Generate outcome under coherent DGM -----
  residual_error <- rnorm(n, mean = 0, sd = residual_error_sd)
  Y <- S + treatment_effect * df$XTreat + residual_error

  return(Y)
}




hist1$y_fixed_progn <- outcome_from_covs_scaled_residuals(hist1, coeffs=coeff$hist1, residual_error_sd=1, target_cor = NULL)
hist2$y_fixed_progn <- outcome_from_covs_scaled_residuals(hist2, coeffs=coeff$hist2, residual_error_sd=1, target_cor = 0.8)


fit <- lm(
  y_fixed_progn ~ X1 + X2 + X3 + X4 + I(X2^2) + sin(X3) + I(X4 > 1) + X1:X5,
  data = hist1
)

summary(fit)$r.squared
fitted_vals <- predict(fit)
cor(fitted_vals, hist1$y_fixed_progn)
fitted_vals_hist2 <- predict(fit, newdata = hist2)
cor(fitted_vals_hist2, hist2$y_fixed_progn)



hist1$y_adjusted_sd <- outcome_from_covs_adjusted_error_sd(hist1, coeffs=coeff$hist1, residual_error_sd=1, target_cor = NULL)
hist2$y_adjusted_sd <- outcome_from_covs_adjusted_error_sd(hist2, coeffs=coeff$hist2, residual_error_sd=1, target_cor = 0.8)



fit <- lm(
  y_adjusted_sd ~ X1 + X2 + X3 + X4 + I(X2^2) + sin(X3) + I(X4 > 1) + X1:X5,
  data = hist1
)

summary(fit)$r.squared
fitted_vals <- predict(fit)
cor(fitted_vals, hist1$y_adjusted_sd)
fitted_vals_hist2 <- predict(fit, newdata = hist2)
cor(fitted_vals_hist2, hist2$y_adjusted_sd)



fit <- lm(
  y_fixed_progn ~ X1 + X2 + X3 + X4 + I(X2^2) + sin(X3) + I(X4 > 1) + X1:X5,
  data = hist1
)

summary(fit)$r.squared
fitted_vals <- predict(fit)
cor(fitted_vals, hist1$y_fixed_progn)
fitted_vals_hist2 <- predict(fit, newdata = hist2)
cor(fitted_vals_hist2, hist2$y_fixed_progn)

current$y <- outcome_from_covs(current, n=sample_size$current, coeffs=coeff$current, residual_error_sd, treatment_effect)

all_trials <- rbind( cbind(hist1, trial= "hist1", current=0), cbind(hist2, trial= "hist2", current=0), cbind(current, trial= "current", current=1))

ggplot(all_trials, aes(x=X1))+geom_histogram()+
  facet_wrap(~trial)
ggplot(all_trials, aes(x=X2))+geom_histogram()+
  facet_wrap(~trial)
ggplot(all_trials, aes(x=X3))+geom_histogram()+
  facet_wrap(~trial)
ggplot(all_trials, aes(x=X4))+geom_histogram()+
  facet_wrap(~trial)

ggplot(all_trials, aes(x=X1, y=y, color=XTreat))+
  geom_point()+
  facet_wrap(~trial)

ggplot(all_trials, aes(x=X2, y=y, color=XTreat))+
  geom_point()+
  facet_wrap(~trial)

ggplot(all_trials, aes(x=X3, y=y, color=XTreat))+
  geom_point()+
  facet_wrap(~trial)

ggplot(all_trials, aes(x=X4, y=y, color=XTreat))+
  geom_point()+
  facet_wrap(~trial)

ggplot(all_trials, aes(x=y))+geom_histogram()+
  facet_wrap(~trial)


