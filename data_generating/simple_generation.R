
source("utils/helper_functions.R")
source("utils/__init__.R")



### configuration 
# start with a base setting - all datasets have the same covariate distribution


covars <- list(
  hist1  = list(mu = c(X1= -2, X2=0, X3= 0), sd = c(X1=1, X2=1, X3=1), pX4 = 0.5),
  hist2  = list(mu = c(X1= -2, X2=0, X3= 0), sd = c(X1=1, X2=1, X3=1), pX4 = 0.5),
  current = list(mu = c(X1= 0, X2=0, X3= 0), sd = c(X1=1, X2=1, X3=1), pX4 = 0.5)
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
  hist1  = list(beta_lin=1, beta_quad=-1, beta_bell=0, beta_sin=2, beta_interaction=0.5),
  hist2  = list(beta_lin=1, beta_quad=-1, beta_bell=0, beta_sin=2, beta_interaction=0.5),
  current = list(beta_lin=1, beta_quad=-1, beta_bell=0, beta_sin=2, beta_interaction=0.5)
)


gen_covs <- function(n, covars) {
  X1 <- rnorm(n, covars$mu["X1"], covars$sd["X1"])
  X2 <- rnorm(n, covars$mu["X2"], covars$sd["X2"])
  X3 <- rnorm(n, covars$mu["X3"], covars$sd["X3"])
  X4  <- rbinom(n, 1, covars$pX4)
  data.frame(X1, X2, X3, X4)
}

hist1 <- gen_covs(sample_size$hist1, covars$hist1) 
hist2 <- gen_covs(sample_size$hist2, covars$hist2) 
current <- gen_covs(sample_size$current, covars$current) 


hist1$XTreat <- rbinom(sample_size$hist1, 1, 0)
hist2$XTreat <- rbinom(sample_size$hist2, 1, 0)
current$XTreat <- rbinom(sample_size$current, 1, randomization_ratio)



outcome_from_covs <- function(df, n, coeffs, residual_error_sd, treatment_effect=0) {
    lin  <- coeffs$beta_lin  * df$X1
    quad <- coeffs$beta_quad * (df$X2^2)
    #bell <- coeffs$beta_bell * exp(- (dt$X3 - bell_c)^2 / (2 * bell_s^2))
    sinus <- coeffs$beta_sin *  sin(df$X3)
    interaction <- coeffs$beta_interaction * df$X4
    base <- lin + quad + interaction + sinus
    treat <- treatment_effect * df$XTreat
    residual_error <- rnorm(n, 0, residual_error_sd)
    return(base + treat + residual_error)
}

hist1$y <- outcome_from_covs(hist1, n=sample_size$hist1, coeffs=coeff$hist1, residual_error_sd)
hist2$y <- outcome_from_covs(hist2, n=sample_size$hist2, coeffs=coeff$hist2, residual_error_sd)
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







# PROCOVA

## create simple prognostic model (linear regression)
linear_progn<-lm(y ~ X1 + X2 + X3 + X4, data=all_trials, subset=(current==0))
summary(linear_progn)
current$lin_progn_score<-predict(linear_progn, newdata = current)



#install.packages("SuperLearner")
#install.packages("glmnet")     # for elastic net / lasso
#install.packages("ranger")     # for random forest
#install.packages("xgboost")    # for gradient boosting
library(SuperLearner)
library(glmnet)
library(ranger)
library(xgboost)


train_data <- subset(all_trials, current == 0)
y <- train_data$y
X <- subset(train_data, select = c(X1, X2, X3, X4))

SL.library <- c("SL.glm",        # standard linear model
                "SL.glmnet",     # elastic net
                "SL.ranger",     # random forest
                "SL.xgboost")    # gradient boosting


set.seed(123)
super_model <- SuperLearner(
  Y = y,
  X = X,
  SL.library = SL.library,
  family = gaussian(),   # use binomial() if y is binary
  method = "method.NNLS" # non-negative least squares ensemble
)

summary(super_model)

X_new <- subset(all_trials, current == 1, select = c(X1, X2, X3, X4))
current$super_progn_score<-predict(super_model, newdata = X_new)$pred


lm(y ~ X1 + X2 + X3 + X4, data=current)%>% summary()

lm(y ~ XTreat, data=current)%>% summary()

lm(y ~ XTreat + lin_progn_score, data=current)%>% summary()
lm(y ~ XTreat + super_progn_score, data=current)%>% summary()

lm(y ~ XTreat + X1 + X2 + X3 + X4, data=current)%>% summary()

lm(y ~ XTreat + super_progn_score, data=current)%>% summary()
lm(y ~ XTreat + X1 + X2 + X3 + X4 + super_progn_score, data=current)%>% summary()


# BDB

## adaptive lasso (Li et al.)

# external trial indicator (B_i in the paper)

  gamma <- 0.1
  lambda <- 0.05

  all_trials$historical<-1-all_trials$current

  inital_model<-lm(y ~ XTreat + historical, data=all_trials, x=TRUE)

  summary(inital_model)

  adaptive_weights <- (abs(coef(inital_model)))^-gamma
  adaptive_weights[1:2] <- 0
  #adaptive_weights[4:7] <- 0


  adaptive_weights <- adaptive_weights[names(adaptive_weights) != "(Intercept)"]


  print(adaptive_weights)
  scaling_factor <- sum(adaptive_weights) / length(adaptive_weights)
  adjusted_adaptive_weights <- adaptive_weights / scaling_factor
  lambda_db_model <- lambda * scaling_factor
    

  predictor_matrix <- inital_model$x
  predictor_matrix <- predictor_matrix[, colnames(predictor_matrix) != "(Intercept)"]

  model.alasso <- glmnet(
        x = predictor_matrix,
        y = all_trials$y,
        family = "gaussian",
        penalty.factor = adjusted_adaptive_weights,
        alpha = 1,
        lambda = lambda_db_model
      )


  coeffs.alasso <- as.matrix(coef(model.alasso))
  print(coeffs.alasso)




## robust MAP 

# install.packages("RBesT")
library(RBesT)

# install.packages("brms")
library(brms)

# Fit a hierarchical model:
# - treatment effect estimated in current trial
# - historical controls inform the prior on control mean (via partial pooling)

# Restrict to control patients for the hierarchical prior on control means
hist_controls <- subset(all_trials, current == 0)

# Hierarchical model for control means across trials
fit_control <- brm(
  y ~ 1 + (1 | trial),   # random intercept per trial (partial pooling)
  data = hist_controls,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "Intercept"),     # vague population mean prior
    prior(exponential(1), class = "sd")            # half-normal(0,1) for between-trial SD
  ),
  chains = 4, iter = 4000, cores = 4
)

# Extract the posterior for the population-level control mean
summary(fit_control)

# Now, in the CURRENT trial, model treatment vs control using that prior as an informative prior
mu_map <- fixef(fit_control)["Intercept", "Estimate"]
sd_map <- fixef(fit_control)["Intercept", "Est.Error"]

# Hierarchical + MAP model for the current trial
fit_current <- brm(
  y ~ XTreat, data = subset(all_trials, current == 1),
  family = gaussian(),
  prior = c(
    prior(normal(mu_map, sd_map), class = "Intercept"),   # MAP prior from historical controls
    prior(normal(0, 10), class = "b")                     # vague for treatment effect
  ),
  chains = 4, iter = 4000, cores = 4
)

summary(fit_current)
