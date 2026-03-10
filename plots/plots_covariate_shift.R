source("utils/__init__.R")
source("utils/helper_functions.R")

tmp <- new.env()
load("results/mild_covariate_shift/simulation_results.RData", envir = tmp)

list2env(
  setNames(as.list(tmp), paste0(ls(tmp), "_mild_shift")),
  envir = .GlobalEnv
)


tmp <- new.env()
load("results/strong_covariate_shift/simulation_results.RData", envir = tmp)

list2env(
  setNames(as.list(tmp), paste0(ls(tmp), "_strong_shift")),
  envir = .GlobalEnv
)

tmp <- new.env()
load("results/no_covariate_shift/simulation_results.RData", envir = tmp)

list2env(
  setNames(as.list(tmp), paste0(ls(tmp), "_no_shift")),
  envir = .GlobalEnv
)


plot_df_combined<-rbind(cbind(plot_df_mild_shift, shift="mild"), cbind(plot_df_strong_shift, shift="strong"),  cbind(plot_df_no_shift, shift="no"))

plot_df <- plot_df_combined %>%
  group_by(method, model, shift) %>%
  summarise(
    bias   = mean(bias),
    MSE    = mean(se),
    reject = mean(reject),
    emp_var = var(est_mean),
    mc_sd = sd(est_mean)
  )%>%ungroup()%>%
  mutate(
    method = unlist(method),
    model  = factor(unlist(model), levels=c("NO_BORROWING", "ANOVA", "ANCOVA", "PROCOVA - superlearner", "PROCOVA - oracle", "A_LASSO", "A_LASSO_ADJ", "COMM_DB", "SAM_PRIOR", "FULL_BORROWING", "POOLING", "Bayesian_PROCOVA")),
  )






head(plot_df_combined)


ggplot(plot_df, aes(x=model, y=bias, shape=method, color=model))+
  geom_point()+
  geom_errorbar(
    aes(
      ymin = bias - mc_sd,
      ymax = bias + mc_sd
    ),
    width = 0.2, 
    size=0.2
  )+
  facet_wrap(~shift, nrow = 1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
