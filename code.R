library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(brms)
library(bayesplot)
library(loo)
library(posterior)
library(knitr)
library(kableExtra)
library(stringr)
data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")
set.seed(121212)

# Data description
###################################################################################

# Rename variables and add new variables
data <- data %>%
  rename(
    sex = Sex,
    beta = `Beta60 (g/kg/h)`,
    age = `Age (years)`,
    Co = `Co (g/Kg)`,
    weight = `Weight (kg)`,
    height = `Height (cm)`,
    AAC = `Amount of Alcohol Consumed (g)`,
    drinkingtime = `Drinking Time (h)`,
    maxBAC = `Maximum BAC (g/Kg)`,
    BACpeaktime = `BAC Peak time (min)`) %>%
  mutate(sex = factor(sex),
         Vd = AAC / (weight * Co), 
         TBW = ifelse(
           sex == "Male",
           2.447 - 0.09516 * age + 0.1074 * height + 0.3362 * weight,   # Male formula
           -2.097 + 0.1069 * height + 0.2466 * weight),
         T_Vd = TBW / weight,
         TV = ifelse(
           sex == "Male",
           1.40864*T_Vd,   # Male formula
           1.29720*T_Vd),
         BMI = weight/(height/100)**2,
         BACpeaktime =  BACpeaktime/60 )                 

data$beta <- abs(data$beta)   # Remove the negative sign (the rate itself is positive)
data$beta <- log(data$beta)   # log-transform makes it closer to normal distribution

# Standardized variables
data <- data %>%
  mutate(
    age_s = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
    weight_s = (weight - mean(weight, na.rm = TRUE)) / sd(weight, na.rm = TRUE),
    height_s = (height - mean(height, na.rm = TRUE)) / sd(height, na.rm = TRUE),
    AAC_s = (AAC - mean(AAC, na.rm = TRUE)) / sd(AAC, na.rm = TRUE),
    drinkingtime_s = (drinkingtime - mean(drinkingtime, na.rm = TRUE)) / sd(drinkingtime, na.rm = TRUE),
    maxBAC_s = (maxBAC - mean(maxBAC, na.rm = TRUE)) / sd(maxBAC, na.rm = TRUE),
    BACpeaktime_s = (BACpeaktime - mean(BACpeaktime, na.rm = TRUE)) / sd(BACpeaktime, na.rm = TRUE))


# Task 1
###################################################################################

# Model for log(beta)
f_beta <- bf(beta ~ 0+ sex + weight_s + height_s + drinkingtime_s)

# Set priors
prior <- c(set_prior("normal(0, 0.5)", class = "b"),
           set_prior("normal(0, 2)", class = "b", coef = "sexmale"),
           set_prior("normal(0, 2)", class = "b", coef = "sexfemale"),
           set_prior("exponential(1)", class = "sigma"))

# Fit Bayesian model for log(beta)
fit_single <- brm(formula = f_beta,
             data = data,
             family = gaussian(),
             prior = prior,
             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

# Task 2
###################################################################################

# Simpler model example
f_beta_example <- bf(beta ~ 0+ sex + weight_s + height_s)

fit_example <- brm(formula = f_beta_example,
             data = data,
             family = gaussian(),
             prior = prior,
             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

# Construct a new subject for prediction
new <- tibble(
  sex = factor("female", levels = levels(data$sex)),
  age = 70,
  weight = 70,
  height = 160
) %>%
  # Standardize using original sample means/SDs
  mutate(
    age_s = (age - mean(data$age, na.rm = TRUE)) / sd(data$age, na.rm = TRUE),
    weight_s = (weight - mean(data$weight, na.rm = TRUE)) / sd(data$weight, na.rm = TRUE),
    height_s = (height - mean(data$height, na.rm = TRUE)) / sd(data$height, na.rm = TRUE),
    sex = 'female',
    TBW = ifelse(
      sex == "Male",
      2.447 - 0.09516 * age + 0.1074 * height + 0.3362 * weight,   # Male formula
      -2.097 + 0.1069 * height + 0.2466 * weight),
    T_Vd = TBW / weight
    )

# Posterior predictive draws for log(beta)
beta_draws <- posterior_predict(fit_example, newdata = new, draws = 12000)

Ct <- 0.15   # current BAC
t <- 2       # time interval (hours)

# exp(beta_draws) converts log(beta) back to elimination rate (g/kg/h)
# C0 = Ct + beta * t gives estimated BAC at t hours earlier
C0_draws <- Ct + exp(beta_draws) * t

# Probability that earlier BAC exceeds 0.47 g/kg
P_over <- mean(C0_draws > 0.47)        

# Task 3
###################################################################################

# Model formulas
f_beta <- bf(beta ~ 0+ sex + weight_s + height_s + drinkingtime_s)
f_Vd<- bf(Vd ~ 0 + T_Vd:sex)

# Set priors
priors_joint <- c(
  set_prior("normal(0, 0.5)", class = "b", resp = "beta"),
  set_prior("normal(0, 2)", class = "b", coef = "sexmale", resp = "beta"),
  set_prior("normal(0, 2)", class = "b", coef = "sexfemale", resp = "beta"),
  set_prior("exponential(1)", class = "sigma", resp = "beta"),
  set_prior("normal(1.49, 2)", class = "b", resp = "Vd", coef = "T_Vd:sexfemale"),
  set_prior("normal(1.51, 2)", class = "b", resp = "Vd", coef = "T_Vd:sexmale"),
  set_prior("exponential(1)", class = "sigma", resp = "Vd")
)

# Joint model for beta and Vd
fit_joint <- brm(
  formula = f_beta + f_Vd + set_rescor(TRUE),
  data = data,
  prior = priors_joint,
  chains = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  seed = 2025,
  save_pars = save_pars(all = TRUE)
)


# Single model test
###################################################################################

# Posterior predictive draws
pp_draws <- exp(posterior_predict(fit_single, ndraws = 12000))

n_draws <- nrow(pp_draws)
n_obs <- ncol(pp_draws)         

# Predictive intervals
alpha_lo_95 <- 0.025; alpha_hi_95 <- 0.975
alpha_lo_50 <- 0.25;  alpha_hi_50 <- 0.75

# Summaries
pp_summary <- tibble(
  obs = seq_len(n_obs),
  beta_obs = exp(data$beta)  
) %>%
  mutate(
    pred_mean = colMeans(pp_draws),
    pred_median = apply(pp_draws, 2, median),
    # 95% PI
    PI95_low = apply(pp_draws, 2, quantile, probs = alpha_lo_95),
    PI95_high = apply(pp_draws, 2, quantile, probs = alpha_hi_95),
    # 50% PI
    PI50_low = apply(pp_draws, 2, quantile, probs = alpha_lo_50),
    PI50_high = apply(pp_draws, 2, quantile, probs = alpha_hi_50),
    # Errors for point predictions
    abs_err_mean = abs(pred_mean - beta_obs),
    sq_err_mean = (pred_mean - beta_obs)^2,
    # PIT value
    PIT = sapply(1:n_obs, function(j) mean(pp_draws[, j] <= beta_obs[j]))
  )

# Coverage rate
covered_95 <- sum(pp_summary$beta_obs >= pp_summary$PI95_low & pp_summary$beta_obs <= pp_summary$PI95_high, na.rm = TRUE)
covered_50 <- sum(pp_summary$beta_obs >= pp_summary$PI50_low & pp_summary$beta_obs <= pp_summary$PI50_high, na.rm = TRUE)

prop_95 <- covered_95 / n_obs
prop_50 <- covered_50 / n_obs

# Point-prediction summary (MAE, RMSE)
MAE <- mean(pp_summary$abs_err_mean, na.rm = TRUE)
RMSE <- sqrt(mean(pp_summary$sq_err_mean, na.rm = TRUE))


# Joint model test
###################################################################################

# Posterior predictive draws
pp_beta <- posterior_predict(fit_joint, resp = "beta", ndraws = 12000)
pp_Vd <- posterior_predict(fit_joint, resp = "Vd", ndraws = 12000)

# Posterior means
Vd_mean <- colMeans(pp_Vd)
beta_mean <- colMeans(pp_beta)
beta_mean <-exp(beta_mean)

# Compute C0 
pp_C0 = (data$AAC)/(data$weight * pp_Vd)
C0_mean <- colMeans(pp_C0)

# Compare with single-model beta predictions
beta_A <- posterior_predict(fit_single, ndraws = 12000)
beta_Amean <- colMeans(pp_beta)
beta_Amean <- exp(beta_Amean)


# Model Comparison
###################################################################################

# Student-t priors
prior_A <- c(set_prior("normal(0, 0.5)", class = "b"),
             set_prior("student_t(3, 0, 2", class = "b", coef = "sexmale"),
             set_prior("student_t(3, 0, 2", class = "b", coef = "sexfemale"),
             set_prior("exponential(1)", class = "sigma"))

# Normal(0,1) priors 
prior_B <- c(set_prior("normal(0, 0.01)", class = "b"),
             set_prior("normal(0, 1)", class = "b", coef = "sexmale"),
             set_prior("normal(0, 1)", class = "b", coef = "sexfemale"),
             set_prior("exponential(1)", class = "sigma"))

# Normal(0,5) priors 
fit_A <- brm(formula = f_beta,
             data = data,
             family = gaussian(),
             prior = prior_A,
             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

# Fit models
fit_B <- brm(formula = f_beta,
             data = data,
             family = gaussian(),
             prior = prior_B,
             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

# LOO
loo_A <- loo(fit_A, moment_match = TRUE)
loo_B <- loo(fit_B, moment_match = TRUE)

#Plot area
###################################################################################

# Density of log(beta)
F1 <- ggplot(data, aes(x = beta)) +
  geom_density() +
  labs(title = "Density Plot", x = "log(Beta)", y = "Density")

###################################################################################

# Short run diagnostic model
fit_d <- brm(formula = f_beta,
             data = data,
             family = gaussian(),
             prior = prior,
             chains = 4, iter = 100, warmup = 30, control = list(adapt_delta = 0.98, max_treedepth = 15))

# Extract draws including warmup
all <- as_draws_df(fit_d, inc_warmup = TRUE)

# Trace plot for b_sexfemale
F2 <- ggplot(all, aes(x = .iteration,
                y = b_sexmale,
                color = factor(.chain))) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  labs(x = "Iteration (including warmup)",
       y = "b_sexfemale",
       color = "Chain") +
  theme_bw()

###################################################################################

# Table of priors
priors_table <- tribble(~'' , ~Prior,
                        "Regression coefficient for male ", "Normal(0, 2)",
                        "Regression coefficient for female", "Normal(0, 2)",
                        "Regression coefficients (others)", "Normal(0, 0.5)",
                        "Residual SD ", "Exponential(1)")

T1 <- knitr::kable(priors_table, 
                   caption = "Weakly informative priors")

###################################################################################

# Posterior density overlay
F3 <- ggplot(data, aes(x = sex, y = beta)) + geom_boxplot() + ggtitle("beta by sex")

###################################################################################

# Trace plots
F4 <- mcmc_plot(fit_single, type = "dens_overlay")

###################################################################################

# MCMC Trace
F5 <- mcmc_trace(fit_single, 
                 pars = c("b_sexfemale", "b_sexmale",
                          "b_weight_s", "b_height_s",
                          "b_drinkingtime_s", "sigma"))

###################################################################################

# Posterior summary table
s <- summary(fit_single)
tab_all <- rbind(as.data.frame(s$fixed),
                 as.data.frame(s$random$ID))
T2 <- knitr::kable(
      tab_all,
      digits = 3,
      caption = "Posterior Summary from BRM")

###################################################################################

# Posterior predictive density overlay
yrep_C <- exp(posterior_predict(fit_single, draws = 2000))
F6 <- ppc_dens_overlay(exp(data$beta), yrep_C[1:2000, ]) + ggtitle("PPC density")

###################################################################################

# LOO PIT QQ-plot
F7 <- pp_check(fit_single, type = "loo_pit_qq")

###################################################################################

# Model checking summary table
model_check_table <- tibble(
  Metric = c("95% predictive interval coverage",
             "50% predictive interval coverage",
             "Mean Absolute Error",
             "Root Mean Square Error"),
  Value = c(prop_95,
            prop_50,
            MAE,
            RMSE))

T3 <- knitr::kable(model_check_table,
                   digits = 3,
                   caption = "Testing table")

###################################################################################

# Posterior distribution of C0
df_C0 <- tibble(C0 = C0_draws)

F8 <- ggplot(df_C0, aes(x = C0)) +
  geom_density(fill = "grey", alpha = 0.4) +
  geom_vline(xintercept = 0.47, color = "black", linetype = "dashed") +
  labs(title = "Posterior Predictive Distribution of C_0",
       x = "C_0 (g/L)",
       y = "Density") +
  theme_minimal(base_size = 14)

###################################################################################
# LOO and K-fold
loo_single <- loo(fit_single, moment_match = TRUE)
kfc_single <- kfold(fit_single, K = 10)

# Extract LOO and K-fold
loo_vals <- loo_single$estimates
k_vals <- kfc_single$estimates

# CV results table
table_cv_single <- tibble(
  Method      = c("LOO", "K-fold"),
  elpd        = c(loo_vals["elpd_loo", "Estimate"],
                  k_vals["elpd_kfold", "Estimate"]),
  elpd_SE     = c(loo_vals["elpd_loo", "SE"],
                  k_vals["elpd_kfold", "SE"]),
  p           = c(loo_vals["p_loo", "Estimate"],
                  k_vals["p_kfold", "Estimate"]),
  p_SE        = c(loo_vals["p_loo", "SE"],
                  k_vals["p_kfold", "SE"]),
  IC          = c(loo_vals["looic", "Estimate"],
                  k_vals["kfoldic", "Estimate"]),
  IC_SE       = c(loo_vals["looic", "SE"],
                  k_vals["kfoldic", "SE"])
)

T4 <- knitr::kable(
  table_cv_single,
  digits = 3,
  caption = "Cross-validation Comparison (LOO vs 10-fold)"
) %>%
  kable_styling(full_width = FALSE)

###################################################################################

# Summary of C0
C0_mean <- mean(C0_draws)
C0_median <- median(C0_draws)
d <- density(C0_draws)
C0_mode <- d$x[which.max(d$y)]

C0_summary <- tibble(
  Statistic = c("Mean", "Median", "Mode", 
                "2.5%", "25%", "75%", "97.5%", "P(C_0 > 0.47)"),
  Value = c(
    C0_mean,
    C0_median,
    C0_mode,
    quantile(C0_draws, 0.025),
    quantile(C0_draws, 0.25),
    quantile(C0_draws, 0.75),
    quantile(C0_draws, 0.975),
    P_over
  )
)
T5 <- knitr::kable(C0_summary,
                   digits = 3,
                   caption = "Example C_0 results")

###################################################################################

# Correlation between beta and Vd
cor_val <- cor(data$beta, data$Vd, use = "complete.obs")
cor_test <- cor.test(data$beta, data$Vd)

F9 <- ggplot(data, aes(x = Vd, y = exp(beta))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", color = "blue", se = TRUE) +
  annotate("text",
           x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("r = ", round(cor_val, 3)),
           size = 5) +
  labs(x = "Vd",
       y = "Beta") +
  theme_minimal(base_size = 14)

###################################################################################

# Density of log(Vd)
F10 <- ggplot(data, aes(x = Vd)) +
  geom_density() +
  labs(title = "Density Plot", x = "V_d", y = "Density")

###################################################################################

# Summary from joint model
s <- summary(fit_joint)

# Fixed effects
fixed_tab <- as.data.frame(s$fixed)
fixed_tab$Parameter <- rownames(fixed_tab)
fixed_tab$Group <- ifelse(grepl("^beta_", fixed_tab$Parameter), "Beta", "Vd")
fixed_tab <- fixed_tab[, c("Parameter","Estimate","Est.Error","Rhat","Bulk_ESS","Tail_ESS","Group")]

# Sigma parameters
sigma_tab <- as.data.frame(s$spec_pars)
sigma_tab$Parameter <- rownames(sigma_tab)
sigma_tab$Group <- ifelse(grepl("beta", sigma_tab$Parameter), 
                          "Beta", "Vd")
sigma_tab <- sigma_tab[, c("Parameter","Estimate","Est.Error","Rhat","Bulk_ESS","Tail_ESS","Group")]

# Residual correlation
rescor_tab <- as.data.frame(s$rescor)
rescor_tab$Parameter <- rownames(rescor_tab)
rescor_tab$Group <- "Correlation"
rescor_tab <- rescor_tab[, c("Parameter","Estimate","Est.Error","Rhat","Bulk_ESS","Tail_ESS","Group")]

# Combine
tab_joint <- dplyr::bind_rows(fixed_tab, sigma_tab, rescor_tab)
tab_joint <- tibble::as_tibble(tab_joint)

T6 <- knitr::kable(
  tab_joint,
  digits = 4,
  caption = "Joint model results"
)

###################################################################################

# Joint posterior samples
beta_post <- exp(posterior_predict(fit_joint, resp = "beta", ndraws = 12000))
Vd_post   <- posterior_predict(fit_joint, resp = "Vd",   ndraws = 12000)

# One sample per draw 
beta_sample <- rowMeans(beta_post)   
Vd_sample   <- rowMeans(Vd_post)    

# Joint posterior density
df_joint <- tibble(
  beta = beta_sample,
  Vd   = Vd_sample
)


F11 <- ggplot(df_joint, aes(x = Vd, y = beta)) +
  geom_bin2d(bins = 40) +
  scale_fill_gradient(low = "white", high = "black") +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Posterior V_d",
    y = "Posterior Beta")

###################################################################################

# LOO and K-fold CV for joint model
loo_joint <- loo(fit_joint, moment_match = TRUE)
kfc <- kfold(fit_joint, K = 10)

# Extract LOO and K-fold
loo_vals <- loo_joint$estimates
k_vals <- kfc$estimates

# CV results table
table_cv <- tibble(
  Method      = c("LOO", "K-fold"),
  elpd        = c(loo_vals["elpd_loo", "Estimate"],
                  k_vals["elpd_kfold", "Estimate"]),
  elpd_SE     = c(loo_vals["elpd_loo", "SE"],
                  k_vals["elpd_kfold", "SE"]),
  p           = c(loo_vals["p_loo", "Estimate"],
                  k_vals["p_kfold", "Estimate"]),
  p_SE        = c(loo_vals["p_loo", "SE"],
                  k_vals["p_kfold", "SE"]),
  IC          = c(loo_vals["looic", "Estimate"],
                  k_vals["kfoldic", "Estimate"]),
  IC_SE       = c(loo_vals["looic", "SE"],
                  k_vals["kfoldic", "SE"])
)

T7 <- knitr::kable(
  table_cv,
  digits = 3,
  caption = "Cross-validation Comparison (LOO vs 10-fold)"
) %>%
  kable_styling(full_width = FALSE)

###################################################################################

# Prediction error metrics
MSE_Vd <- mean(abs(Vd_mean - data$Vd) , na.rm = TRUE )
MSE_beta<- mean(abs(beta_mean - exp(data$beta)) , na.rm = TRUE )
MSE_betaA <- mean(abs(beta_Amean - exp(data$beta)) , na.rm = TRUE )
pp_C0 = (data$AAC)/(data$weight * pp_Vd)
C0_mean <- colMeans(pp_C0)
MSE_C0<- mean(abs(C0_mean - data$Co) , na.rm = TRUE )

# Vd coverage rates
# 95%
Vd_low_95  <- apply(pp_Vd, 2, quantile, 0.025)
Vd_high_95 <- apply(pp_Vd, 2, quantile, 0.975)
Vd_cov_95  <- mean(data$Vd >= Vd_low_95 & data$Vd <= Vd_high_95)

# 50%
Vd_low_50  <- apply(pp_Vd, 2, quantile, 0.25)
Vd_high_50 <- apply(pp_Vd, 2, quantile, 0.75)
Vd_cov_50  <- mean(data$Vd >= Vd_low_50 & data$Vd <= Vd_high_50)

# C0 predictive matrix
pp_C0 <- (data$AAC) / (data$weight * pp_Vd)

# C0 coverage rates
# 95%
C0_low_95  <- apply(pp_C0, 2, quantile, 0.025)
C0_high_95 <- apply(pp_C0, 2, quantile, 0.975)
C0_cov_95  <- mean(data$Co >= C0_low_95 & data$Co <= C0_high_95)

# 50%
C0_low_50  <- apply(pp_C0, 2, quantile, 0.25)
C0_high_50 <- apply(pp_C0, 2, quantile, 0.75)
C0_cov_50  <- mean(data$Co >= C0_low_50 & data$Co <= C0_high_50)

# Summary table
table_mse <- tibble(
  Metric = c(
    "MSE_Vd",
    "MSE_beta (joint model)",
    "MSE_beta (model C)",
    "MSE_C0",
    "Vd Coverage (95%)",
    "Vd Coverage (50%)",
    "C0 Coverage (95%)",
    "C0 Coverage (50%)"
  ),
  Estimate = c(
    MSE_Vd,
    MSE_beta,
    MSE_betaA,
    MSE_C0,
    Vd_cov_95,
    Vd_cov_50,
    C0_cov_95,
    C0_cov_50
  )
)

T8 <- knitr::kable(
  table_mse,
  digits = 4,
  caption = "Prediction Error and Coverage Rate"
) %>%
  kable_styling(full_width = FALSE)

###################################################################################

# Comparison priors table
priors_comparison <- tribble(
  ~Parameter, ~Original,            ~Prior_A,             ~Prior_B,       
  
  "Regression coefficient: male",
  "Normal(0, 2)",               "Student-t(3, 0, 2)", "Normal(0, 1)", 
  
  "Regression coefficient: female",
  "Normal(0, 2)",               "Student-t(3, 0, 2)", "Normal(0, 1)",  
  
  "Regression coefficients (others)",
  "Normal(0, 0.5)",             "Normal(0, 0.5)",     "Normal(0, 0.01)", 
  
  "Residual SD (sigma)",
  "Exponential(1)",             "Exponential(1)",     "Exponential(1)"
)

T9 <-knitr::kable(priors_comparison,
                  caption = "Comparison priors with original prior")


###################################################################################

# LOO comparison across models
comp <- loo_compare(loo_A, loo_B, loo_single)

comp_tbl <- comp |> 
  as.data.frame() |> 
  tibble::rownames_to_column("Model") |> 
  dplyr::select(Model, elpd_diff, se_diff)

T10 <-knitr::kable(comp_tbl,
                   digits = 3,
                   caption = "LOO Comparison of different priors' model")

###################################################################################

# Compute posterior mean
Cbeta_A = exp(posterior_predict(fit_A))
Cbeta_B = exp(posterior_predict(fit_B))
Cbeta_A1 <- colMeans(Cbeta_A)
Cbeta_B1 <- colMeans(Cbeta_B)

# Plot mean posterior mean vs true beta
F12 <- ggplot(data) +
          geom_point(aes(x = exp(data$beta), y = Cbeta_B1, color = 'Prior B results')) +
          geom_point(aes(x = exp(data$beta), y = colMeans(pp_draws), color = 'Original prior results')) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
          theme_minimal(base_size = 14) +
          labs(x = "True Beta",
          y = "Posterior Mean Beta",
          title = "Posterior Mean vs True Beta")

###################################################################################

pandoc <- rmarkdown::pandoc_exec()
cmd <- sprintf('"%s" P2.Rmd --from=markdown --to=plain | wc -w', pandoc)
system(cmd)
