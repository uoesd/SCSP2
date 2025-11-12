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

data <- read_excel("SCS_BAC_and_BrAC_split_TOP.xlsx")

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
  mutate(sex = factor(sex))



num_summary <- data %>%
  summarise(
    n = n(),
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    q2.5_beta = quantile(beta, 0.025, na.rm = TRUE),
    q97.5_beta = quantile(beta, 0.975, na.rm = TRUE),
    min_beta = min(beta, na.rm = TRUE),
    max_beta = max(beta, na.rm = TRUE),
    prop_neg = mean(beta < 0, na.rm = TRUE)
  )

data <- data %>%
  mutate(
    age_s = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
    weight_s = (weight - mean(weight, na.rm = TRUE)) / sd(weight, na.rm = TRUE),
    height_s = (height - mean(height, na.rm = TRUE)) / sd(height, na.rm = TRUE),
    AAC_s = (AAC - mean(AAC, na.rm = TRUE)) / sd(AAC, na.rm = TRUE),
    drinkingtime_s = (drinkingtime - mean(drinkingtime, na.rm = TRUE)) / sd(drinkingtime, na.rm = TRUE),
    maxBAC_s = (maxBAC - mean(maxBAC, na.rm = TRUE)) / sd(maxBAC, na.rm = TRUE),
    BACpeaktime_s = (BACpeaktime - mean(BACpeaktime, na.rm = TRUE)) / sd(BACpeaktime, na.rm = TRUE)
  )
data$beta <- abs(data$beta)  # because Beta60 in file is negative of slope. :contentReference[oaicite:2]{index=2}

print(num_summary)
ggplot(data, aes(x = beta)) + geom_histogram(bins = 40) + ggtitle("Histogram of beta")
ggplot(data, aes(y = beta)) + geom_boxplot() + ggtitle("Boxplot of beta")
ggplot(data, aes(x = sex, y = beta)) + geom_boxplot() + ggtitle("beta by sex")
ggplot(data, aes(x = weight, y = beta)) + geom_point() + geom_smooth(method = "loess") + ggtitle("beta vs weight")

formula_main <- bf(beta ~ 1 + sex + age_s + weight_s + height_s)
emp_mean <- mean(data$beta, na.rm = TRUE)
emp_sd   <- sd(data$beta, na.rm = TRUE)
chains <- 2
iter <- 1000
warmup <- 200
control <- list(adapt_delta = 0.98, max_treedepth = 15)
seed <- 2025

# Student-t intercept prior (robust)
prior_A <- c(set_prior(paste0("student_t(3, ", signif(emp_mean, 3), ", ", signif(emp_sd * 2, 3), ")"), class = "Intercept"),
             set_prior("normal(0, 0.01)", class = "b"),
             set_prior("exponential(1)", class = "sigma"))


# Student-t likelihood prior includes nu
prior_B <- c(prior_A, set_prior("gamma(2, 0.1)", class = "nu"))


# Normal intercept prior (sensitivity)
prior_C <- c(set_prior(paste0("normal(", signif(emp_mean,3), ", ", signif(emp_sd * 2, 3), ")"), class = "Intercept"),
             set_prior("normal(0, 0.01)", class = "b"),
             set_prior("exponential(1)", class = "sigma"))



message("Fitting Model A (gaussian likelihood + student-t intercept prior)...")
fit_A <- brm(formula = formula_main,
             data = data,
             family = gaussian(),
             prior = prior_A,
             chains = chains, iter = iter, warmup = warmup, control = control, seed = seed)

message("Fitting Model B (student-t likelihood + student-t prior)...")
fit_B <- brm(formula = formula_main,
             data = data,
             family = student(),
             prior = prior_B,
             chains = chains, iter = iter, warmup = warmup, control = control, seed = seed)

message("Fitting Model C (gaussian likelihood + normal intercept prior - sensitivity)...")
fit_C <- brm(formula = formula_main,
             data = data,
             family = gaussian(),
             prior = prior_C,
             chains = chains, iter = iter, warmup = warmup, control = control, seed = seed)

priors <- c(
  prior(normal(0,2), class = "b"),
  prior(cauchy(0,2), class = "sigma", lb = 0),
  prior(constant(8), class = "nu"))

fit <- brm(formula_main, data = data, 
           prior = priors,
           family = student(),
           chains = 2, iter = 1000, warmup = 200, seed = 42,
           control = list(adapt_delta = 0.98))

print(summary(fit_A))
print(summary(fit_B))
print(summary(fit_C))
#plot(fit_A)
#plot(fit_B)install.packages("tinytex")
tinytex::install_tinytex()   # installs TinyTeX (no admin usually required)

#plot(fit_C)

yrep_A <- posterior_predict(fit_A, draws = 200)
ppc_dens_overlay(data$beta, yrep_A[1:200, ]) + ggtitle("PPC density - Model A")
yrep_B <- posterior_predict(fit_B, draws = 200)
ppc_dens_overlay(data$beta, yrep_B[1:200, ]) + ggtitle("PPC density - Model B")
yrep_C <- posterior_predict(fit_C, draws = 200)
ppc_dens_overlay(data$beta, yrep_C[1:200, ]) + ggtitle("PPC density - Model C")

loo_A <- loo(fit_A, moment_match = TRUE)
loo_B <- loo(fit_B, moment_match = TRUE)
loo_C <- loo(fit_C, moment_match = TRUE)
print(loo_compare(loo_A, loo_B, loo_C))


new <- tibble(
  sex = factor("female", levels = levels(data$sex)),
  age = 70,
  weight = 70,
  height = 160
) %>%
  mutate(
    age_s = (age - mean(data$age, na.rm = TRUE)) / sd(data$age, na.rm = TRUE),
    weight_s = (weight - mean(data$weight, na.rm = TRUE)) / sd(data$weight, na.rm = TRUE),
    height_s = (height - mean(data$height, na.rm = TRUE)) / sd(data$height, na.rm = TRUE),
    sex = 'female'
  )

beta_draws <- posterior_predict(fit_A, newdata = new, draws = 1600)
Ct <- 0.15   
t <- 2      
C0_draws <- Ct + beta_draws * t
P_over <- mean(C0_draws > 0.47)        
quantile(C0_draws, prob = c(0.025,0.25, 0.5, 0.975))
P_over
# empirical 2.5th percentile of beta distribution:
emp_beta_2.5 <- quantile(data$beta, probs = 0.025, na.rm = TRUE)
quantile(C0_draws, probs = 0.025, na.rm = TRUE)

# Evaluate predictive performance on training set (observed beta)
# Assumes: 'data' is your dataframe (with observed beta) and 'fit_A' is your brms fit object
library(brms)
library(posterior)
library(dplyr)
library(ggplot2)

# 1) Draw posterior predictive samples for observed rows
# posterior_predict returns matrix draws x observations (includes observation noise)
pp_draws <- posterior_predict(fit_B, ndraws = 1600)  # 2000 draws is enough; increase if needed
# posterior_epred returns expected mean draws (no obs noise) if you prefer
epred_draws <- posterior_epred(fit_B, ndraws = 1600)

n_draws <- nrow(pp_draws)
n_obs <- ncol(pp_draws)           # should be 100 in your case

# 2) For each observation compute 50% and 95% predictive intervals (from pp_draws)
alpha_lo_95 <- 0.025; alpha_hi_95 <- 0.975
alpha_lo_50 <- 0.25;  alpha_hi_50 <- 0.75

pp_summary <- tibble(
  obs = seq_len(n_obs),
  beta_obs = data$beta  # observed true beta for each row
) %>%
  mutate(
    pred_mean = colMeans(epred_draws),
    pred_median = apply(epred_draws, 2, median),
    # 95% PI from posterior_predict (includes obs noise)
    PI95_low = apply(pp_draws, 2, quantile, probs = alpha_lo_95),
    PI95_high = apply(pp_draws, 2, quantile, probs = alpha_hi_95),
    # 50% PI
    PI50_low = apply(pp_draws, 2, quantile, probs = alpha_lo_50),
    PI50_high = apply(pp_draws, 2, quantile, probs = alpha_hi_50),
    # point error metrics
    abs_err_mean = abs(pred_mean - beta_obs),
    sq_err_mean = (pred_mean - beta_obs)^2,
    # PIT value (proportion of predictive draws <= observed)
    PIT = sapply(1:n_obs, function(j) mean(pp_draws[, j] <= beta_obs[j]))
  )

# 3) Coverage counts
covered_95 <- sum(pp_summary$beta_obs >= pp_summary$PI95_low & pp_summary$beta_obs <= pp_summary$PI95_high, na.rm = TRUE)
covered_50 <- sum(pp_summary$beta_obs >= pp_summary$PI50_low & pp_summary$beta_obs <= pp_summary$PI50_high, na.rm = TRUE)

prop_95 <- covered_95 / n_obs
prop_50 <- covered_50 / n_obs

# 4) Point-prediction summary (MAE, RMSE) using posterior mean
MAE <- mean(pp_summary$abs_err_mean, na.rm = TRUE)
RMSE <- sqrt(mean(pp_summary$sq_err_mean, na.rm = TRUE))

# 5) Bayesian R^2 (using brms helper)
r2_draws <- as.numeric(bayes_R2(fit_A))  # ensures numeric vector
r2_median <- median(r2_draws)
r2_CI <- quantile(r2_draws, c(0.025, 0.975))

# 6) Print results
cat("Observations (n):", n_obs, "\n")
cat("95% predictive interval coverage:", covered_95, "/", n_obs, " = ", round(100*prop_95,2), "%\n")
cat("50% predictive interval coverage:", covered_50, "/", n_obs, " = ", round(100*prop_50,2), "%\n")
cat("MAE (using posterior mean):", signif(MAE,4), "\n")
cat("RMSE (using posterior mean):", signif(RMSE,4), "\n")
cat("Bayesian R^2 median:", signif(r2_median,4), " (95% CI:", signif(r2_CI[1],4), "-", signif(r2_CI[2],4), ")\n")

# 7) PIT histogram (should be flat if predictive distribution calibrated)
ggplot(pp_summary, aes(x = PIT)) +
  geom_histogram(bins = 20, color = "black", fill = "skyblue") +
  labs(title = "PIT histogram (calibration): ideally flat", x = "PIT value", y = "count")

# 8) Optional: list the observation indices that fall outside 95% PI
outliers_95_idx <- pp_summary %>% filter(!(beta_obs >= PI95_low & beta_obs <= PI95_high)) %>% pull(obs)
cat("Indices outside 95% PI (count):", length(outliers_95_idx), "\n")
if (length(outliers_95_idx) > 0) cat("Examples:", head(outliers_95_idx, 20), "\n")

# 9) If you want a flexible "how many are correctly predicted" with tolerance:
#    Define tolerance tol (absolute error), e.g., tol = 0.01 g/kg/h
tol <- 0.01
n_within_tol <- sum(pp_summary$abs_err_mean <= tol, na.rm = TRUE)
cat("Number of obs with abs error <= ", tol, ": ", n_within_tol, "/", n_obs,
    " (", round(100 * n_within_tol / n_obs,1), "% )\n", sep = "")

lm(beta~1, data)
summary(lm(beta~1 + sex + weight_s + age_s+height_s, data))
# Example: after fitting fit_A
mcmc_trace(as.array(fit_A), pars = c("b_Intercept", "b_weight_s"))






###
library(knitr)

#Before fitting the Bayesian regression, appropriate priors were specified for the intercept, regression coefficients, and residual standard deviation.
#We use weak information priors, allowing the observed data to dominate inference process.
#The chosen priors are summarised below:

priors_table <- tribble(
  ~Parameter, ~Prior, ~Role,
  "Intercept", "Normal(0.184, 0.066)", "Average alcohol elimination rate",
  "Coefficients (b)", "Normal(0, 0.01)", "Covariate effects",
  "Residual SD (σ)", "Exponential(1)", "Random variation")

kable(priors_table,
      align = "lcl")

ggplot(data, aes(x = beta)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = emp_mean, sd = emp_sd * 2), color = "red", size = 1.2) +
  labs(title = "Empirical β distribution vs Normal prior for intercept",
       x = expression(beta~" elimination rate (g/kg/h)"),
       y = "Density") +
  theme_minimal()

#The intercept prior was specified as a Normal(0.184, 0.066) distribution, 
#where the mean corresponds to the empirical average of β estimated from the observed data
#and the standard deviation (0.066 ≈ 2 × empirical SD) provides sufficient uncertainty to cover the plausible physiological range.
#As shown in Figure, the red curve (prior) fits well with the empirical histogram, this ensures the prior reflects realistic biology without constraining the data.

#For the regression coefficients, a Normal(0, 0.01) prior was used.
#This setting assumes that each predictor may have a limited but non-zero influence on β after standardisation.
#This approach follows recommendations from Gelman et al. (2008)

#Reference:
#Gelman, A., Jakulin, A., Pittau, M. G., & Su, Y.-S. (2008).
#A weakly informative default prior distribution for logistic and other regression models.
#Annals of Applied Statistics, 2(4), 1360–1383.
#https://doi.org/10.1214/08-AOAS191

#The residual standard deviation σ was assigned an Exponential(1) prior.
#This prior keeps σ positive and allows the model to capture unexplained variation in β.
#This choice is consistent with the guidance of Gelman (2006) and Gelman et al. (2020).

#References:
#Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models.
#Bayesian Analysis, 1(3), 515–533.
#https://doi.org/10.1214/06-BA117A

#Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B., Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020). 
#Bayesian workflow. arXiv preprint arXiv:2011.01808. 
#https://arxiv.org/abs/2011.01808


r2_A <- median(as.numeric(bayes_R2(fit_A)))
r2_B <- median(as.numeric(bayes_R2(fit_B)))
r2_C <- median(as.numeric(bayes_R2(fit_C)))

model_comp_final <- tribble(
  ~Model, ~`elpd_diff`, ~`se_diff`, ~`Bayesian R²`,
  "A",  0.0,  0.0, 0.1315,
  "B", -0.7,  0.3, 0.1307,
  "C",  0.0,  0.0, 0.1316
)

kable(model_comp_final,
      caption = "Model comparison based on LOO and Bayesian R²",
      align = "lrrr",
      digits = 4)
