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

formula_main <- bf(beta ~ 1 + sex + age_s + weight_s + height_s + AAC_s + drinkingtime_s + maxBAC_s + BACpeaktime_s)
emp_mean <- mean(data$beta, na.rm = TRUE)
emp_sd   <- sd(data$beta, na.rm = TRUE)
chains <- 4
iter <- 4000
warmup <- 1500
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
print(summary(fit_A))
print(summary(fit_B))
print(summary(fit_C))
plot(fit_A)
plot(fit_B)
plot(fit_C)

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

# formula: change covariates as appropriate
f <- bf(Beta60 ~ 1 + sex + age + weight)

priors <- c(
  prior(normal(0,2), class = "b"),
  prior(cauchy(0,2), class = "sigma", lb = 0),
  prior(constant(8), class = "nu"))

fit <- brm(f, data = data, 
           prior = priors,
           family = student(),
           chains = 4, iter = 4000, warmup = 1000, seed = 42,
           control = list(adapt_delta = 0.98))

fit
plot(fit)
new <fitnew <- data.frame(sex= 'female', age = 70,
                  weight = 70)

new2 <- data.frame(sex= data$sex, age = data$age,
                   weight = data$weight)

beta_draws <- posterior_predict(fit, newdata = new, draws = 2000)
beta_draws2 <- posterior_predict(fit, newdata = new2, draws = 2)

Ct <- 0.15   
t <- 2      
C0_draws <- Ct + beta_draws * t



P_over <- mean(C0_draws > 0.47)        
quantile(C0_draws, prob = c(0.025,0.25, 0.5, 0.975))
P_over


# empirical 2.5th percentile of beta distribution:
emp_beta_2.5 <- quantile(dat$beta, probs = 0.025, na.rm = TRUE)
emp_beta_97.5 <- quantile(dat$beta, probs = 0.975, na.rm = TRUE)

lm(beta~1+sex+age_s, data)
summary(lm(beta~1 + sex + weight_s + age_s, data))
