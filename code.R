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
data$beta <- abs(data$beta) 
data$beta <- log(data$beta) 


data <- data %>%
  mutate(
    age_s = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
    weight_s = (weight - mean(weight, na.rm = TRUE)) / sd(weight, na.rm = TRUE),
    height_s = (height - mean(height, na.rm = TRUE)) / sd(height, na.rm = TRUE),
    AAC_s = (AAC - mean(AAC, na.rm = TRUE)) / sd(AAC, na.rm = TRUE),
    drinkingtime_s = (drinkingtime - mean(drinkingtime, na.rm = TRUE)) / sd(drinkingtime, na.rm = TRUE),
    maxBAC_s = (maxBAC - mean(maxBAC, na.rm = TRUE)) / sd(maxBAC, na.rm = TRUE),
    BACpeaktime_s = (BACpeaktime - mean(BACpeaktime, na.rm = TRUE)) / sd(BACpeaktime, na.rm = TRUE))

cov_xy <- cov(data$Vd, data$beta)
cor_xy <- cor(data$Vd, data$beta)
cor.test(data$Vd, data$beta)

model <- lm(beta~ 0+ sex+ weight_s + height_s + drinkingtime_s, data)
summary(model)


f_beta <- bf(beta ~ 0+ sex + weight_s + height_s + drinkingtime_s)

# Student-t intercept prior (robust)
prior_A <- c(set_prior("normal(0, 0.01)", class = "b"),
             set_prior("student_t(3, 0.15, 0.02", class = "b", coef = "sexmale"),
             set_prior("student_t(3, 0.18, 0.02", class = "b", coef = "sexfemale"),
             set_prior("exponential(1)", class = "sigma"))

# Student-t likelihood prior includes nu
prior_B <- c(prior_A, set_prior("constant(8)", class = "nu"))


# Normal intercept prior (sensitivity)
prior_C <- c(set_prior("normal(0, 1)", class = "b"),
             set_prior("normal(0, 1)", class = "b", coef = "sexmale"),
             set_prior("normal(0, 1)", class = "b", coef = "sexfemale"),
             set_prior("exponential(1)", class = "sigma"))

prior <- c(set_prior("normal(0, 0.5)", class = "b"),
           set_prior("normal(0, 2)", class = "b", coef = "sexmale"),
           set_prior("normal(0, 2)", class = "b", coef = "sexfemale"),
           set_prior("exponential(1)", class = "sigma"))

message("Fitting Model A (gaussian likelihood + student-t intercept prior)...")
#fit_A <- brm(formula = f_beta,
#             data = data,
#             family = gaussian(),
#             prior = prior_C,
#             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

message("Fitting Model B (student-t likelihood + student-t prior)...")
#fit_B <- brm(formula = f_beta,
#             data = data,
#             family = student(),
#             prior = prior_B,
#             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

message("Fitting Model C (gaussian likelihood + normal intercept prior - sensitivity)...")
fit_C <- brm(formula = f_beta,
             data = data,
             family = gaussian(),
             prior = prior,
             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))


#loo_A <- loo(fit_A, moment_match = TRUE)
#loo_B <- loo(fit_B, moment_match = TRUE)
#loo_C <- loo(fit_C, moment_match = TRUE)
#print(loo_compare(loo_A, loo_B, loo_C))


#########2

f_beta1 <- bf(beta ~ 0+ sex + weight_s + height_s)

fit_example <- brm(formula = f_beta1,
             data = data,
             family = gaussian(),
             prior = prior,
             chains = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.98, max_treedepth = 15))

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
    sex = 'female',
    TBW = ifelse(
      sex == "Male",
      2.447 - 0.09516 * age + 0.1074 * height + 0.3362 * weight,   # Male formula
      -2.097 + 0.1069 * height + 0.2466 * weight),
    T_Vd = TBW / weight
  )

beta_draws <- posterior_predict(fit_example, newdata = new, draws = 12000)
Ct <- 0.15   
t <- 2      
C0_draws <- Ct + exp(beta_draws) * t
P_over <- mean(C0_draws > 0.47)        
quantile(C0_draws, prob = c(0.025,0.25, 0.5, 0.975))
P_over
# empirical 2.5th percentile of beta distribution:
emp_beta_2.5 <- quantile(data$beta, probs = 0.025, na.rm = TRUE)
quantile(C0_draws, probs = 0.025, na.rm = TRUE)




##########3


summary(lm(Vd~ 0 + T_Vd:sex, data))

f_beta <- bf(beta ~ 0+ sex + weight_s + height_s + drinkingtime_s)
f_Vd<- bf(Vd ~ 0 + T_Vd:sex)

ggplot(data, aes(x = Vd)) + geom_histogram(bins = 30) + ggtitle("Histogram of raw Vd (naive calc)")
ggplot(data, aes(x = Vd, y = beta)) + geom_point() + geom_smooth(method = "loess") +
  labs(title = "beta vs log(Vd)")



priors_joint <- c(
  set_prior("normal(0, 0.5)", class = "b", resp = "beta"),
  set_prior("normal(0, 2)", class = "b", coef = "sexmale", resp = "beta"),
  set_prior("normal(0, 2)", class = "b", coef = "sexfemale", resp = "beta"),
  set_prior("exponential(1)", class = "sigma", resp = "beta"),
  set_prior("normal(1/0.838, 0.2)", class = "b", resp = "Vd", coef = "T_Vd:sexfemale"),
  set_prior("normal(1/0.825, 0.2)", class = "b", resp = "Vd", coef = "T_Vd:sexmale"),
  set_prior("exponential(1)", class = "sigma", resp = "Vd")
)

# Fit the joint model (this may take time). Use student family if you prefer robust errors.
fit_joint <- brm(
  formula = f_beta + f_Vd + set_rescor(TRUE),
  data = data,
  prior = priors_joint,
  chains = 4, iter = 4000, warmup = 1500,
  control = list(adapt_delta = 0.98, max_treedepth = 15),
  seed = 2025,
  save_pars = save_pars(all = TRUE)
)
print(summary(fit_joint))

pp_check(fit_joint, resp = "beta", type = "dens_overlay", ndraws = 200)
pp_check(fit_joint, resp = "Vd", type = "dens_overlay", ndraws = 200)
loo_joint <- loo(fit_joint, moment_match = TRUE)
print(loo_joint)

epred_beta  <- posterior_epred(fit_joint, resp = "beta",  ndraws = 1000)  # 200 x N matrix
epred_Vd    <- posterior_epred(fit_joint, resp = "Vd",    ndraws = 1000)

pred_mean_beta <- colMeans(epred_beta)
pred_mean_Vd   <- colMeans(epred_Vd)
plot(pred_mean_Vd, data$Vd, xlab = "Predicted mean Vd", ylab = "Observed Vd")
abline(0,1, col = "red")

plot(exp(pred_mean_beta), exp(data$beta), xlab = "Predicted mean beta", ylab = "Observed beta")
abline(0,1, col = "red")

post <- as_draws_df(fit_joint)

#pp_new_beta  <- posterior_epred(fit_joint, newdata = new, resp = "beta", draws = 4000, re_formula = NULL)
#pp_new_Vd <- posterior_epred(fit_joint, newdata = new, resp = "Vd", draws = 4000, re_formula = NULL)
# Note: epred returns expected mean (no residual). If you want individual draws including residual noise:
#pp_new_beta_ind  <- posterior_predict(fit_joint, newdata = new, resp = "beta", draws = 4000)
#pp_new_Vd_ind <- posterior_predict(fit_joint, newdata = new, resp = "Vd", draws = 4000)
#beta_draws  <- as.vector(pp_new_beta_ind)
#Vd_draws <- as.vector(pp_new_Vd_ind)

pp_Vd <- posterior_predict(fit_joint, resp = "Vd", ndraws = 4000)
pp_Vd1 <- posterior_epred(fit_joint, resp = "Vd", ndraws = 4000)

Vd_obs <- data$Vd



Vd_low  <- apply(pp_Vd, 2, quantile, 0.025)
Vd_high <- apply(pp_Vd, 2, quantile, 0.975)
mean(data$Vd >= Vd_low & data$Vd <= Vd_high)


Vd_low  <- apply(pp_Vd1, 2, quantile, 0.25)
Vd_high <- apply(pp_Vd1, 2, quantile, 0.75)
mean(data$Vd >= Vd_low & data$Vd <= Vd_high)
data <- data %>% mutate(TV = ifelse(
                         sex == "Male",
                         1.40864*T_Vd,   # Male formula
                         1.29720*T_Vd))            


Vd_mean <- colMeans(pp_Vd)
plot(Vd_mean, data$Vd, pch=19)
abline(0,1,col="red")
plot(data$TV, data$Vd, pch=19)
abline(0,1,col="red")
PIT_Vd <- sapply(1:ncol(pp_Vd), function(i) mean(pp_Vd[,i] <= Vd_obs[i]))
hist(PIT_Vd, breaks=20)


pp_beta <- posterior_predict(fit_joint, resp = "beta", ndraws = 4000)
beta_obs <- data$beta
beta_low  <- apply(pp_beta, 2, quantile, 0.025)
beta_high <- apply(pp_beta, 2, quantile, 0.975)
mean(beta_obs >= beta_low & beta_obs <= beta_high)

beta_mean <- colMeans(pp_beta)
plot(beta_mean, beta_obs, pch=19)
abline(0,1,col="red")
PIT_beta <- sapply(1:ncol(pp_beta), function(i) mean(pp_beta[,i] <= beta_obs[i]))
hist(PIT_beta, breaks=20)

MSE_Vd <- mean(abs(Vd_mean - data$Vd) , na.rm = TRUE )
MSE_Vd
MSE_TV <- mean(abs(data$TV - data$Vd) , na.rm = TRUE )
MSE_TV

beta_A <- posterior_predict(fit_C, ndraws = 4000)

MSE_beta<- mean(abs(beta_mean - data$beta) , na.rm = TRUE )
MSE_beta
MSE_betaA <- mean(abs(beta_A - data$beta) , na.rm = TRUE )
MSE_betaA
#prior noninformative or theoratical for coefficie

###################

pp_draws <- exp(posterior_predict(fit_C, ndraws = 12000))

n_draws <- nrow(pp_draws)
n_obs <- ncol(pp_draws)         

# 2) For each observation compute 50% and 95% predictive intervals (from pp_draws)
alpha_lo_95 <- 0.025; alpha_hi_95 <- 0.975
alpha_lo_50 <- 0.25;  alpha_hi_50 <- 0.75

pp_summary <- tibble(
  obs = seq_len(n_obs),
  beta_obs = exp(data$beta)  # observed true beta for each row
) %>%
  mutate(
    pred_mean = colMeans(pp_draws),
    pred_median = apply(pp_draws, 2, median),
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
r2_draws <- as.numeric(bayes_R2(fit_C))  # ensures numeric vector
r2_median <- median(r2_draws)
r2_CI <- quantile(r2_draws, c(0.025, 0.975))



#Plot area
###################################################################################

F1 <- ggplot(data, aes(x = beta)) +
  geom_density() +
  labs(title = "Density Plot", x = "log(Beta)", y = "Density")

###################################################################################

fit_d <- brm(formula = f_beta,
             data = data,
             family = gaussian(),
             prior = prior_C,
             chains = 4, iter = 100, warmup = 30, control = list(adapt_delta = 0.98, max_treedepth = 15))

all <- as_draws_df(fit_d, inc_warmup = TRUE)

F2 <- ggplot(all, aes(x = .iteration,
                y = b_sexmale,
                color = factor(.chain))) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  labs(x = "Iteration (including warmup)",
       y = "b_sexfemale",
       color = "Chain") +
  theme_bw()

###################################################################################

priors_table <- tribble(~Parameter, ~Prior,
                        "Regression coefficient for male ", "Normal(0, 2)",
                        "Regression coefficient for female", "Normal(0, 2)",
                        "Regression coefficients (others)", "Normal(0, 0.5)",
                        "Residual SD ", "Exponential(1)")

T1 <- knitr::kable(priors_table, 
                   caption = "Weakly informative priors used in the Bayesian regression model")

###################################################################################

F3 <- mcmc_plot(fit_C, type = "dens_overlay")

###################################################################################

F4 <- mcmc_trace(fit_C, 
                 pars = c("b_sexfemale", "b_sexmale",
                          "b_weight_s", "b_height_s",
                          "b_drinkingtime_s", "sigma"))

###################################################################################

s <- summary(fit_C)
tab_all <- rbind(as.data.frame(s$fixed),
                 as.data.frame(s$random$ID))
T2 <- knitr::kable(
      tab_all,
      digits = 3,
      caption = "Posterior Summary from BRM")

###################################################################################

yrep_C <- exp(posterior_predict(fit_C, draws = 2000))
F5 <- ppc_dens_overlay(exp(data$beta), yrep_C[1:2000, ]) + ggtitle("PPC density")

###################################################################################

F6 <- pp_check(fit_C, type = "loo_pit_qq")

###################################################################################

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

df_C0 <- tibble(C0 = C0_draws)

F7 <- ggplot(df_C0, aes(x = C0)) +
  geom_density(fill = "grey", alpha = 0.4) +
  geom_vline(xintercept = 0.47, color = "black", linetype = "dashed") +
  labs(title = "Posterior Predictive Distribution of C_0",
       x = "C_0 (g/L)",
       y = "Density") +
  theme_minimal(base_size = 14)

###################################################################################

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
T4 <- knitr::kable(C0_summary,
                   digits = 3,
                   caption = "Example results")

###################################################################################


F1000 <- ggplot(data, aes(x = exp(data$beta), y = colMeans(pp_draws), colour = sex)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Observed vs Predicted", x = "Observed", y = "Predicted") + 
  theme(aspect.ratio = 1)

