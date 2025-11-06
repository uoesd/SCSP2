library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
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
    BACpeaktime = `BAC Peak time (min)`
  )


data$Beta60 <- abs(data$beta)  # because Beta60 in file is negative of slope. :contentReference[oaicite:2]{index=2}
library(brms)
# scale continuous covariates for stability
data$Age_z    <- scale(data$age)
data$Weight_z <- scale(data$weight)
data$Dose_z   <- scale(data$AAC)

# formula: change covariates as appropriate
f <- bf(Beta60 ~ sex + Age_z + Weight_z)

priors <- c(
  prior(normal(0,10), class = "b"),
  prior(cauchy(0,2), class = "sigma", lb = 0),
  prior(constant(7), class = "nu"))

fit <- brm(f, data = data, 
           prior = priors,
           family = student(),
           chains = 4, iter = 4000, warmup = 1000, seed = 42,
           control = list(adapt_delta = 0.98))


new <- data.frame(sex= 'female', Age_z=(70-mean(data$age))/sd(data$age),
                  Weight_z=(70-mean(data$weight))/sd(data$weight))

beta_draws <- posterior_predict(fit, newdata = new, draws = 2000)

Ct <- 0.15   
t <- 2      
C0_draws <- Ct + beta_draws * t

P_over <- mean(C0_draws > 0.47)        
quantile(C0_draws, prob = c(0.025,0.25, 0.5, 0.975))
P_over
