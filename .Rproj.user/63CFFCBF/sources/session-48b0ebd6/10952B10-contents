#Project1####
library(tidyverse)
getwd()  # Check current working directory
list.files("Data")  # Check if 'carotene.xlsx' is inside the "Data" folder
install.packages("readxl")  # Only run once if not installed
library(readxl)
library(dplyr)  # Load dplyr for select()
plasmaB <- read_excel("Data/carotene.xlsx")
plasmaB_bmi <- select(plasmaB, bmi, betaplasma)
plasmaB_bmi

#1 Plasma β-carotene and body mass index
#1(a)
#Linear model
summary(plasmaB)

head(plasmaB)

plasmaB_bmi_lm <- lm(betaplasma ~ bmi, data = plasmaB_bmi)
plasmaB_bmi_lm
plasmaB_bmi <- mutate(plasmaB_bmi, yhat = predict(plasmaB_bmi_lm))
head(plasmaB_bmi)
plasmaB_bmi <- mutate(plasmaB_bmi, e = plasmaB_bmi_lm$residuals)
glimpse(plasmaB_bmi)
plasmaB_elims <- max(abs(plasmaB_bmi$e))*c(-1, 1)
plasmaB_elims
ggplot(data = plasmaB_bmi, aes(x = yhat, y = e)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0) +
  geom_smooth() +
  expand_limits(y = plasmaB_elims) +
  xlab("Predicted plasma β-carotene (ng/ml)") +
  ylab("Residual") +
  labs(tag = "B") +
  labs(title = "Residuals vs predicted values Y-hat")

# Make a Q-Q-plot and a histogram for the residuals.
ggplot(data = plasmaB_bmi, aes(sample = e)) +
  geom_qq(size = 3) + geom_qq_line() +
  labs(tag = "C") +
  labs(title = "Normal Q-Q-plot of the residuals")

ggplot(data = plasmaB_bmi, aes(x = e)) + 
  geom_bar() + scale_x_binned()

#Logarithm model
ggplot(plasmaB_bmi, aes(x = bmi, y = log(betaplasma))) + geom_point()
logplasmaB_bmi_lm <- lm(log(betaplasma) ~ bmi, data = plasmaB_bmi)
logplasmaB_bmi_lm
logplasmaB_bmi <- plasmaB_bmi %>% mutate(logbetaplasma = log(plasmaB_bmi$betaplasma))
logplasmaB_bmi$yhat <- NULL
logplasmaB_bmi$e <- NULL
logplasmaB_bmi
logplasmaB_bmi <- mutate(logplasmaB_bmi, yhat = predict(logplasmaB_bmi_lm))
head(logplasmaB_bmi)
logplasmaB_bmi <- mutate(logplasmaB_bmi, e = logplasmaB_bmi_lm$residuals)
glimpse(logplasmaB_bmi)
logplasmaB_elims <- max(abs(logplasmaB_bmi$e))*c(-1, 1)
logplasmaB_elims
ggplot(data = logplasmaB_bmi, aes(x = yhat, y = e)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0) +
  geom_smooth() +
  expand_limits(y = logplasmaB_elims) +
  xlab("Logarithm predicted plasma β-carotene (ng/ml)") +
  ylab("Residual") +
  labs(tag = "B") +
  labs(title = "Residuals vs log predicted values Y-hat")

# Make a Q-Q-plot and a histogram for the residuals.
ggplot(data = logplasmaB_bmi, aes(sample = e)) +
  geom_qq(size = 3) + geom_qq_line() +
  labs(tag = "C") +
  labs(title = "Normal Q-Q-plot of the residuals")

ggplot(data = logplasmaB_bmi, aes(x = e)) + 
  geom_bar() + scale_x_binned()


# 1(b)
CI_coeffs <- confint(logplasmaB_bmi_lm)
logplasmaB_bmi_lm
CI_coeffs <- mutate(logplasmaB_bmi_lm$Coeffcients, CI_coeffs)
#Plot the data together with the fitted line, 
#a 95 % confidence interval for the fitted line and 
#a 95 % prediction interval for new observations.
c(min(plasmaB_bmi$bmi), max(plasmaB$bmi))
bmi_seq <- data.frame(bmi = seq(16.5, 50))
bmi_seq
bmi_seq |> mutate(
  fit = predict(logplasmaB_bmi_lm, newdata = bmi_seq),
  conf = predict(logplasmaB_bmi_lm, newdata = bmi_seq, interval = "confidence"),
  pred = predict(logplasmaB_bmi_lm, newdata = bmi_seq, interval = "prediction")) -> 
  logplasmaB_bmi_ints
glimpse(logplasmaB_bmi_ints)
head(logplasmaB_bmi_ints)
logplasmaB_bmi_ints$pred[, "fit"]

ggplot(logplasmaB_bmi_ints, aes(x = bmi)) + 
  geom_point(data = logplasmaB_bmi, aes(y = logbetaplasma), size = 3) +
  geom_line(aes(y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = conf[, "lwr"], ymax = conf[, "upr"]), alpha = 0.2) +
  geom_line(aes(y = pred[, "lwr"]), color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = pred[, "upr"]), color = "red", linetype = "dashed", linewidth = 1) +
  xlab("bmi") +
  ylab("Log Plasma β-carotene (ng/ml)") +
  labs(title = "Log Plasma β-carotene level varies with BMI",
       caption = "data, fitted line, 95% confidence and prediction intervals")

#Transform the relationship back to betaplasma =e^Beta0*e^Beta1*bmi*e^epsiloni
logplasmaB_bmi_ints
backlogplasmaB_bmi_ints <- logplasmaB_bmi_ints %>% 
  mutate(across(-bmi, exp))
backlogplasmaB_bmi_ints
ggplot(backlogplasmaB_bmi_ints, aes(x = bmi)) + 
  geom_point(data = plasmaB_bmi, aes(y = betaplasma), size = 3) +
  geom_line(aes(y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = conf[, "lwr"], ymax = conf[, "upr"]), alpha = 0.2) +
  geom_line(aes(y = pred[, "lwr"]), color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = pred[, "upr"]), color = "red", linetype = "dashed", linewidth = 1) +
  xlab("bmi") +
  ylab("Plasma β-carotene (ng/ml)") +
  labs(title = "Plasma β-carotene level varies with BMI",
       caption = "data, fitted line, 95% confidence and prediction intervals")

# 1(c)
# (i) When BMI is increased by 1 unit 
plasmaB_bmiplus1_lm <- lm(betaplasma ~ I(bmi + 1), data = plasmaB_bmi)
plasmaB_bmiplus1_lm
confint(plasmaB_bmiplus1_lm)

# (ii) When BMI is decreased by 1 unit 
plasmaB_bmiminus1_lm <- lm(betaplasma ~ I(bmi - 1), data = plasmaB_bmi)
confint(plasmaB_bmiminus1_lm)
confint(plasmaB_bmi_lm)

# (iii) When BMI is decreased by 10 unit 
plasmaB_bmiminus10_lm <- lm(betaplasma ~ I(bmi - 10), data = plasmaB_bmi)
confint(plasmaB_bmiminus10_lm)

# Obs: CI for B1 is unchanged while the CI for intercept B0 changes due to the change in BMI
