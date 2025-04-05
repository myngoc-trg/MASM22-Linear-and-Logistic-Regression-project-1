#Project1####
library(tidyverse)
getwd()  # Check current working directory
list.files("Data")  # Check if 'carotene.xlsx' is inside the "Data" folder
install.packages("readxl")  # Only run once if not installed
library(readxl)
library(dplyr)  # Load dplyr for select()
# Lecture 3 - Elasticity - multiple linear regression

library(car)
library(rstatix)
library(GGally)
plasmaB <- read_excel("Data/carotene.xlsx")
plasmaB_bmi <- select(plasmaB, bmi, betaplasma)
plasmaB_bmi

#1 Plasma β-carotene and body mass index
##1(a)
###Linear model
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

#### Make a Q-Q-plot and a histogram for the residuals.
ggplot(data = plasmaB_bmi, aes(sample = e)) +
  geom_qq(size = 3) + geom_qq_line() +
  labs(tag = "C") +
  labs(title = "Normal Q-Q-plot of the residuals model 1")

ggplot(data = plasmaB_bmi, aes(x = e)) + 
  geom_bar() + scale_x_binned() +
  labs(title = "Histogram for the residuals model 1")

###Logarithm model
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

#### Make a Q-Q-plot and a histogram for the residuals.
ggplot(data = logplasmaB_bmi, aes(sample = e)) +
  geom_qq(size = 3) + geom_qq_line() +
  labs(tag = "C") +
  labs(title = "Normal Q-Q-plot of the residuals")

ggplot(data = logplasmaB_bmi, aes(x = e)) + 
  geom_bar() + scale_x_binned()


## 1(b)
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

###Transform the relationship back to betaplasma =e^Beta0*e^Beta1*bmi*e^epsiloni
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

## 1(c)
### (i) When BMI is increased by 1 unit 
plasmaB_bmiplus1_lm <- lm(betaplasma ~ I(bmi + 1), data = plasmaB_bmi)
plasmaB_bmiplus1_lm
confint(plasmaB_bmiplus1_lm)

### (ii) When BMI is decreased by 1 unit 
plasmaB_bmiminus1_lm <- lm(betaplasma ~ I(bmi - 1), data = plasmaB_bmi)
confint(plasmaB_bmiminus1_lm)
confint(plasmaB_bmi_lm)

### (iii) When BMI is decreased by 10 unit 
plasmaB_bmiminus10_lm <- lm(betaplasma ~ I(bmi - 10), data = plasmaB_bmi)
confint(plasmaB_bmiminus10_lm)

# Obs: CI for B1 is unchanged while the CI for intercept B0 changes due to the change in BMI

## 1(d) Lec.4
# Test whether there is a significant (α = 5 %) linear relationship between log plasma β-carotene BMI
# t-test
# Hypotheses. H0: B_1 = 0, H1: B_1 not = 0
glimpse(logplasmaB_bmi)
head(logplasmaB_bmi)
logplasmaB_bmi_sum <- summary(logplasmaB_bmi_lm)
logplasmaB_bmi_sum
# Distibution of B_1 when H0 true is t(315-1)
# p-value very small, can reject H0 with very high significance level

#2 Plasma β-carotene and smoking habit
## 2(a), turn the categorical variable smokstat into a factor variable
glimpse(plasmaB)
mutate(plasmaB,
       smokstat = factor(smokstat,
                     levels = c(1, 2, 3),
                     labels = c("never", "former", "current"))) -> plasmaB
plasmaB$soil <- NULL
glimpse(plasmaB)
plasmaB_lm <- lm(betaplasma ~ smokstat, data = plasmaB)
plasmaB_sum <- summary(plasmaB_lm)
plasmaB_sum
confint(plasmaB_lm)

plasmaB_x0 <- data.frame(smokstat = c("never", "former", "current"))
cbind(plasmaB_x0,
      predict(plasmaB_lm, plasmaB_x0, se.fit = TRUE),
      conf = predict(plasmaB_lm, plasmaB_x0, interval = "confidence"),
      pred = predict(plasmaB_lm, plasmaB_x0, interval = "prediction")) |>
  mutate(df = NULL, residual.scale = NULL,
         conf.fit = NULL, pred.fit = NULL,
         se.pred = sqrt(plasmaB_sum$sigma^2 + se.fit^2))
logplasmaB_lm <- lm(log(betaplasma) ~ smokstat, data = plasmaB)
logplasmaB_lm
logplasmaB_sum <- summary(logplasmaB_lm)
logplasmaB_sum

cbind(plasmaB_x0,
      predict(logplasmaB_lm, plasmaB_x0, se.fit = TRUE),
      conf = predict(logplasmaB_lm, plasmaB_x0, interval = "confidence"),
      pred = predict(logplasmaB_lm, plasmaB_x0, interval = "prediction")) |>
  mutate(df = NULL, residual.scale = NULL,
         conf.fit = NULL, pred.fit = NULL,
         se.pred = sqrt(logplasmaB_sum$sigma^2 + se.fit^2))
# ANSWER: which category would be suited to use as reference category
# Never??? Maybe 
?boxplot
boxplot(plasmaB)

# Create boxplots & violin plots for plasma β-carotene
p1 <- ggplot(plasmaB, aes(x = smokstat, y = betaplasma, fill = smokstat)) +
  geom_violin(alpha = 0.5) +  # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Boxplot overlay
  theme_minimal() +
  labs(title = "Plasma β-Carotene by Smoking Status",
       x = "Smoking Status",
       y = "Plasma β-Carotene (mg/kg)") +
  scale_fill_brewer(palette = "Pastel1")
p1

p2 <- ggplot(plasmaB, aes(x = smokstat, y = log(betaplasma), fill = smokstat)) +
  geom_violin(alpha = 0.5) +  # Violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA) +  # Boxplot overlay
  theme_minimal() +
  labs(title = "Log Plasma β-Carotene by Smoking Status",
       x = "Smoking Status",
       y = "Log Plasma β-Carotene (mg/kg)") +
  scale_fill_brewer(palette = "Pastel1")
p2
# p1 is highly skewed, p2 is more symmetric
# So we should still use log plasma B as dependent variable

## 2(b)
# "never" as reference
plasmaB <- mutate(plasmaB, smokstat = relevel(smokstat, "never"))
# B0 is "never", B1 "former", B2 "current"
logplasmaB_smokstat_lm <- lm(log(betaplasma) ~ smokstat, data = plasmaB)
logplasmaB_smokstat_sum <- summary(logplasmaB_smokstat_lm)
logplasmaB_smokstat_sum

# "current" as reference
# B0 is "current", B1 is "never", B2 is "never"
plasmaB_current <- mutate(plasmaB, smokstat = relevel(smokstat, "current"))
# B0 is "never", B1 "former", B2 "current"
logplasmaB_current_lm <- lm(log(betaplasma) ~ smokstat, data = plasmaB_current)
logplasmaB_current_sum <- summary(logplasmaB_current_lm)
logplasmaB_current_sum
# Model with "never" reference has lower std. errors
# EXPLAIN?


## 2(c)
cbind(plasmaB_x0,
      predict(logplasmaB_lm, plasmaB_x0, se.fit = TRUE),
      conf = predict(logplasmaB_lm, plasmaB_x0, interval = "confidence"),
      pred = predict(logplasmaB_lm, plasmaB_x0, interval = "prediction")) |>
  mutate(df = NULL, residual.scale = NULL,
         conf.fit = NULL, pred.fit = NULL,
         se.pred = sqrt(logplasmaB_sum$sigma^2 + se.fit^2))

cbind(plasmaB_x0,
      predict(logplasmaB_current_lm, plasmaB_x0, se.fit = TRUE),
      conf = predict(logplasmaB_current_lm, plasmaB_x0, interval = "confidence"),
      pred = predict(logplasmaB_current_lm, plasmaB_x0, interval = "prediction")) |>
  mutate(df = NULL, residual.scale = NULL,
         conf.fit = NULL, pred.fit = NULL,
         se.pred = sqrt(logplasmaB_current_sum$sigma^2 + se.fit^2))

# Relate these expected values to the corresponding means in 2(a
# Explain why the predictions and their confidence intervals are the same regardless of which model version you used

## 2(d)
# Perform a suitable test for whether there are significant (α = 5 %) differences in log plasma β-caroteneLec between any of the smokstat categories
logplasmaB_sum
# Global F-test
# Hypotheses. H0: B_1 = B_2 = 0, H1: B_j not 0
# F_obs = 5.75 = MS(Regr) / MS(Error)
# F(2, 312)
# p < 0.05, reject H0
# There is significant difference in logplasmaB between any of the smokstat categories 

#3 Multiple linear regression
## 3(a)
# Turn sex and vituse into factor variables
# 
glimpse(plasmaB)
mutate(plasmaB,
       sex = factor(sex,
                         levels = c(1, 2),
                         labels = c("male", "female"))) -> plasmaB

mutate(plasmaB,
       vituse = factor(vituse,
                         levels = c(1, 2, 3),
                         labels = c("often", "sometimes", "no"))) -> plasmaB
plasmaB_x1 <- data.frame(sex = c("male", "female"))


plasmaB_sex_lm <- lm(betaplasma ~ sex, data = plasmaB)
plasmaB_sex_sum <- summary(plasmaB_sex_lm)
cbind(plasmaB_x1,
      predict(plasmaB_sex_lm, plasmaB_x1, se.fit = TRUE),
      conf = predict(plasmaB_sex_lm, plasmaB_x1, interval = "confidence"),
      pred = predict(plasmaB_sex_lm, plasmaB_x1, interval = "prediction")) |>
  mutate(df = NULL, residual.scale = NULL,
         conf.fit = NULL, pred.fit = NULL,
         se.pred = sqrt(plasmaB_sex_sum$sigma^2 + se.fit^2))

plasmaB_x2 <- data.frame(vituse = c("often", "sometimes", "no"))
plasmaB_vituse_lm <- lm(betaplasma ~ vituse, data = plasmaB)
plasmaB_vituse_sum <- summary(plasmaB_vituse_lm)
cbind(plasmaB_x2,
      predict(plasmaB_vituse_lm, plasmaB_x2, se.fit = TRUE),
      conf = predict(plasmaB_vituse_lm, plasmaB_x2, interval = "confidence"),
      pred = predict(plasmaB_vituse_lm, plasmaB_x2, interval = "prediction")) |>
  mutate(df = NULL, residual.scale = NULL,
         conf.fit = NULL, pred.fit = NULL,
         se.pred = sqrt(plasmaB_vituse_sum$sigma^2 + se.fit^2))

table(plasmaB$sex)
table(plasmaB$vituse)

# "female" as sex will be reference as bigger category
# "often" as vitus: reference

# 3(b)
# Calculate pairwise correlations between the continuous x-variables bmi, age, calories, fat, cholesterol, fiber, alcohol, and betadiet
plasmaB |> select(bmi, age, calories, fat, cholesterol, fiber, alcohol, betadiet) |> 
  cor_test() |> filter(var1 < var2)

plasmaB |> select(bmi, age, calories, fat, cholesterol, fiber, alcohol, betadiet) |> 
  cor_test() |> filter(var1 < var2, abs(cor) > 0.6) 

# Compute correlation matrix and filter strong correlations
highcor <- plasmaB |> 
  select(bmi, age, calories, fat, cholesterol, fiber, alcohol, betadiet) |> 
  cor_test() |> 
  as.data.frame()  
  filter(var1 < var2, abs(cor) > 0.6)  # Keep only strong correlations
highcor
# Extract variable names with strong correlation
highcor <- unique(c(highcor$var1, highcor$var2))

# Plot pairwise relationships for selected variables
plasmaB |> 
  select(all_of(highcor)) |> 
  ggpairs(lower = list(continuous = wrap("points", size = 0.5)))

# Is the person consuming the equivalent 200 alcoholic drinks, corresponding to 12 bottles of Absolut Vodka, each week, extreme in any other variables as well?
# Alcohol and calories should be correlated
# Plot alcohol vs calories
ggplot(data = plasmaB, mapping = aes(x = alcohol, y = calories)) +
  geom_point(size = 1)
# There is some significant relationship between alcohol vs calories
# It seems like alcohol consumption increases calrories consumption very largely


# 3(c)
# Ignore any potential problems 
# a model where log plasma β-carotene depends on all the other variables, bmi, age, calories, fat, cholesterol, fiber, alcohol, betadiet, smokstat, sex, and vituse
head(plasmaB)
glimpse(plasmaB)
?lm
logplasmaB_all_lm <- lm(log(betaplasma) ~ bmi + age + calories + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, data = plasmaB)
logplasmaB_all_sum <- summary(logplasmaB_all_lm)
logplasmaB_all_sum
# Present the VIF/GVIF-values for the variables
# dicate any variables where more than 80 % of the variablility can be explained using the other x-variables
vif(logplasmaB_all_lm)
# calories: GVIF^(1/(2*Df)) = 3.63 > 3.16
# calories has large dependency on other variables

# Remove the most problematic x-variable and refit the model without it 
# (Model.3(c)). 
# Present, and comment on, the new VIF/GVIF-values.
logplasmaB_wocalories_lm <- lm(log(betaplasma) ~ bmi + age + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, data = plasmaB)
logplasmaB_wocalories_sum <- summary(logplasmaB_wocalories_lm)
logplasmaB_wocalories_sum
vif(logplasmaB_wocalories_lm)
# The model looks good with GVIF-values


# 3(d)
logplasmaB_wocalories_sum
# x1 = bmi, x2 = age, x3 = fat, x4 = cholesterol, x5 = fiber, x6 = alcohol, x7 = betadiet, 
# Dummy variales: x8 = smokstatformer, x9 = smokstatcurrent, x10 = smokstatcurrent
# Dummy variables: x11 = sexfemale
# Dummy variables: x12 = vitusesometimes, x12 = vituseno
# Model: Y_i = x_i.T dot B + epsilon_i

# 3(d) i
# Is there a significant relationship between log plasma β-carotene and BMI, given the other variables in the model?
# p-value for B_bmi  is very small, can reject H0 with a very high significance level
# i.e. there is a significant relationship between log plasma and bmi

# 3(d) ii
# Is this model significantly better than Model 1(b),which used only bmi?
# Partial F-test
logplasmaB_bmi_lm
anova(logplasmaB_bmi_lm, logplasmaB_wocalories_lm)
# p-value very small, reject H0
# This new model is significantly better

# 3(d) iii
# Is this model significantly better than Model 2(b),which used only smokstat?
# Partial F-test
anova(logplasmaB_smokstat_lm, logplasmaB_wocalories_lm)
# p-value very small, reject H0
# This new model is significantly better


# 3(e)
# Make a visual inspection of the studentized residuals for Model 3(C)
# looking for outliers, non-constant variance, and non-normality, using suitable plots with suitable reference line
plasmaB <- mutate(plasmaB, sex = relevel(sex, "female"))
plasmaB_wocalories_lm <- lm(betaplasma ~ bmi + age + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, data = plasmaB)
plasmaB_wocalories_sum <- summary(plasmaB_wocalories_lm)
plasmaB_wocalories_sum
logplasmaB_wocalories_sum
plasmaB_pred <- mutate(plasmaB, 
                       yhat_linear = predict(plasmaB_wocalories_lm),
                       r_linear = rstudent(plasmaB_wocalories_lm),
                       yhat_log = predict(logplasmaB_wocalories_lm),
                   r_log = rstudent(logplasmaB_wocalories_lm),
                   v_log = hatvalues(logplasmaB_wocalories_lm),
                   D_log = cooks.distance(logplasmaB_wocalories_lm))
glimpse(plasmaB_pred)
logplasmaB_wocalories_sum
highlightcolors <- c("|r*|>3" = "red")

ggplot(plasmaB_pred, aes(x = yhat_log, y = r_log)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_point(data = filter(plasmaB_pred, abs(r_log) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Studentized residuals vs log predictor",
       subtitle = "Log-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")
# Looks good-ish


# 3(f)
# Calculate the leverage for Model 3(c)
# plot them against the linear predictor, with suitable reference lines.
# with 1/n and 2(p+1)/n horizontal lines:
# p+1 = 
pplus1 <- length(logplasmaB_wocalories_lm$coefficients)
n <- nobs(logplasmaB_wocalories_lm)
glimpse(plasmaB_pred)
highlightcolors <- c("|r*|>3" = "red",
                     "largest leverage" = "magenta", 
                     "all data" = "orange")
ggplot(cbind(plasmaB_pred), aes(x = yhat_linear, y = v_log)) +
  geom_point(size = 2) +
  #geom_point(data = filter(plasmaB_pred, calories < 30), 
             #aes(color = "calories<30"), size = 3) +
  geom_hline(yintercept = 1/n) +
  geom_hline(yintercept = 2*pplus1/n, color = "red") +
  labs(title = "Leverage for log-lin model vs predicted plasmaB",
       caption = "y = 1/n (black) and 2(p+1)/n (red)",
       color = "Highlight") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)

which_max_leverage <- which.max(plasmaB_pred$v_log)
which_max_leverage
which_max_leverage <- as.data.frame(plasmaB_pred[which_max_leverage, ])
glimpse(plasmaB_pred)
# The one has max leverage has extremely low calories intake (27.9), extremely high cholesterol, 
# and relatively high fat 64.3 > 3rd qu., high alcohol intake, betadiet
# This is potentially influential obs, far from the centre of gravity of X-space
# May not be actually influential
# Let's do more test, exclusion?
# Using suitable plots? what?


ggplot(plasmaB_pred, aes(calories, betaplasma)) + geom_point() +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  labs(title = "Betaplasma by calories",
       sutitle = "including the strange one",
       color = "Highlight") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)

ggplot(plasmaB_pred, aes(alcohol, betaplasma)) + geom_point() +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  labs(title = "Betaplasma by alcohol intake",
       sutitle = "including the strange one",
       color = "Highlight") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)


# Maybe not necessary to include this residual plot to show where the largest leverage one locates
# Not sure
ggplot(plasmaB_pred, aes(x = yhat_log, y = r_log)) +
  geom_point() +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  geom_hline(yintercept = c(-2, 0, 2)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_point(data = filter(plasmaB_pred, abs(r_log) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Studentized residuals vs log predictor",
       subtitle = "Log-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")
# residual for that weird person is fine tho?
# what plots?

ggplot(plasmaB_pred, aes(x = yhat_log, y = sqrt(abs(r_log)))) +
  geom_point() +
  geom_hline(yintercept = c(0, sqrt(qnorm(0.75)), sqrt(2))) +
  geom_hline(yintercept = sqrt(3), linetype = 2) +
  geom_point(data = filter(cod_pred, abs(r_linear) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Sqrt absolute studentized residuals vs l predictor",
       subtitle = "Lin-lin model",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")


# 3(g)
# Calculate Cook’s distance for Model.3(c) and plot against the linear predictor
f1.plasmaB <- pplus1
f2.plasmaB <- logplasmaB_wocalories_lm$df.residual
cook.limit.plasmaB <- qf(0.5, f1.plasmaB, f2.plasmaB)
glimpse(plasmaB_pred)
ggplot(plasmaB_pred, aes(yhat_linear, D_log)) + 
  geom_point(size = 1) +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  geom_hline(yintercept = cook.limit.plasmaB, color = "red") +
  geom_hline(yintercept = 4/n, linetype = 2, color = "red") +
  xlab("Fitted values for linear predictor") +
  ylab("D_i") +
  labs(title = "Pike: Cook's D",
       caption = "4/n (dashed), F_0.5, p+1, n-(p+1) (solid)",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors)
# The one with largest leverage doesn't seem to have a large influence

which_max_cook <- which.max(plasmaB_pred$D_log)
which_max_cook
which_max_cook <- as.data.frame(plasmaB_pred[which_max_cook,])
summary(plasmaB_pred)

# Calculate the DFBETAS
glimpse(dfbetas(logplasmaB_wocalories_lm))
# Identify the β-parameters that have been affected
glimpse(logplasmaB_wocalories_sum)
plasmaB_pred <- mutate(
  plasmaB_pred,
  df0 = dfbetas(logplasmaB_wocalories_lm)[, "(Intercept)"],
  df1 = dfbetas(logplasmaB_wocalories_lm)[, "bmi"],
  df2 = dfbetas(logplasmaB_wocalories_lm)[, "age"],
  df3 = dfbetas(logplasmaB_wocalories_lm)[, "fat"],
  df4 = dfbetas(logplasmaB_wocalories_lm)[, "cholesterol"],
  df5 = dfbetas(logplasmaB_wocalories_lm)[, "fiber"],
  df6 = dfbetas(logplasmaB_wocalories_lm)[, "alcohol"],
  df7 = dfbetas(logplasmaB_wocalories_lm)[, "betadiet"],
  df8 = dfbetas(logplasmaB_wocalories_lm)[, "smokstatformer"],
  df9 = dfbetas(logplasmaB_wocalories_lm)[, "smokstatcurrent"],
  df10 = dfbetas(logplasmaB_wocalories_lm)[, "sexmale"],
  df11 = dfbetas(logplasmaB_wocalories_lm)[, "vitusesometimes"],
  df12 = dfbetas(logplasmaB_wocalories_lm)[, "vituseno"])
highlightshapes <- c("Cook's D>0.1" = 24)
ggplot(plasmaB_pred, aes(x = yhat_linear, y = df4)) +
  geom_point(size = 2) +
  geom_point(data = filter(plasmaB_pred, abs(r_log) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(plasmaB_pred, D_log > 0.1),
             aes(shape = "Cook's D>0.1"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.plasmaB)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_0(i)") +
  xlab("Fitted values for log predictor") +
  labs(title = "Pike: DFBETAS_0: impact on the intercept",
       #subtitle = "without the strange fish",
       caption = "y = sqrt(F_0.5) and 2/sqrt(n)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

#df4,df1, df10, df11, 
# df4, cholesterol seems to bed affected by the yhat_log = 4.3
which_max_df4 <- which.max(plasmaB_pred$df4)
which_max_df4
which_max_df4 <- as.data.frame(plasmaB_pred[which_max_df4,])

# Ilustrate by plotting the relevant DFBETAS against the corresponding x-variable, with suitable reference lines
ggplot(plasmaB_pred, aes(x = cholesterol, y = df4)) +
  geom_point(size = 2) +
  geom_point(data = filter(plasmaB_pred, abs(r_log) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(plasmaB_pred, D_log > 0.1),
             aes(shape = "Cook's D>0.1"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.plasmaB)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_0(i)") +
  xlab("x-values for cholesterol") +
  labs(title = "Pike: DFBETAS_0: impact on the intercept",
       #subtitle = "without the strange fish",
       caption = "y = sqrt(F_0.5) and 2/sqrt(n)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

# High cholesterol value decreases the Beta_4 parameters
# Explain why the observation has had a large influence on the estimate, with the help of suitable plots