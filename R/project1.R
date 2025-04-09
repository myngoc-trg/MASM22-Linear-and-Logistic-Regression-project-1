#Project1####
library(tidyverse); theme_set(theme_bw() + theme(text = element_text(size = 18)))
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
                         labels = c("often", "notoften", "no"))) -> plasmaB
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
# Make a Q-Q-plot and a histogram for the residuals.
ggplot(data = plasmaB_pred, aes(sample = r_log)) +
  geom_qq(size = 2) + geom_qq_line() +
  labs(tag = "C") +
  labs(title = "Normal Q-Q-plot of the residuals")

ggplot(data = plasmaB_pred, aes(x = r_log)) + 
  geom_bar() + scale_x_binned()

# 3(f)
# Calculate the leverage for Model 3(c)
# plot them against the linear predictor, with suitable reference lines.
# with 1/n and 2(p+1)/n horizontal lines:
# p+1 = 
pplus1 <- length(logplasmaB_wocalories_lm$coefficients)
n <- nobs(logplasmaB_wocalories_lm)
glimpse(plasmaB_pred)
highlightcolors <- c("|r*|>3" = "red",
                     "largest leverage" = "purple", 
                     "all data" = "orange",
                     "largest cook's distance" = "pink")
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
summary(plasmaB_pred)
# The one has max leverage has extremely low calories intake (27.9), extremely high cholesterol, 
# and relatively high fat 64.3 > 3rd qu., high alcohol intake, betadiet
# This is potentially influential obs, far from the centre of gravity of X-space
# May not be actually influential
# Let's do more test, exclusion?
# Using suitable plots? what?

glimpse(plasmaB_pred)
ggplot(plasmaB_pred, aes(calories, betaplasma)) + geom_point() +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  labs(title = "Betaplasma by calories",
       sutitle = "including the strange one",
       color = "Highlight") +
  xlab("calories (MJ)") +
  ylab("Betaplasma (mg/ml)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)

ggplot(plasmaB_pred, aes(alcohol, betaplasma)) + geom_point() +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  labs(title = "Betaplasma by alcohol intake",
       #sutitle = "including the strange one",
       color = "Highlight") +
  xlab("alcohol, number of alcoholic drinks per week") +
  ylab("Betaplasma (mg/ml)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)

ggplot(plasmaB_pred, aes(cholesterol, betaplasma)) + geom_point() +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  labs(title = "Betaplasma by alcohol intake",
       sutitle = "including the strange one",
       color = "Highlight") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)

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
# df4 seems to be largest, beta_4 seems to be affected the most without obs i

# Calculate the DFBETAS
glimpse(dfbetas(logplasmaB_wocalories_lm))
# Identify the β-parameters that have been affected
glimpse(logplasmaB_wocalories_sum)
plasmaB_pred <- mutate(
  plasmaB_pred,
  df0 = dfbetas(logplasmaB_wocalories_lm)[, "(Intercept)"],
  df1 = dfbetas(logplasmaB_wocalories_lm)[, "age"],
  df2 = dfbetas(logplasmaB_wocalories_lm)[, "bmi"],
  df3 = dfbetas(logplasmaB_wocalories_lm)[, "fat"],
  df4 = dfbetas(logplasmaB_wocalories_lm)[, "fiber"],
  df5 = dfbetas(logplasmaB_wocalories_lm)[, "alcohol"],
  df6 = dfbetas(logplasmaB_wocalories_lm)[, "cholesterol"],
  df7 = dfbetas(logplasmaB_wocalories_lm)[, "betadiet"],
  df8 = dfbetas(logplasmaB_wocalories_lm)[, "smokstatformer"],
  df9 = dfbetas(logplasmaB_wocalories_lm)[, "smokstatcurrent"],
  df10 = dfbetas(logplasmaB_wocalories_lm)[, "sexmale"],
  df11 = dfbetas(logplasmaB_wocalories_lm)[, "vitusesometimes"],
  df12 = dfbetas(logplasmaB_wocalories_lm)[, "vituseno"])
highlightshapes <- c("Cook's D>0.1" = 24)
ggplot(plasmaB_pred, aes(x = yhat_linear, y = df6)) +
  geom_point(size = 2) +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  geom_point(data = which_max_cook, 
             aes(color = "largest cook's distance"), size = 3) +
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
  ylab("DFBETAS_6(i)") +
  xlab("Fitted values for log predictor") +
  labs(title = "DFBETAS_6: impact on the beta_6",
       #subtitle = "without the strange fish",
       caption = "y = sqrt(F_0.5) and 2/sqrt(n)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

#df4,df1, df10, df11, 
# df4, cholesterol seems to bed affected by the yhat_log = 4.3
which_max_df6 <- which.max(abs(plasmaB_pred$df6))
which_max_df6
which_max_df6 <- as.data.frame(plasmaB_pred[which_max_df6,])

# Ilustrate by plotting the relevant DFBETAS against the corresponding x-variable, with suitable reference lines
ggplot(plasmaB_pred, aes(x = bmi, y = df6)) +
  geom_point(size = 2) +
  geom_point(data = which_max_leverage, 
             aes(color = "largest leverage"), size = 3) +
  geom_point(data = which_max_cook, 
             aes(color = "largest cook's distance"), size = 3) +
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
  ylab("DFBETAS_6(i)") +
  xlab("cholesterol (mg)") +
  labs(title = "DFBETAS_6: impact on the cholesterol",
       #subtitle = "without the strange fish",
       caption = "y = sqrt(F_0.5) and 2/sqrt(n)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

glimpse(plasmaB_pred)




# High cholesterol value decreases the Beta_4 parameters
# Explain why the observation has had a large influence on the estimate, with the help of suitable plots
as.data.frame(filter(plasmaB_pred, abs(r_log) > 3, cholesterol > 900))
summary(plasmaB_pred)
# Outlier: extremely low calories (=12.6) even though bmi = 31.2 (overweight), highest cholesterol, lowest betaplasma
# High fat


# Investigate the observation with the largest leverage, identified in 3(f )
# determine whether it has had any alarming influence on the estimates
which_max_leverage
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
# Not really influential


# 4 Removing the influential observation
# 4(a)
# Create a new dataset where the influential observation identified in 3(g) has been removed 
# estimate a new version of Model.3(c) on this reduced data set

logplasmaB_wocalories_excl_lm <- update(logplasmaB_wocalories_lm, 
                                        data = filter(plasmaB, cholesterol < 900.6))
plasmaB_pred <- mutate(plasmaB_pred,
                    yhat_excl = predict(logplasmaB_wocalories_excl_lm, 
                                        newdata = plasmaB))
# Calculate and plot Cook’s distance against the linear predictor, with suitable reference lines, for Model.3(c) on the reduced data set
# compare with the corresponding plot for the full data set in 3.(g)
plasmaB_excl_pred <- cbind(
  filter(plasmaB, cholesterol < 900.6),
  fit = predict(logplasmaB_wocalories_excl_lm),
  r = rstudent(logplasmaB_wocalories_excl_lm),
  D = cooks.distance(logplasmaB_wocalories_excl_lm))
f1.excl <- length(logplasmaB_wocalories_excl_lm$coefficients)
f2.excl <- logplasmaB_wocalories_excl_lm$df.residual
cook.limit.excl <- qf(0.5, f1.excl, f2.excl)

n_excl <- nobs(logplasmaB_wocalories_excl_lm)
filter(plasmaB_excl_pred, D > 0.1)
glimpse(plasmaB_excl_pred)

highlightshapes <- c("Cook's D>0.1" = 24)

##Plot D vs yhat####
ggplot(plasmaB_excl_pred, aes(fit, D)) + 
  geom_point(size = 2) +
  geom_point(data = filter(plasmaB_excl_pred, abs(r) > 3),
             aes(color = "|r*|>3"), size = 2) +
  geom_point(data = filter(plasmaB_excl_pred, D > 0.1),
             aes(shape = "Cook's D>0.1"), size = 1) +
  #geom_hline(yintercept = cook.limit.excl, color = "red") +
  geom_hline(yintercept = 4/n_excl, 
             linetype = 2, color = "red") +
  xlab("Fitted values") +
  ylab("D_i") +
  labs(title = "Cook's D",
       subtitle = "without the strange person",
       caption = "4/n (dashed), F_0.5, p+1, n-(p+1) (solid)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors) +
  scale_shape_manual(values = highlightshapes)

vif(logplasmaB_wocalories_excl_lm)

# 4(b)
# For both data sets, full and reduced, perform a stepwise variable selection
# arting with the null model, using the null model as the smallest model allowed and the full model as the largest model allowed,
# with BIC as criterion

# For reduced dataset
null_lm <- lm(log(betaplasma) ~ 1, data = plasmaB_excl_pred)
null_sum <- summary(null_lm)
null_sum
logplasmaB_wocalories_lm
anova(null_lm, logplasmaB_wocalories_excl_lm)
logplasmaB_all_excl_lm <- lm(log(betaplasma) ~ bmi + age + calories + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, data = plasmaB_excl_pred)
anova(null_lm, logplasmaB_all_excl_lm)
logplasmaB_bmi_lm <- lm(log(betaplasma) ~bmi, data = plasmaB_excl_pred)

logplasmaB_age_lm <- lm(log(betaplasma) ~ age, data = plasmaB_excl_pred)
logplasmaB_age_lm

logplasmaB_calories_lm <- lm(log(betaplasma) ~ calories, data = plasmaB)
logplasmaB_calories_lm

BIC(null_lm, logplasmaB_bmi_lm)
logplasmaB_bmiage_lm <- lm(log(betaplasma) ~ bmi + age, data = plasmaB_excl_pred)

BIC(logplasmaB_bmi_lm, logplasmaB_bmiage_lm)
# Not really significant change

logplasmaB_bmicalories_lm <- lm(log(betaplasma) ~ bmi + calories, data = plasmaB_excl_pred)

BIC(logplasmaB_bmicalories_lm, logplasmaB_bmi_lm)
# nope, increase

logplasmaB_bmifat_lm <- lm(log(betaplasma) ~ bmi + fat, data = plasmaB_excl_pred)

BIC(logplasmaB_bmi_lm, logplasmaB_bmifat_lm)
# litle increase, maybe

logplasmaB_bmicholesterol_lm <- lm(log(betaplasma) ~ bmi + cholesterol, data = plasmaB_excl_pred)
BIC(logplasmaB_bmi_lm, logplasmaB_bmicholesterol_lm)
# litle increase, maybe

logplasmaB_bmifiber_lm <- lm(log(betaplasma) ~ bmi + fiber, data = plasmaB_excl_pred)
BIC(logplasmaB_bmi_lm, logplasmaB_bmifiber_lm)

logplasmaB_bmifiberalcohol_lm <- lm(log(betaplasma) ~ bmi + fiber +alcohol, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifiberalcohol_lm, logplasmaB_bmifiber_lm)

logplasmaB_bmifiberbetadiet_lm <- lm(log(betaplasma) ~ bmi + fiber + betadiet, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifiberbetadiet_lm, logplasmaB_bmifiber_lm)

logplasmaB_bmifibersmokstat_lm <- lm(log(betaplasma) ~ bmi + fiber + smokstat, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstat_lm, logplasmaB_bmifiber_lm)

logplasmaB_bmifibersmokstatsex_lm <- lm(log(betaplasma) ~ bmi + fiber + smokstat +sex, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstat_lm, logplasmaB_bmifibersmokstatsex_lm)

# FINAL
logplasmaB_bmifibersmokstatvituse_lm <- lm(log(betaplasma) ~ bmi + fiber + smokstat +vituse, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvituse_lm, logplasmaB_bmifibersmokstat_lm)
#

logplasmaB_bmifibersmokstatvituseage_lm <- lm(log(betaplasma) ~ bmi + fiber + smokstat +vituse +age, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvituseage_lm, logplasmaB_bmifibersmokstatvituse_lm )

logplasmaB_bmicholesterolfiber_lm <- lm(log(betaplasma) ~ bmi + cholesterol + fiber, data = plasmaB_excl_pred)
BIC(logplasmaB_bmicholesterolfiber_lm, logplasmaB_bmifibersmokstatvituse_lm)

logplasmaB_bmicholesterolfiberalcohol_lm <- lm(log(betaplasma) ~ bmi + cholesterol + fiber +alcohol, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvituse_lm, logplasmaB_bmicholesterolfiberalcohol_lm)

logplasmaB_bmicholesterolfiberbetadiet_lm <- lm(log(betaplasma) ~ bmi + cholesterol + fiber + betadiet, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvituse_lm, logplasmaB_bmicholesterolfiberbetadiet_lm)

logplasmaB_bmicholesterolfibersmokstat_lm <- lm(log(betaplasma) ~ bmi + cholesterol + fiber + smokstat, data = plasmaB_excl_pred)
logplasmaB_bmicholesterolfibersmokstat_lm
BIC(logplasmaB_bmifibersmokstatvituse_lm, logplasmaB_bmicholesterolfibersmokstat_lm)

logplasmaB_bmicholesterolfibersex_lm <- lm(log(betaplasma) ~ bmi + cholesterol + fiber +sex, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvituse_lm, logplasmaB_bmicholesterolfibersex_lm)

# SOMETHING HERE
logplasmaB_bmicholesterolfibervituse_lm <- lm(log(betaplasma) ~ bmi + cholesterol + fiber + vituse, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvituse_lm, logplasmaB_bmicholesterolfibervituse_lm)

# HERE
logplasmaB_bmicholesterolfibervituseage_lm <- lm(log(betaplasma) ~ bmi +age+ cholesterol + fiber + vituse, data = plasmaB_excl_pred)
BIC(logplasmaB_bmicholesterolfibervituseage_lm, logplasmaB_bmifibersmokstatvituse_lm)

logplasmaB_bmicholesterolfibervituseagefat_lm <- lm(log(betaplasma) ~ bmi +age+ cholesterol + fiber + vituse +fat, data = plasmaB_excl_pred)
BIC(logplasmaB_bmicholesterolfibervituseage_lm, logplasmaB_bmicholesterolfibervituseagefat_lm)


logplasmaB_bmicholesterolfibervitusefat_lm <- lm(log(betaplasma) ~ bmi + fat+cholesterol + fiber + vituse, data = plasmaB_excl_pred)
BIC(logplasmaB_bmicholesterolfibervitusefat_lm, logplasmaB_bmifibersmokstatvituse_lm)
logplasmaB_bmifibersmokstatvituse_lm

logplasmaB_bmifibersmokstatvitusecholesterol <- lm(log(betaplasma) ~ bmi + fiber + smokstat +vituse + cholesterol, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvitusecholesterol, logplasmaB_bmifibersmokstatvituse_lm)

logplasmaB_bmifibersmokstatvitusecholesterolage <- lm(log(betaplasma) ~ bmi + fiber + smokstat +vituse + cholesterol + age, data = plasmaB_excl_pred)
BIC(logplasmaB_bmifibersmokstatvitusecholesterolage, logplasmaB_bmifibersmokstatvituse_lm)
## 1st attempt: Final model: log(betaplasma) ~ bmi + fat+cholesterol + fiber + vituse
logplasmaB_finalmodel_lm <- logplasmaB_bmicholesterolfibervituseage_lm
logplasmaB_finalmodel_sum <- summary(logplasmaB_finalmodel_lm)
logplasmaB_finalmodel_sum
logplasmaB_wocalories_sum
BIC(logplasmaB_bmicholesterolfibervituseage_lm, logplasmaB_wocalories_lm)
cbind(beta = logplasmaB_finalmodel_lm$coefficients, 
      confint(logplasmaB_finalmodel_lm))

# logplasmaB_bmicholesterolfibervituseage_lm  dof= 8, BIC=679.7964




# 4(b) TRY AGAIN :(((((((()))))))), I AM CONFUSED
# REDUCED dataset, 1
# BIC since we supplied k=log(n)
null1_lm <- lm(log(betaplasma) ~ 1, data = plasmaB_excl_pred)
null1_sum <- summary(null1_lm)
null1_sum
all1_lm <- lm(log(betaplasma) ~ bmi + age + calories + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, data = plasmaB_excl_pred)
step(null1_lm, 
     scope = list(lower = null1_lm, upper = all1_lm),
     direction = "both",
     k = log(nobs(all1_lm)))
glimpse(plasmaB_excl_pred)
# log(betaplasma) ~ bmi + fiber + calories + vituse, data = plasmaB_excl_pred

reduced_lm <- lm(formula = log(betaplasma) ~ bmi + fiber + calories + vituse, 
                 data = plasmaB_excl_pred)
model_j_lm <- lm(log(betaplasma) ~ bmi + fiber + fat + vituse, data = plasmaB_excl_pred)
BIC(reduced_lm, model_j_lm)
reduced_sum <- summary(reduced_lm)
reduced_sum

# FULL dataset, 2
# BIC since we supplied k=log(n)
null2_lm <- lm(log(betaplasma) ~ 1, data = plasmaB)
null2_sum <- summary(null2_lm)
null2_sum
all2_lm <- lm(log(betaplasma) ~ bmi + age + calories + fat + cholesterol + fiber + alcohol + betadiet + smokstat + sex + vituse, data = plasmaB)
step(null2_lm, 
     scope = list(lower = null2_lm, upper = all2_lm),
     direction = "both",
     k = log(nobs(all2_lm)))
full_lm <- lm(formula = log(betaplasma) ~ bmi + fiber + cholesterol + vituse, 
              data = plasmaB)
full_sum <- summary(full_lm)
full_sum
# Full model: Cholesterol	Significant (p < 0.001), calories Removed from model
# Reduced: Cholesterol not in model, calories	added; significant (p < 0.001)


# 5 Fine-tuning the model
# vitusesometimes seems not significant
mutate(plasmaB_excl_pred, 
       vituseno = factor(vituse == "no",
                         levels = c(FALSE, TRUE),
                         labels = c("yes use", "no use"))) -> plasmaB
glimpse(plasmaB_excl_pred)

reduced_vituseno_lm <- lm(log(betaplasma) ~ bmi + fiber + calories + vituseno, 
                          data = plasmaB_excl_pred)
reduced_vituseno_sum <- summary(reduced_vituseno_lm)
reduced_vituseno_sum
anova(reduced_lm, reduced_vituseno_lm)
# There are no difference in these two models with categories vituse

reduced_vituseno_sum$adj.r.squared
reduced_sum$adj.r.squared
BIC(reduced_lm, reduced_vituseno_lm)
AIC(reduced_lm, reduced_vituseno_lm)
# BIC, AIC lower for vituseno

cbind(beta = reduced_vituseno_lm$coefficients, 
      confint(reduced_vituseno_lm))


# 5(b)
# Refit Model.1(b) and Model.2(b) on the reduced data set as well
bmi_lm <- lm(log(betaplasma) ~ bmi, data = plasmaB_excl_pred)
bmi_lm
# "never" as reference
plasmaB_excl_pred <- mutate(plasmaB_excl_pred, smokstat = relevel(smokstat, "never"))
smokstat_lm <- lm(log(betaplasma) ~ smokstat, data = plasmaB_excl_pred)
smokstat_lm

# Present a table with the number of β-parameters, the residual, standard deviation, the R2, adjusted R2, AIC and BIC 
# for the all five models (1(b), 2(b), 3(c), 4(b), and 5(a))
# Load necessary packages
library(broom)
library(purrr)

# Assume your five models are named as follows:
# model_1b, model_2b, model_3c, model_4b, model_5a
logplasmaB_wocalories_excl_lm
# Create a named list of models
models <- list(
  `Model 1(b)` = bmi_lm,
  `Model 2(b)` = smokstat_lm,
  `Model 3(c)` = logplasmaB_wocalories_excl_lm,
  `Model 4(b)` = reduced_lm,
  `Model 5(a)` = reduced_vituseno_lm
)

# Extract summary metrics
model_summaries <- imap_dfr(models, function(model, model_name) {
  sum_mod <- summary(model)
  
  tibble(
    Model = model_name,
    Parameters = length(coef(model)),
    Residual_SD = sum_mod$sigma,
    R2 = sum_mod$r.squared,
    Adjusted_R2 = sum_mod$adj.r.squared,
    AIC = AIC(model),
    BIC = BIC(model)
  )
})

bmi_sum <- summary(bmi_lm)
smokstat_sum <- summary(smokstat_lm)
logplasmaB_wocalories_excl_sum <- summary(logplasmaB_wocalories_excl_lm)
reduced_sum
model_summaries
print(model_summaries, n = Inf)
collect.R2AIC <- data.frame(
  nr = seq(1, 5),
  model = c("1(b)", "2(b)", "3(c)", "4(b)", "5(a)"),
  numb_of_paramters = c(length(coef(bmi_lm)),
                         length(coef(smokstat_lm)),
                         length(coef(logplasmaB_wocalories_excl_lm)),
                         length(coef(reduced_lm)),
                         length(coef(reduced_vituseno_lm))),
  sigma = c(bmi_sum$sigma,
            smokstat_sum$sigma,
            logplasmaB_wocalories_excl_sum$sigma,
            reduced_sum$sigma,
            reduced_vituseno_sum$sigma),
  R2 = c(bmi_sum$r.squared,
         smokstat_sum$r.squared,
         logplasmaB_wocalories_excl_sum$r.squared,
         reduced_sum$r.squared,
         reduced_vituseno_sum$r.squared),
  R2.adj = c(bmi_sum$adj.r.squared,
             smokstat_sum$adj.r.squared,
             logplasmaB_wocalories_excl_sum$adj.r.squared,
             reduced_sum$adj.r.squared,
             reduced_vituseno_sum$adj.r.squared))
collect.R2AIC |> arrange(desc(R2.adj))
ggplot(collect.R2AIC) +
  geom_point(aes(x = model, y = R2, color = "R2"), size = 3) + 
  geom_point(aes(x = model, y = R2.adj, color = "R2-adjusted"), size = 3) + 
  geom_line(aes(x = nr, y = R2, color = "R2"), linewidth = 1) +
  geom_line(aes(x = nr, y = R2.adj, color = "R2-adjusted"), 
            linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 1) +
  labs(title = "betaplasma: R2 and R2-adjusted",
       color = "") +
  ylab("R2") +
  theme(legend.position = "bottom")
# adj. r square highest for model 3, best
cbind(collect.R2AIC,
      AIC = AIC(bmi_lm, smokstat_lm, logplasmaB_wocalories_excl_lm, reduced_lm, reduced_vituseno_lm),
      BIC = BIC(bmi_lm, smokstat_lm, logplasmaB_wocalories_excl_lm, reduced_lm, reduced_vituseno_lm)) -> collect.R2AIC
glimpse(collect.R2AIC)
collect.R2AIC |>
  rename(df = AIC.df, AIC = AIC.AIC, BIC = BIC.BIC) |>
  mutate(BIC.df = NULL) -> collect.R2AIC
collect.R2AIC |> arrange(AIC)
# AIC: best for 3

collect.R2AIC |> arrange(BIC)
# BIC : Best for model 5
collect.R2AIC |> arrange(R2)

# model 5: best
# 3: highest adj. R^2, lowestAIC, number of paramters 13, second worst BIC
# 5: 2nd highest adj R^2, second lowest AIC, number of paramters 5, best BIC