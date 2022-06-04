###########################################################################################
###                                  SURVIVAL ANALYSIS                                  ###
###########################################################################################
library(survival)
library(survminer)
library(dplyr)
# Surv()    : creates a survival object
# survfit() : fits a survival curve using a formula
# survdiff(): log-rank test for differences in survival between two or more groups

data('lung')
data = lung

# time  : survival time 
# status: censoring status 1=censored, 2=dead

### ---------------------------------------------------------------------------------------
### Plot (reduced)
### ---------------------------------------------------------------------------------------
data_red    = head(lung)
data_red$ID = factor(seq(1:nrow(data_red)))
ggplot(data = data_red, aes(x=ID, y=time)) + 
  geom_bar(stat='identity', width=0.2) + 
  geom_point(aes(color=status, shape=as.factor(status)), size=6) + 
  coord_flip()

### ---------------------------------------------------------------------------------------
### Survival Object
### ---------------------------------------------------------------------------------------
# The function Surv(time, event) create a survival object
Surv(data$time, data$status)

### ---------------------------------------------------------------------------------------
### KAPLAN-MEIER
### ---------------------------------------------------------------------------------------
### |-------------------------------------------|
### | Kaplan-Meier estimator for survival curve |
### |-------------------------------------------|
fit = survfit(Surv(time, status) ~ 1, data=data)
summary(fit)

# n           : total number of subjects
# time        : the event time points on the curve (t=t*_j)
# n.risk      : the number of subjects at risk at time t
# n.event     : the number of events that occurred at time t
# n.censor    : the number of censored subjects, who exit the risk set at time t
# surv        : the kaplan-meier estimator for survival S(t)
# std.err     : the standard error for S(t)
# lower, upper: lower and upper confidence limits for the survival curve S(t), respectively
# cumhaz      : the cumulative hazard curve H(t) = - log(S(t))
# std.err     : the standard error for H(t)

### |----------------------|
### | Median Survival time |
### |----------------------|
# Time at which the survival probability (S(t)) is 0.5
median_St = fit$time[fit$surv<=0.5][1]
median_St
# or
surv_median(fit)
# or
print(fit)

### |-------------------------|
### | Kaplan-Meier curve plot |
### |-------------------------|
plot(fit, conf.int=T, xlab='Time', ylab='Survival Probability', main='Kaplan-Meier curve')
# or
ggsurvplot(fit, data = data,
           risk.table = TRUE, 
           risk.table.col = 'strata', 
           surv.median.line='hv', 
           ggtheme = theme_bw(),
           break.time.by = 90, 
           title = 'Kaplan-Meier curve')

### |--------------------------------------|
### | Cumulative incidence function /      |
### | Cumulative failure probability (CFP) |
### |--------------------------------------|
# It shows the cumulative probabilities of experiencing the event of interest and 
# it is computed as CFP(t) = P(T<t), so can be estimated as 1-S(t)
cumulative_incidence = 1 - fit$surv
ggsurvplot(fit, data = data,
           fun='event',
           risk.table = TRUE, 
           risk.table.col = 'strata', 
           surv.median.line='hv', 
           ggtheme = theme_bw(),
           break.time.by = 90, 
           title = 'Cumulative Incidence curve')

### |-----------------------------------|
### | Cumulative hazard function (H(t)) |
### |-----------------------------------|
# It can be interpreted as the cumulative force of mortality.
H = fit$cumhaz
ggsurvplot(fit, data = data,
           fun='cumhaz',
           risk.table = TRUE, 
           risk.table.col = 'strata', 
           ggtheme = theme_bw(),
           break.time.by = 90, 
           title = 'Cumulative Hazard curve')

### |-------|
### | Recap |
### |-------|
curves = data.frame('time'          = fit$time,
                    'Survival'      = fit$surv,
                    'Cum_incidence' = 1-fit$surv,
                    'Cum_hazard'    = fit$cumhaz)
head(curves)

### |------------------------------------|
### | Survival probability at given time |
### |------------------------------------|
t_0 = 180
summary(fit, times=t_0)

# Survival probability every six months
t_0 = seq(0, 365*3, 182.5)
summary(fit, times=t_0)

### ---------------------------------------------------------------------------------------
### KAPLAN-MEIER by a feature
### ---------------------------------------------------------------------------------------
### |-------------------------------------------|
### | Kaplan-Meier estimator for survival curve |
### |-------------------------------------------|
feature     = data$sex
fit.feature = survfit(Surv(time, status) ~ feature, data=data)
print(fit.feature)
summary(fit.feature)$table
summary(fit.feature)

### |---------------------------------------|
### | Kaplan-Meier curves plot by a feature |
### |---------------------------------------|
plot(fit.feature, conf.int=T, xlab='Time', ylab='Survival Probability',
     col=c('blue', 'red'), main='Kaplan-Meier curve')
legend('topright', legend=c('Feat_1', 'Feat_2'), lty=c(1,1), col=c('blue', 'red'))
# or
ggsurvplot(fit.feature, data = data,
           conf.int = T,
           risk.table = TRUE,
           risk.table.col = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(),
           break.time.by = 90,
           legend.labs = c("Feat_1", "Feat_2"), 
           legend.title = "Feature",
           palette = c("blue", "red"),
           title = 'Kaplan-Meier curves by a feature')

### |---------------------------------------------------|
### | Cumulative incidence function /                   |
### | Cumulative failure probability (CFP) by a feature |
### |---------------------------------------------------|
ggsurvplot(fit.feature, data = data,
           fun = 'event',
           conf.int = T,
           risk.table = TRUE,
           risk.table.col = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(),
           break.time.by = 90,
           legend.labs = c("Feat_1", "Feat_2"), 
           legend.title = "Feature",
           palette = c("blue", "red"),
           title = 'Cumulative Incidence curves by a feature')

### |------------------------------------------------|
### | Cumulative hazard function (H(t)) by a feature |
### |------------------------------------------------|
ggsurvplot(fit.feature, data = data,
           fun = 'cumhaz',
           conf.int = T,
           risk.table = TRUE,
           risk.table.col = "strata", 
           ggtheme = theme_bw(),
           break.time.by = 90,
           legend.labs = c("Feat_1", "Feat_2"), 
           legend.title = "Feature",
           palette = c("blue", "red"),
           title = 'Cumulative Hazard curves by a feature')

### |-----------------------------|
### | Log-rank test for a feature |
### |-----------------------------|
# Comparison of two or more survival curves.
# The null hypothesis is that there is no difference in survival between the two groups.
# The log rank test is a non-parametric test, which makes no assumptions about the survival
# distributions. The test compares the observed number of events in each group to what 
# would be expected if the null hypothesis were true (the survival curves were identical).
# The log rank statistic is approximately distributed as a chi-square test statistic.
log_rank = survdiff(Surv(time, status) ~ feature, data=data)
log_rank

# small p-value => the curves are different
# high p-value  => the curves are the same

ggsurvplot(fit.feature, data = data,
           conf.int = T,
           risk.table = TRUE,
           risk.table.col = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(),
           break.time.by = 90,
           legend.labs = c("Feat_1", "Feat_2"), 
           legend.title = "Feature",
           palette = c("blue", "red"),
           title = 'Kaplan-Meier curves by a feature',
           pval = T)

### |--------------|
### | Hazard Ratio |
### |--------------|
# To quantify the difference in the survivals we can compute the hazard ratio, i.e.
# the ratio between the death hazard of the first group vs the other one.
# The risk of death in the first group is hazard_ratio * the risk of death in the second.
observed     = log_rank$obs
expected     = log_rank$exp
hazard_ratio = (observed[1]/expected[1])/(observed[2]/expected[2])
hazard_ratio

# = 1: no effect 
# > 1: decrease in the hazard (increase of the survival) -> the numerator is a protective factor
# < 1: increase in the hazard (decrease of the survival) -> the numerator is a risk factor

### |------------------------------------|
### | Survival probability at given time |
### |------------------------------------|
t_0 = 180
summary(fit.feature, times=t_0)

# Survival probability every six months
t_0 = seq(0, 365*3, 182.5)
summary(fit.feature, times=t_0)

### ---------------------------------------------------------------------------------------
### COX PH MODEL - Univariate
### ---------------------------------------------------------------------------------------
# Compute Cox proportional-hazards regression models: coxph(formula, data, method)
# formula: linear model with a survival object as response variable
# data   : data frame containing the variables
# method : specifies how to handle ties
x1    = data$age
cox_m = coxph(Surv(time, status) ~ x1, data=data)
cox_m
summary(cox_m)

# 1) STATISTICAL SIGNIFICANCE
#    The column marked "z" gives the Wald statistic value. It corresponds to the
#    ratio of each regression coefficient to its standard error (z = coef/se(coef)).
#    The wald statistic evaluates, whether the beta coefficient of a given
#    variable is statistically significantly different from 0.
#
# 2) THE REGRESSION COEFFICIENTS
#    The second feature to note in the Cox model results is the the sign of the regression 
#    coefficients (coef). A positive sign means that the hazard (risk of death) is higher,
#    and thus the prognosis is worse, for subjects with higher values of that variable.
#
# 3) HAZARD RATIO & CONFIDENCE INTERVAL
#    The exponentiated coefficients (exp(coef) = exp(0.0187) = 1.019), also known as hazard
#    ratios, give the effect size of covariates. For example, the increase of 1 unit in the 
#    covariate increase the hazard of 1.9%. The summary output also gives upper and lower 
#    95% confidence intervals for the hazard ratio (exp(coef)).
#
# 4) GLOBAL STATISTICAL SIGNIFICANCE OF THE MODEL
#    Finally, the output gives p-values for three alternative tests for overall significance
#    of the model: The likelihood-ratio test, Wald test, and score logrank statistics.
#    These three methods are asymptotically equivalent. For large enough N, they will give
#    similar results. For small N, they may differ somewhat. The Likelihood ratio test has
#    better behavior for small sample sizes, so it is generally preferred.

### |------------|
### | Curve plot |
### |------------|
plot(survfit(cox_m, data=data), lwd=2, lty=1, xlab='Time', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()

### |----------------------------------------|
### | How the covariate influences the curve |
### |----------------------------------------|
# We take some values in the range of the covariate x1, for example:
x1_new = data.frame(x1 = c(50, 65, 80))
x1_new

fit_m = survfit(cox_m, newdata=x1_new)
fit_m

plot(fit_m, conf.int=T, col=c('red', 'green', 'blue'), lwd=2, lty=1, 
     xlab='Time', ylab='Survival Probablity', main='Adjusted Survival Probability')
grid()
legend('topright', c('x1=50', 'x1=65', 'x1=80'), lty=c(1,1,1), lwd=c(2,2,2), col=c('red', 'green', 'blue'))

### ---------------------------------------------------------------------------------------
### COX PH MODEL - Multivariate
### ---------------------------------------------------------------------------------------
x1 = data$age
x2 = as.factor(data$sex)
x3 = data$ph.karno
x4 = data$wt.loss

cox_mult = coxph(Surv(time, status) ~ x1 + x2 + x3 + x4, data=data)
cox_mult
summary(cox_mult)

# 1. The last three p-values indicate the significance of the model. These tests evaluate
#    the null hypothesis that all of the betas are 0 (low p-val => the model is valid).
# 2. We can look at the p-values of the singular covariates to see which are significant.
# 3. We can look at the HR for the covariates (exp(coef)) to see if they're >/< 1.
#    3.1. Consider the categorical variable, say x2, with exp(coef) = exp(0.514) = 1.67.
#         This means that being the (displayed) group of x2 increases the hazard by a factor
#         of 1.67, which means 67%. In this case being of that group is a bad prognostic.
#    3.2. Consider the continuous variable, say x3, with exp(coef) = exp(-0.013) = 0.987.
#         This means that, holding the other variables constant, a higher value fo x3 is
#         associated with good survival. 

### |------------------------|
### | Plot for Hazard Ratios |
### |------------------------|
# It doesn't read the coxph with x1, .., xk
data$sex = as.factor(data$sex)
cox_plot = coxph(Surv(time, status) ~ age + sex + ph.karno + wt.loss, data=data)
ggforest(cox_plot, data=data)

### |------------|
### | Curve plot |
### |------------|
plot(survfit(cox_mult, data=data), lwd=2, lty=1, xlab='Time', ylab='Survival Probability',
     main='Baseline estimated survival probability')
grid()

### |---------------------------|
### | Cox Model Goodness of fit |
### |---------------------------|
### Martingale residuals
# We want the martingale residuals to have 0 mean along time
ggcoxdiagnostics(cox_mult, type='martingale')

### Deviance residuals
# These are a normalized transform of the martingale residuals.
# They should be symmetrically distributed about 0 with standard deviation of 1.
#   * Positive values correspond to individuals that died too soon compare to the expected survival times
#   * Negative values correspond to individuals that lieved too long
ggcoxdiagnostics(cox_mult, type='deviance')

### Shoenfeld residuals
# These residuals represent the difference between the observed covariate and the expected
# given the risk set at that time. They should be flat, centered about 0.
# A plot showing a non-random pattern against time is evidence of violation of the assumptions.
ggcoxdiagnostics(cox_mult, type='schoenfeld')

### Log(-log(KM)) - only for categorical variables
# We plot log(-log(KM(t))) vs. t or log(t) and look for parallelism.
# If we find parallelism then the assumptions are satisfied.
x_cat = as.factor(data$sex)
km    = survfit(Surv(time, status) ~ x_cat, data=data)
plot(km, fun='cloglog')

### |--------------|
### | COX.ZPH test |
### |--------------|
# The function cox.zph() is used to test the proportional hazards assumption for each
# covariate included in a Cox refression model fit. For each covariate, cox.zph() correlates 
# the corresponding set of scaled Schoenfeld residuals with time, to test for independence
# between residuals and time. Also, it performs a global test for the model as a whole.
#
# The proportional hazard assumption is supported by a non-significant relationship between
# residuals and time, and refuted by a significant relationship.
# 
# Test for PH using scaled Schoenfeld test for PH:
#    H0: Hazards are proportional
#    H1: Hazards are NOT proportional
# cox.zph() return tests for each X and for the global model
test.ph = cox.zph(cox_mult)
test.ph