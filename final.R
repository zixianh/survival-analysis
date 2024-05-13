rm(list=ls())
library(survival)
library(tidyverse)
library(MASS)
library(ggplot2)
library(survminer)

pbcseq
data(pbc, package="survival")

## (b)
S = Surv(pbcseq$day, pbcseq$futime, pbcseq$status)


## (c)
### (i)
# The Kaplan-Meier estimator plus the empirical distribution function of X is 1. 
# S(t) = 1 - F(t)

### (ii)
# Build the Kaplan-Meier estimator for each level of the covariates drug and sex
fit <- survfit(Surv(futime, status) ~ trt + sex, data = pbcseq)

# Plot the estimator for all of the levels of each covariate on the same figure
ggsurvplot(fit, xscale="d_y", xlab="Year", break.time.by = 365.25*2,
           legend = c(0.85, 0.75), ylim=c(0,1))


### (iii)
S2 = Surv(pbcseq$futime,  pbcseq$status)
survdiff(Surv(futime, status) ~ trt + sex, data = pbcseq)
# Since the p-value is 5e-08 which is quite smaller than 0.05, we will reject 
# the null hypothesis and conclude that sex and drug have significantly 
# different Kaplan-Meier survival functions.
# The log-rank test statistic follows chi-square distribution, 
# and it is calculated by comparing the observed and expected number of 
# events in each group at each time point.


## (e)
### (i)
pbcseq$log_alk_pho = log(pbcseq$alk.phos)
pbcseq$log_platelet = log(pbcseq$platelet)

### (ii)
colSums(is.na(pbcseq))
str(pbcseq)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[is.na(uniqv)==FALSE]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#### mode
pbcseq$ascites <- replace_na(pbcseq$ascites, getmode(pbcseq$ascites))
pbcseq$hepato <- replace_na(pbcseq$hepato, getmode(pbcseq$hepato))
pbcseq$spiders <- replace_na(pbcseq$spiders, getmode(pbcseq$spiders))
pbcseq$chol <- replace_na(pbcseq$chol, getmode(pbcseq$chol))

#### mean
pbcseq$log_alk_pho <- replace_na(pbcseq$log_alk_pho, mean(pbcseq$log_alk_pho, na.rm=T))
pbcseq$log_platelet <- replace_na(pbcseq$log_platelet, mean(pbcseq$log_platelet, na.rm=T))


### (iii)
initial_model <- coxph(S2 ~ trt + bili, data = pbcseq)
summary(initial_model)


### (iv)
Scope = list(upper=~trt+bili+age+sex+ascites+hepato+spiders+edema+chol+albumin+ast+protime+stage+log_alk_pho+log_platelet, 
             lower=~trt+bili)
final_m <- stepAIC(initial_model, direction="forward", 
                       scope=Scope, trace = FALSE)
final_m$anova


### (v)
summary(final_m)
anova(final_m)
# The type of the drug is used to predict the survival of the patients
# since its coefficient is significant at 5% significance level. 


### (vi)
anova(final_m)
summary(final_m)
cox.zph(final_m)
# the global chi-square test (on 9 degrees of freedom) has a p-value 
# significantlt less than 0.05.


### (vii)
CPH1 <- coxph(S2 ~ trt+bili+age+log_alk_pho+stage+hepato+chol+sex+spiders, data=pbcseq)
CPH2 <- coxph(S2 ~ trt+bili+age+log_alk_pho+stage+hepato+chol+spiders, data=pbcseq)

extractAIC(CPH1)
extractAIC(CPH2)
# The model with sex has smaller AIC value which indicates the model with sex
# is better. 


### (viii)
fit1 <- survfit(CPH1, data=pbcseq)
fit2 <- survfit(CPH2, data=pbcseq)

ggsurvplot(fit1, xscale="d_y", xlab="Year", break.time.by = 365.25*2,
           legend = c(0.85, 0.75), ylim=c(0,1))

ggsurvplot(fit2, xscale="d_y", xlab="Year", break.time.by = 365.25*2,
           legend = c(0.85, 0.75), ylim=c(0,1))

