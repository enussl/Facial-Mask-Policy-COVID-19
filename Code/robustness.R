# (I) OVERVIEW

# This script runs all the separate robustness checks using the helper function that are sourced
# from helperfunctions.R.

################################################################################

# (II) ENVIRONMENT AND PACKAGES

# Empty environment
rm(list = ls())

# Reproducibility (bootstrap and more)
set.seed(42) 

# Set working directory to the root node of folder structure
#setwd(".\\Mask_Project\\Final")
setwd("C:/Users/eminu/OneDrive/Desktop/Facial-Mask-Policy-COVID-19")

# Read helper functions. Note that they are in the same directory. Add the corresponding path
# otherwise.
#source(".\\Scripts\\helperfunctions.R")
source("./Code/helperfunctions.R")

################################################################################

# (III) OVERVIEW OF ROBUSTNESS CHECKS

#     - (1) Adjust time-level to cluster on for the cluster robust standard error estimators (month) 
#     - (2) Add additional information variables to case growth to match the heuristics of Chernozhukov (2020)
#     - (3) Merge the half-cantons of Switzerland together
#     - (4) Change the parameter r.infovar to d = 14 to examine the timing of the information dynamics for the r-value
#     - (5) Change the timing regarding case-growth by estimating d = 21, which is motivated by Chernozhukov (2020)
#     - (6) Examine the results after removing observations for which cooks distance > 4/n
#     - (7) Briefly touch on the daily models as done in Chernozhukov (2020)
#     - (8) Briefly touch on the alternative method to compute the weekly case growth to stress the similarity with the R value (only in text form)
#     - (9) Short period where the cantonal heterogeneity was most present (August to October)
#     - (10) Double Machine Learning approaches

# Below, we report the point estimate, std.error and p.value concerning W for all these alterations.

# (1) Clustering on months. We obviously only do it for the fixed effects models and we further look at the total effects although
# there is no change compared to the direct effect. 

# Result matrix
results.month = matrix(NA, nrow = 2*2, ncol = 3)
colnames(results.month) = c("estimate", "std.error", "p.val")

# Parameters
resp.poss = c("median_R_mean", "casegrowth")
cluster.poss = c("month")
frequency = "weekly"
infovar = FALSE
type.effect = "total"
model = "within"

# Allocation
j = 1

for (response in resp.poss) {
  
  # Store formula
  #formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
  
  # Store data; change configuration of data processing here to examine lag, shift and r.infovar. Add monthly indicators.
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  mo = rep(c(1,2,3,4), each = 4)
  mo = mo[-c(16,16)]
  data$month = rep(mo, 26)
  
  # Fixed effects model. We use the feols package as the plm package only allows clustering at the same levek as fixed effects, i.e weeks
  # instead of months.
  
  # fit = plm(formula, 
  #           data = data, model = model, 
  #           index = c("X.Canton_3", "X.oneweek"), 
  #           effect = "twoways")
  
  if (response == "casegrowth") {
    fit = feols(Y ~ X.casegrowthlag_we + X.casegrowthlag1_we + X.sre000d0_we + 
                  X.tre200d0_we + X.ure200d0_we + X.restGatherings_we + X.cancEvents_we + 
                  X.workClosing2a_we + W + X.ferien_we | X.Canton_3 + X.oneweek, data = data,
                cluster = c("month"))
  } else if (response == "median_R_mean") {
    fit = feols(Y ~ X.median_R_mean.lag_we + X.median_R_mean.lag1_we + X.sre000d0_we + 
                  X.tre200d0_we + X.ure200d0_we + X.restGatherings_we + X.cancEvents_we + 
                  X.workClosing2a_we + W + X.ferien_we | X.Canton_3 + X.oneweek, data = data,
                cluster = c("month"))
  }
  
  
  # res = coeftest(fit, vcov = function(x) 
  #   clubSandwich::vcovCR(x, type = "CR0", cluster = c(data$month)))
  
  # Results
  results.month[j,1] = fit[["coeftable"]]["W",1]
  results.month[j,2] = fit[["coeftable"]]["W",2]
  results.month[j,3] = fit[["coeftable"]]["W",4]
  j = j + 1 
}

# Do the double clustering separately
for (response in resp.poss) {
  
  # Store data; change configuration of data processing here to examine lag, shift and r.infovar. Add monthly indicators.
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  mo = rep(c(1,2,3,4), each = 4)
  mo = mo[-c(16,16)]
  data$month = rep(mo, 26)
  
  # Fixed effects model. We use the feols package as the plm package only allows clustering at the same levek as fixed effects, i.e weeks
  # instead of months.
  
  # fit = plm(formula, 
  #           data = data, model = model, 
  #           index = c("X.Canton_3", "X.oneweek"), 
  #           effect = "twoways")
  
  if (response == "casegrowth") {
    fit = feols(Y ~ X.casegrowthlag_we + X.casegrowthlag1_we + X.sre000d0_we + 
                  X.tre200d0_we + X.ure200d0_we + X.restGatherings_we + X.cancEvents_we + 
                  X.workClosing2a_we + W + X.ferien_we | X.Canton_3 + X.oneweek, data = data,
                cluster = c("month", "X.Canton_3"))
  } else if (response == "median_R_mean") {
    fit = feols(Y ~ X.median_R_mean.lag_we + X.median_R_mean.lag1_we + X.sre000d0_we + 
                  X.tre200d0_we + X.ure200d0_we + X.restGatherings_we + X.cancEvents_we + 
                  X.workClosing2a_we + W + X.ferien_we | X.Canton_3 + X.oneweek, data = data,
                cluster = c("month", "X.Canton_3"))
  }
  
  
  # res = coeftest(fit, vcov = function(x) 
  #   clubSandwich::vcovCR(x, type = "CR0", cluster = c(data$month)))
  
  # Results
  results.month[j,1] = fit[["coeftable"]]["W",1]
  results.month[j,2] = fit[["coeftable"]]["W",2]
  results.month[j,3] = fit[["coeftable"]]["W",4]
  j = j + 1 
  
}

# Naming
results.month = as.data.frame(results.month, row.names = c("Month Clustered FE R", "Month Clustered FE Casegrowth", 
                                                           "Canton-Month Clustered FE R", "Canton-Month Clustered FE Casegrowth"))

################################################################################

# (2) Additional information variables as in Chernozhukov (2020)

# Run the models for response = case-growth and additional information variables. Note that we can only
# run it using the random effects approach as some of the additional information variables are at a national
# level which is not compatible with the fixed effects approach due to time fixed effects. We use the HC3 standard errors
# as always with the random effects models.

# Result matrix
results.infovar = matrix(NA, nrow = 2, ncol = 3)
colnames(results.infovar) = c("estimate", "std.error", "p.val")

# Parameters
type.effect.poss = c("direct", "total")
frequency = "weekly"
response = "casegrowth"
model = "random"
infovar = TRUE

# Allocation
j = 1

# Run the models
for (type.effect in type.effect.poss) {
  
  results.infovar[j,1] = estimation(frequency = frequency, response = response, model = model, infovar = infovar, type.effect = type.effect)$results[[1]]["W",1]
  results.infovar[j,2] = estimation(frequency = frequency, response = response, model = model, infovar = infovar, type.effect = type.effect)$results[[1]]["W",2]
  results.infovar[j,3] = estimation(frequency = frequency, response = response, model = model, infovar = infovar, type.effect = type.effect)$results[[1]]["W",4]
  j = j + 1
}

# Naming
results.infovar = as.data.frame(results.infovar, row.names = c("Add. Infovars RE Casegrowth direct", "Add. Infovars RE Casegrowth total"))


# (3) Merge the half cantons of Switzerland together as one can argue that they do not constitute separate observations due to
# being strongly interconnected and the small number of inhabitants and the geographical proximity. We alter the construction of the 
# data to make that happen. Note that the cantons are not ordered for the own estimator for the standard errors which is not a
# problem as we do not use it for this robustness check. We use time-canton standard errors for the fixed effects models
# as they are usually most conservative and HC3 standard errors for the random effects models as always.

# Result matrix
results.halfcanton = matrix(NA, nrow = 4, ncol = 3)
colnames(results.halfcanton) = c("estimate", "std.error", "p.val")

# Parameters
resp.poss = c("casegrowth", "median_R_mean")
model.poss = c("within", "random")
type.effect = "total"
frequency = "weekly"
infovar = FALSE

# Allocation
j = 1

# Run the models
for (model in model.poss) {
  
  for (response in resp.poss) {
    
    # Store formula
    formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
    
    # Store data; change configuration of data processing here to examine lag, shift and r.infovar
    data = data.prep.halfcantons(shift = 14, response = response, r.infovar = 21)$data
    
    # Models
    fit = plm(formula, 
              data = data, model = model, 
              index = c("X.Canton_3", "X.oneweek"), 
              effect = "twoways")
    
    if (model == "within") {
      
      # Covariance matrix
      res = coeftest(fit, vcov = function(x) 
        clubSandwich::vcovCR(x, type = "CR0", cluster = c(data$X.oneweek, data$X.Canton_3)))
      
      # Results
      results.halfcanton[j,1] = res["W",1]
      results.halfcanton[j,2] = res["W",2]
      results.halfcanton[j,3] = res["W",4]
      
    } else {
      
      # Covariance matrix
      res = coeftest(fit, vcov = function(x) 
        plm::vcovHC(x, method = "white1", type = "HC3"))
      
      # Results
      results.halfcanton[j,1] = res["W",1]
      results.halfcanton[j,2] = res["W",2]
      results.halfcanton[j,3] = res["W",4]
    }
    j = j + 1
  }
}

# Naming
results.halfcanton = as.data.frame(results.halfcanton, row.names = c("Halfcant. FE Casegrowth", "Halfcant. FE R", "Halfcant. RE Casegrowth", 
                                                                     "Halfcant. RE R"))

# Half-cantons for DFE
multiple_split_half = function(response, frequency, infovar, type.effect) {
  
  # response: response variable in {median_R_mean, casegrowth}
  # frequency: frequency in {daily, weekly}
  # infovar: add additional information variables
  # type.effect: direct or total to include behavior variables or not
  
  # Output: bias.corrected estimate of mask policy
  
  # Data
  data = data.prep.halfcantons(response = response, r.infovar = 21, shift = 14)$data
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = "within", infovar = infovar, type.effect = type.effect)
  
  # Number of splits along the cross-section
  s = 20
  
  # Fit non-debiased model daily or weekly and save time and unit for sampling later
  if (frequency == "daily") {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.day"), effect = "twoways")
    time   = as.double(data$X.day)
    unit   = as.double(data$X.Canton_3)
    
  } else {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.oneweek"), effect = "twoways")
    time   = as.double(data$X.oneweek)
    unit   = as.double(data$X.Canton_3)
  }
  
  # Start computation 
  across = 0 * coef(uc)
  
  for (k in 1:s) {
    
    # Sampling process
    sample1 = sample(unique(unit), ceiling(length(unique(unit))/2), replace = FALSE)
    
    subsample1 = ((unit %in% sample1) & (time <= median(time))) |
      ((!(unit %in% sample1)) & (time > median(time)))
    
    subsample2 = ((unit %in% sample1) & (time > median(time))) |
      ((!(unit %in% sample1)) & (time <= median(time)))
    
    if (frequency == "daily") {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.day"), effect = "twoways")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.day"),effect = "twoways")
      
      
    } else {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"), effect = "twoways")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"),effect = "twoways")
      
    }
    across = across + ((coef(cross1) + coef(cross2))/2)/s
  }
  
  # Average cross over corrected
  acbc = 2 * coef(uc) - across
  return(acbc)
}


# Run it for the DFE approach
multiple_split_half(response = "median_R_mean", frequency = "weekly", infovar = FALSE, type.effect = "total")
multiple_split_half(response = "casegrowth", frequency = "weekly", infovar = FALSE, type.effect = "total")

################################################################################

# (4) Change the parameter of r.infovar to 14. We look at the most conservative estimator for the standard errors (canton-time), which
# corresponds to the 4th entry of the estimation output. There is no real difference between the direct and the total effect so we look
# at the direct effect. For the random effects model, we use the HC3 estimator as always.

# Result matrix
results.rinfovar = matrix(NA, nrow = 2, ncol = 3)
colnames(results.rinfovar) = c("estimate", "std.error", "p.val")

# Parameters
model.poss = c("within", "random")
response = "median_R_mean"
frequency = "weekly"
type.effect = "total"
infovar = FALSE

# Allocation
j = 1

# Run the models
for (model in model.poss) {
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
  
  # Store data; change configuration of data processing here to examine lag, shift and r.infovar
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 14, frequency = frequency)$data
  
  # Models
  fit = plm(formula, 
            data = data, model = model, 
            index = c("X.Canton_3", "X.oneweek"), 
            effect = "twoways")
  
  if (model == "within") {
    
    # Covariance matrix
    res = coeftest(fit, vcov = function(x) 
      clubSandwich::vcovCR(x, type = "CR0", cluster = c(data$X.oneweek, data$X.Canton_3)))
    
    # Results
    results.rinfovar[j,1] = res["W",1]
    results.rinfovar[j,2] = res["W",2]
    results.rinfovar[j,3] = res["W",4]
    
  } else {
    
    # Covariance matrix
    res = coeftest(fit, vcov = function(x) 
      sandwich::vcovHC(x, type = "HC3"))
    
    # Results
    results.rinfovar[j,1] = res["W",1]
    results.rinfovar[j,2] = res["W",2]
    results.rinfovar[j,3] = res["W",4]
  }
  j = j + 1
}

# Naming
results.rinfovar = as.data.frame(results.rinfovar, row.names = c("R-Infovar 14 FE R", "R-Infovar 14 RE R"))

################################################################################

# (5) Change the timing regarding case growth. We use d = 21 as explained above. We again use 
# the canton-time standard errors for the fixed effects model and the HC3 estimator for the random effects model. We again
# look at the direct effects but there is no difference compared to the total effects.

# Result matrix
results.casetiming = matrix(NA, nrow = 2, ncol = 3)
colnames(results.casetiming) = c("estimate", "std.error", "p.val")

# Parameters
model.poss = c("within", "random")
response = "casegrowth"
frequency = "weekly"
type.effect = "total"
infovar = FALSE

# Allocation
j = 1

for (model in model.poss) {
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
  
  # Store data; change configuration of data processing here to examine lag, shift and r.infovar
  data = data.prep(lag = 7, shift = 21, response = response, r.infovar = 21, frequency = frequency)$data
  
  # Fixed effects model
  fit = plm(formula, 
            data = data, model = model, 
            index = c("X.Canton_3", "X.oneweek"), 
            effect = "twoways")
  
  if (model == "within") {
    
    # Covariance matrix
    res = coeftest(fit, vcov = function(x) 
      clubSandwich::vcovCR(x, type = "CR0", cluster = c(data$X.oneweek, data$X.Canton_3)))
    
    # Results
    results.casetiming[j,1] = res["W",1]
    results.casetiming[j,2] = res["W",2]
    results.casetiming[j,3] = res["W",4]
    
  } else {
    
    # Covariance matrix
    res = coeftest(fit, vcov = function(x) 
      sandwich::vcovHC(x, type = "HC3"))
    
    # Results
    results.casetiming[j,1] = res["W",1]
    results.casetiming[j,2] = res["W",2]
    results.casetiming[j,3] = res["W",4]
  }
  j = j + 1
}

# Naming
results.casetiming = as.data.frame(results.casetiming, row.names = c("Case-Infovar 21 FE Casegrowth", "Case-Infovar 21 RE Casegrowth"))


# Do the timing for the de-biased fixed effects approach
multiple_split_timing = function(response, frequency, infovar, type.effect) {
  
  # response: response variable in {median_R_mean, casegrowth}
  # frequency: frequency in {daily, weekly}
  # infovar: add additional information variables
  # type.effect: direct or total to include behavior variables or not
  
  # Output: bias.corrected estimate of mask policy
  
  # Data
  data = data.prep.halfcantons(response = response, r.infovar = 14, shift = 21)$data
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = "within", infovar = infovar, type.effect = type.effect)
  
  # Number of splits along the cross-section
  s = 20
  
  # Fit non-debiased model daily or weekly and save time and unit for sampling later
  if (frequency == "daily") {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.day"), effect = "twoways")
    time   = as.double(data$X.day)
    unit   = as.double(data$X.Canton_3)
    
  } else {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.oneweek"), effect = "twoways")
    time   = as.double(data$X.oneweek)
    unit   = as.double(data$X.Canton_3)
  }
  
  # Start computation 
  across = 0 * coef(uc)
  
  for (k in 1:s) {
    
    # Sampling process
    sample1 = sample(unique(unit), ceiling(length(unique(unit))/2), replace = FALSE)
    
    subsample1 = ((unit %in% sample1) & (time <= median(time))) |
      ((!(unit %in% sample1)) & (time > median(time)))
    
    subsample2 = ((unit %in% sample1) & (time > median(time))) |
      ((!(unit %in% sample1)) & (time <= median(time)))
    
    if (frequency == "daily") {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.day"), effect = "twoways")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.day"),effect = "twoways")
      
      
    } else {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"), effect = "twoways")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"),effect = "twoways")
      
    }
    across = across + ((coef(cross1) + coef(cross2))/2)/s
  }
  
  # Average cross over corrected
  acbc = 2 * coef(uc) - across
  return(acbc)
}

multiple_split_timing(response = "casegrowth", frequency = "weekly", infovar = FALSE, type.effect = "total")
multiple_split_timing(response = "median_R_mean", frequency = "weekly", infovar = FALSE, type.effect = "total")



################################################################################

# (6) Remove observations according to the cooks distance rule of thumb being CD > 4/n. We also report the change in R^2
# induced via the removal of these observations. This only works for the fixed effects approach. We again use the direct effect due 
# to there being no difference when comparing to the total effect. We use the canton-time standard errors. Note that the degrees of freedom 
# are not correct yet.

# Result matrix
results.cooks = matrix(NA, nrow = 2, ncol = 7)
colnames(results.cooks) = c("estimate", "std.error", "p.val", "r.squared", "deleted", "diff.point", "diff.rsquared")

# Parameters
resp.poss = c("casegrowth", "median_R_mean")
model = "within"
frequency = "weekly"
type.effect = "total"
infovar = FALSE

# Allocation
j = 1

for (response in resp.poss) {
  
  # Two-way demean data to use in linear model (lm) 
  data = tw.demean(data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data,
                   response = response, frequency = frequency)
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
  
  # Run fixed effects regression via lm
  fit.original = lm(formula, data = data)
  
  # Save orginal R^2
  r.squared.original = r.sq(fit = fit.original, frequency = frequency, response = response)
  
  # Compute cooks distance and determine influential observations
  cooksD = cooks.distance(fit.original)
  influential = cooksD[(cooksD > 4/nrow(data))]
  influential.05 = cooksD[(cooksD > 0.5)]
  
  # Remove outlaying observations
  names_of_influential = names(influential)
  outliers = data[names_of_influential,]
  data_without_outliers = data %>% 
    anti_join(outliers)
  
  # Run fixed effects regression via lm on reduced data
  fit.reduced = lm(formula, data = data_without_outliers) 
  
  # Covariance matrix
  res = coeftest(fit.reduced, vcov = function(x)
    sandwich::vcovCL(x, cluster = ~ X.oneweek + X.Canton_3))
  
  # Save new R^2
  r.squared.reduced = r.sq(fit = fit.reduced, frequency = frequency, response = response)
  
  # Results
  results.cooks[j,1] = res["W",1]
  results.cooks[j,2] = res["W",2]
  results.cooks[j,3] = res["W",4]
  results.cooks[j,4] = r.squared.reduced
  results.cooks[j,5] = length(influential)
  results.cooks[j,6] = summary(fit.reduced)$coefficients["W",1] - summary(fit.original)$coefficients["W",1]
  results.cooks[j,7] = r.squared.reduced - r.squared.original
  j = j + 1
}

# Naming
results.cooks = as.data.frame(results.cooks, row.names = c("CooksD FE Casegrowth", "CooksD FE R"))

# Removing observations from panel data is not justified easily. We therefore also look at a different method of dealing with
# influential observations through robust regression as opposed to deleting observations depending on the cooks distance. This is further
# justified as the distribution of our residuals show heavy-tailed patterns. For more information, see the summary under
# https://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch13.pdf.

# Result matrix
results.robust = matrix(NA, nrow = 2, ncol = 2)
colnames(results.robust) = c("estimate", "std.error")

# Parameters
resp.poss = c("casegrowth", "median_R_mean")
model = "within"
frequency = "weekly"
type.effect = "total"
infovar = FALSE

# Allocation
j = 1

for (response in resp.poss) {
  
  # Two-way demean data to use in linear model (lm) 
  data = tw.demean(data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data,
                   response = response, frequency = frequency)
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
  
  # Run robust linear regression via rlm
  fit = rlm(formula, data = data)
  
  # Results
  results.robust[j,1] = summary(fit)$coef["W",1]
  results.robust[j,2] = summary(fit)$coef["W",2]
  j = j + 1
}

# Naming
results.robust = as.data.frame(results.robust, row.names = c("Robust FE Casegrowth", "Robust FE R"))


################################################################################

# (7) & (8) Quickly discuss in the paper. See results-05 in the data folder for 
# the results.

################################################################################

# (9) Short period where the treatment was most heterogeneous. We take the Canton-Time estimator for the fixed effects
# and the HC3 estimator for the random effects. We look at the direct effects. The weekly approach does not work for the random effects
# due to a shortage of observations.


# Result matrix
results.shortperiod = matrix(NA, nrow = 4, ncol = 3)
colnames(results.shortperiod) = c("estimate", "std.error", "p.val")

# Parameters
resp.poss = c("casegrowth", "median_R_mean")
model.poss = c("within", "random")
frequency = "weekly"
type.effect = "total"
infovar = FALSE
startdate = "2020-08-21"
enddate = "2020-10-18"

# Allocation
j = 1

for (model in model.poss) {
  
  for (response in resp.poss) {
    
    # Store formula
    formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
    
    # Store data; change configuration of data processing here to examine lag, shift and r.infovar
    data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency,
                     startdate = startdate, enddate = enddate)$data
    
    # Fixed effects model
    fit = plm(formula, 
              data = data, model = model, 
              index = c("X.Canton_3", "X.oneweek"), 
              effect = "twoways")
    
    if (model == "within") {
      
      # Covariance matrix
      res = coeftest(fit, vcov = function(x) 
        clubSandwich::vcovCR(x, type = "CR0", cluster = c(data$X.oneweek, data$X.Canton_3)))
      
      # Results
      results.shortperiod[j,1] = res["W",1]
      results.shortperiod[j,2] = res["W",2]
      results.shortperiod[j,3] = res["W",4]
      
    } else {
      
      # Covariance matrix
      res = coeftest(fit, vcov = function(x) 
        sandwich::vcovHC(x, type = "HC3"))
      
      # Results
      results.shortperiod[j,1] = res["W",1]
      results.shortperiod[j,2] = res["W",2]
      results.shortperiod[j,3] = res["W",4]
    }
    
    j = j + 1
  }
}

# Naming
results.shortperiod = as.data.frame(results.shortperiod, row.names = c("Short Per. FE Daily Casegrowth", "Short Per. FE Daily R",
                                                                       "Short Per. FE Weekly Casegrowth", "Short Per. FE Weekly R"))

# Do that for the DFE
multiple_split_short = function(response, frequency, infovar, type.effect) {
  
  # response: response variable in {median_R_mean, casegrowth}
  # frequency: frequency in {daily, weekly}
  # infovar: add additional information variables
  # type.effect: direct or total to include behavior variables or not
  
  # Output: bias.corrected estimate of mask policy
  
  # Data
  data = data.prep.halfcantons(response = response, r.infovar = 21, shift = 14, startdate = "2020-08-21",
                               enddate = "2020-10-18")$data
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = "within", infovar = infovar, type.effect = type.effect)
  
  # Number of splits along the cross-section
  s = 20
  
  # Fit non-debiased model daily or weekly and save time and unit for sampling later
  if (frequency == "daily") {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.day"), effect = "twoways")
    time   = as.double(data$X.day)
    unit   = as.double(data$X.Canton_3)
    
  } else {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.oneweek"), effect = "twoways")
    time   = as.double(data$X.oneweek)
    unit   = as.double(data$X.Canton_3)
  }
  
  # Start computation 
  across = 0 * coef(uc)
  
  for (k in 1:s) {
    
    # Sampling process
    sample1 = sample(unique(unit), ceiling(length(unique(unit))/2), replace = FALSE)
    
    subsample1 = ((unit %in% sample1) & (time <= median(time))) |
      ((!(unit %in% sample1)) & (time > median(time)))
    
    subsample2 = ((unit %in% sample1) & (time > median(time))) |
      ((!(unit %in% sample1)) & (time <= median(time)))
    
    if (frequency == "daily") {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.day"), effect = "twoways")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.day"),effect = "twoways")
      
      
    } else {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"), effect = "twoways")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"),effect = "twoways")
      
    }
    across = across + ((coef(cross1) + coef(cross2))/2)/s
  }
  
  # Average cross over corrected
  acbc = 2 * coef(uc) - across
  return(acbc)
}

multiple_split_short(response = "median_R_mean", frequency = "weekly", infovar = FALSE, type.effect = "total")


################################################################################

# (10) Double Machine Learning. See results.csv file for the results.

################################################################################

# Construct the main table and safe the results as a csv file
results.robustnesscheck = rbind(results.month[c(1,3,5,6),], results.infovar, results.halfcanton,
                                results.rinfovar, results.casetiming,
                                results.cooks[,1:3], results.shortperiod)

write.csv(results.robustnesscheck,".\\Data\\results.robustnesschecks.csv", row.names = TRUE)



conf_int = function(point_est, sd, df) {
  return(round(c(point_est, point_est - qt(0.975, df)*sd, point_est + qt(0.975, df)*sd),2))
}


