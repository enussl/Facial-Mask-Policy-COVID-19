# (I) OVERVIEW

# This script contains all the helper functions used for this project. Note that
# there is a detailed description of the functions and the dependencies in the read-me
# file on github.

################################################################################

# (II) ENVIRONMENT AND PACKAGES

# Empty environment
rm(list = ls())

# Set working directory to the root node of folder structure
#setwd(".\\Mask_Project\\Final")
setwd("C:/Users/eminu/OneDrive/Desktop/Facial-Mask-Policy-COVID-19")

# Load used packages
library(readr)
library(dplyr)
library(tidyverse)
library(stargazer)
library(grf)  
library(ranger)
library(lfe)
library(corrplot)
library(fastDummies)
library(faraway)
library(xtable)
library(lme4)
library(plm)
library(lmtest)
library(sandwich)
library(texreg)
library(car)
library(tidyquant)
library(modelr)
library(DoubleML)
library(mlr3)
library(collapse)
library(clusterSEs)
library(boot)
library(faraway)
library(nlme)
library(dynpanel)
library(abind)
library(blockmatrix)
library(dotwhisker)
library(tikzDevice)
library(extrafont)
library(rcompanion)
library(RColorBrewer)
library(ggcorrplot)
library(plot.matrix)
library(clubSandwich)

################################################################################

# (III) FUNCTIONS

panellag = function(x, i, t, lag=1) {
  
  # x: Vector or matrix to get lags of; sorted by id and time
  # i: unit index
  # t: time index
  # lag: How many time periods to lag. Can be negative if leading values are desired.
  
  # Output: returns lagged copy of x
  
  x = as.vector(x)
  i = as.vector(i)
  t = as.vector(t)
  
  if (!identical(order(i,t),1:length(i))) {
    stop("inputs not sorted.")
  }
  if (is.matrix(x)) {
    return(apply(x,MARGIN=2,FUN=function(c) { panel.lag(c,i,t,lag) }))
  }
  if (length(i) != length(x) || length(i) != length(t) ) {
    stop("Inputs not same length")
  }
  if (lag>0) {
    x.lag = x[1:(length(x)-lag)]
    x.lag[i[1:(length(i)-lag)]!=i[(1+lag):length(i)] ] = NA
    x.lag[(t[1:(length(i)-lag)]+lag)!=t[(1+lag):length(i)] ] = NA
    val = (c(rep(NA,lag),x.lag))
  } else if (lag<0) {
    lag = abs(lag)
    x.lag = x[(1+lag):length(x)]
    x.lag[i[1:(length(i)-lag)]!=i[(1+lag):length(i)] ] = NA
    x.lag[t[1:(length(i)-lag)]+lag!=t[(1+lag):length(i)] ] = NA
    val = (c(x.lag,rep(NA,lag)))
  } else { # lag=0
    return (x)
  }
  if (inherits(x,"Date") & class(val)=="numeric") {
    stopifnot(0==as.numeric(as.Date("1970-01-01")))
    val = as.Date(val, origin="1970-01-01")
  }
  return(val)
}

paneldiff = function(x, id, t, lag = 1) {
  
  # x: Vector or matrix to get differences of; sorted by id and time
  # id: unit index
  # t: time index
  # lag: How many time periods to lag. Can be negative if leading values are desired.
  
  # Output: returns differenced copy of x
  
  return(x - panellag(x,id,t,lag = lag))
}

panelma = function(x, id, t, len) {
  
  # x: vector to compute moving average of; sorted by id and time
  # id: unit index
  # t: time index
  # len: length of MA (corresponds to h for us)
  
  # Output: returns average of x[t:ken(t-len)]
  
  stopifnot(len>=0)
  ma  = rep(0, length(x))
  for(i in 0:len) {
    ma  = ma + panellag(x, id,t, lag=i)
  }
  return(ma/(len+1))
}

lognozero = function(x, minval = 0.1) {
  
  # x: Vector or matrix to apply function on
  
  # Output: returns max{log(0.1), log(y)} for y in the real numbers to avoid log(0)
  
  x[x < minval] = minval
  return(log(x))
}

dlogd = function(x, id, t, lag, minval = 0.1) {
  
  # x: Vector or matrix to apply function on; sorted by id and time
  # id: unit index
  # t: time index
  # lag: How many time periods to lag. Can be negative if leading values are desired.
  
  # Output: returns response regarding casegrowth as defined in the paper 
  
  dx  = paneldiff(x, id, t, lag = lag)
  dx[dx < minval] = minval
  dlogdx = paneldiff(log(dx),id,t,lag = lag)
} 

data.prep = function(lag, shift, response, startdate = "2020-07-06",
                     enddate = "2020-10-18", r.infovar, frequency) {
  
  # 2020-12-21
  
  # lag: smoothing of covariates of lag days
  # shift: shifting of response when using casegrowth as defined 
  # response: in {casegrowth, median_R_mean}
  # startdate, enddate: to be adjusted for different time periods
  # r.infovar: lag length of information variables
  # frequency: in {daily, weekly}
  
  # Output: data frame used for the estimation (Y, W, X)
  
  # Note that the baseline mask policy of mandatory mask wearing was introduced nationally on "2020-07-06" and the 
  # first vaccinations took place on "2020-12-21". This is the reason for our sample selection.
  
  # Read in data
  DataCovid = read.csv(".\\Data\\DataCovid.csv") 
  dat.new = read.csv(".\\Data\\policy_stuff.csv", header = TRUE, sep = ",")
  
  # Merge data on common days
  data = merge(DataCovid, dat.new)
  
  # Delete unnecessary variables to might lead to conflicts
  data = data %>%
    dplyr::select(-c(deaths_1d, recovered_1d, tests_1d, hosp_1d, value.total_1d, sre000d0_1d, tre200d0_1d, ure200d0_1d, O65P, restGatheringsCH)) %>%
    rename(testingPolicy = testing_policyCH)
  
  # Transform percentage to be a real percentage 
  data = data %>%
    mutate(percentage_age = percentage_age*100)
  
  # Get the daily date and weekly date as integer. Further construct disjoint weeks
  data$day = cut.Date(as.Date(data$datum), breaks = "1 day", labels = FALSE) 
  data$oneweek = cut.Date(as.Date(data$datum), breaks = "1 week", start.on.monday = TRUE, labels = FALSE) 
  
  # Extract which value of oneweek corresponds to startdate, enddate to apply the time selection
  # via the vector oneweek later
  start.day = which(data$datum == startdate)[1]
  start.week = data[start.day, "oneweek"]
  end.day   = which(data$datum == enddate)[1]
  end.week = data[end.day, "oneweek"]
  
  if (response == "casegrowth" && frequency == "daily") {
    
    # Construct target variable and lead target variable
    data$casegrowth = dlogd(x = data$sumTotal, id = data$geoRegion, t = as.Date(data$datum), lag = lag)
    data$casegrowth = panellag(x = data$casegrowth, i = data$geoRegion, t = as.Date(data$datum), lag = - shift)
    
    # Apply the MA-operator to (X, W)
    to.be.smoothed = c("grocery_and_pharmacy", "transit_stations", "Number.of.transactions", "workplaces",
                       "schoolClosing", "cancEvents", "restGatherings","testingPolicy", "workClosing2a",
                       "sre000d0", "tre200d0", "ure200d0", "mask_treat")
    
    for (v in to.be.smoothed) {
      data[,v] = panelma(x = data[,v], id = data$geoRegion, t = as.Date(data$datum), len = lag - 1)
    }
    
    # Lagged target variable
    data$casegrowthlag = dlogd(x = data$sumTotal, id = data$geoRegion, t = as.Date(data$datum), lag = lag)
    
    # log Delta Cit with Cit; new cases
    data$logdiffcases = lognozero(paneldiff(x = data$sumTotal, id = data$geoRegion, t = as.Date(data$datum), lag = lag))
    
    # Construct national aggregates
    sumTotalCH = aggregate(data$sumTotal, by = list(data$datum), FUN  = sum)
    sumTotalCH = rename(sumTotalCH, datum = Group.1)
    sumTotalCH = rename(sumTotalCH, sumTotalCH = x)
    data = merge(sumTotalCH, data)
    data = data %>%
      arrange(geoRegion)
    
    # Lagged target variable on national level; works with id argument
    data$casegrowthlag.national = dlogd(x = data$sumTotalCH, id = data$geoRegion, t = as.Date(data$datum), lag = lag)
    
    # log Delta Cit; same as above but national level
    data$logdiffcases.national = lognozero(paneldiff(x = data$sumTotalCH, id = data$geoRegion, t = as.Date(data$datum), lag = lag))
    
    # Delta log T_it; already new tests so we do not have to apply the difference operator
    data$difflogtests = dlogd(x = data$tests, id = data$geoRegion, t = as.Date(data$datum), lag = lag)
    
    # Delta log transactions (one day rate)
    data$Number.of.transactions[data$Number.of.transactions < 0.1] = 0.1
    data$transactiongrowth = paneldiff(x = log(data$Number.of.transactions), id = data$geoRegion, t = as.Date(data$datum), lag = 1)
    
    # Correct period
    data = data %>%
      filter(datum >= startdate & datum <= enddate)
    
    # Change the ordering of cantons to adjust for the ordering of the adjacency matrix: we work with Canton_3 as the
    # cantonal indicator later so we order the data along the list target where we input the ordering of the adjacency matrix
    # used later so the data is ordered.
    target = c(2, 4, 5, 11, 3, 9, 17, 16, 22, 24, 20, 12, 8, 10, 6, 1, 15, 13, 14, 23, 18, 26, 21, 19, 7, 25)
    data = data %>% 
      slice(order(factor(Canton_3, levels = target)))
    
    # Construct X daily
    X = data[,c("casegrowthlag", "logdiffcases", "casegrowthlag.national", "logdiffcases.national", "difflogtests",
                "percentage_age", "Density", "population",
                "schoolClosing", "restGatherings", "cancEvents",
                "testingPolicy", "workClosing2a",
                "grocery_and_pharmacy", "transit_stations", "transactiongrowth", "workplaces",
                "sre000d0", "tre200d0", "ure200d0", "ferien", "Canton_3", "day")]
    
    W = data$mask_treat
    Y = data$casegrowth
    
  }
  
  if (response == "casegrowth" && frequency == "weekly") {
    
    # Construct target variable and lagged target variable. We do that by computing the new cases per canton and week. We can then
    # compute the log-growth rate using that series, which constructs 7 times that same growth rate per canton and week. When aggregating, 
    # we can then take the mean over disjoint weeks to obtain the desired variable. We do that for all growth rates.
    
    # Construct national aggregates
    sumTotalCH = aggregate(data$sumTotal, by = list(data$datum), FUN  = sum)
    sumTotalCH = rename(sumTotalCH, datum = Group.1)
    sumTotalCH = rename(sumTotalCH, sumTotalCH = x)
    data = merge(sumTotalCH, data)
    data = data %>%
      arrange(geoRegion)
    
    # Get new cases per day in each canton. Do it for the new cases on a national level as well.
    data = data %>%
      group_by(geoRegion) %>%
      mutate(NewCases = c(0,diff(sumTotal)),
             NewCases.national = c(0, diff(sumTotalCH))) %>%
      ungroup()
    
    data = data %>%
      relocate(sumTotal, .before = NewCases)
    
    # We now have all variables constant across oneweek and canton which need to be constant. The other ones
    # are aggregated to weekly data by simply averaging the values per week and canton. Note that we take the sum
    # for transactions and tests as we need the new number of transactions and tests per canton and week. These variables 
    # already come supplied as new transactions and tests so we do not have to treat them as the cases where we had to explicitely 
    # calculate the new cases as done in the above code.
    data = data %>%
      group_by(Canton_3, oneweek) %>%
      summarize(newcases_we           = sum(NewCases),
                newcases.national_we  = sum(NewCases.national),
                newtransactions_we    = sum(Number.of.transactions),
                newtests_we           = sum(tests),
                percentage_age_we = mean(percentage_age),
                Density_we        = mean(Density),
                population_we     = mean(population),
                grocery_and_pharmacy_we   = mean(grocery_and_pharmacy),
                transit_stations_we       = mean(transit_stations),
                workplaces_we             = mean(workplaces),
                ferien_we                 = mean(ferien), 
                sre000d0_we               = mean(sre000d0),
                tre200d0_we               = mean(tre200d0),
                ure200d0_we               = mean(ure200d0),
                schoolClosing_we          = mean(schoolClosing), 
                restGatherings_we         = mean(restGatherings), 
                cancEvents_we             = mean(cancEvents), 
                testingPolicy_we          = mean(testingPolicy), 
                workClosing2a_we          = mean(workClosing2a), 
                mask_treat_we             = mean(mask_treat)) %>%
      ungroup()
    
    # Using newcases_we, newcases.national_we, newtransactions_we and newtests_we, we can now compute the weekly growth rates.
    
    # Create the response variable and lead response variable. 
    # Check if values are below the threshold of 0.1 and set it to 0.1 if so. Does not happen
    data$newcases_we[data$newcases_we < 0.1] = 0.1
    data$casegrowth_we = paneldiff(x = log(data$newcases_we), id = data$Canton_3, t = as.Date(data$oneweek), l = 1)
    data$casegrowth_we = panellag(x = data$casegrowth_we, i = data$Canton_3, t = as.Date(data$oneweek), lag = -shift/7)
    
    # Lagged target variable
    data$casegrowthlag_we =  panellag(x = data$casegrowth_we, i = data$Canton_3, t = data$oneweek, lag = 2)
    data$casegrowthlag1_we = panellag(x = data$casegrowth_we, i = data$Canton_3, t = data$oneweek, lag = 1)
    
    # log Delta Cit
    data$logdiffcases_we = lognozero(paneldiff(x = data$newcases_we, id = data$Canton_3, t = as.Date(data$oneweek), lag = 1))
    
    # Lagged target variable on national level; works with id argument
    data$newcases.national_we[data$newcases.national_we < 0.1] = 0.1
    data$casegrowthlag.national_we = paneldiff(x = log(data$newcases.national_we), id = data$Canton_3, t = as.Date(data$oneweek), l = 1)
    
    # log Delta Cit; same as above but national level
    data$logdiffcases.national_we = lognozero(paneldiff(x = data$newcases.national_we, id = data$Canton_3, t = as.Date(data$oneweek), lag = 1))
    
    # Delta log Tit; already new tests so we do not have to apply the difference operator
    data$newtests_we[data$newtests_we < 0.1] = 0.1
    data$difflogtests_we = paneldiff(x = log(data$newtests_we), id = data$Canton_3, t = as.Date(data$oneweek), l = 1)
    
    # Delta log transactions 
    data$newtransactions_we[data$newtransactions_we < 0.1] = 0.1
    data$transactiongrowth_we = paneldiff(x = log(data$newtransactions_we), id = data$Canton_3, t = as.Date(data$oneweek), l = 1)
    
    # Correct period
    data = data %>%
      filter(oneweek >= start.week & oneweek <= end.week)
    
    # Change the ordering of cantons to adjust for the ordering of the adjacency matrix: we work with Canton_3 as the
    # cantonal indicator later so we order the data along the list target where we input the ordering of the adjacency matrix
    # used later so the data is ordered.
    target = c(2, 4, 5, 11, 3, 9, 17, 16, 22, 24, 20, 12, 8, 10, 6, 1, 15, 13, 14, 23, 18, 26, 21, 19, 7, 25)
    data = data %>% 
      slice(order(factor(Canton_3, levels = target)))
    
    # Construct X weekly
    X = data[,c("casegrowthlag_we", "casegrowthlag1_we", "logdiffcases_we", "casegrowthlag.national_we", "logdiffcases.national_we", "difflogtests_we",
                "percentage_age_we", "Density_we", "population_we",
                "schoolClosing_we", "restGatherings_we", "cancEvents_we",
                "testingPolicy_we", "workClosing2a_we",
                "grocery_and_pharmacy_we", "transit_stations_we", "transactiongrowth_we", "workplaces_we",
                "sre000d0_we", "tre200d0_we", "ure200d0_we", "ferien_we", "Canton_3", "oneweek")]
    W = data$mask_treat_we
    Y = data$casegrowth_we
    
  }
  
  if (response == "median_R_mean" && frequency == "daily") {
    
    # lagged target variable
    data$median_R_mean.lag = panellag(x = data$median_R_mean, i = data$geoRegion, t = as.Date(data$datum), lag = r.infovar)
    #data$median_R_mean.lag1 = panellag(x = data$median_R_mean, i = data$geoRegion, t = as.Date(data$datum), lag = 7)
    
    # Delta log transactions (one day rate)
    data$Number.of.transactions[data$Number.of.transactions < 0.1] = 0.1
    data$transactiongrowth = paneldiff(x = log(data$Number.of.transactions), id = data$geoRegion, t = as.Date(data$datum), lag = 1)
    
    # Correct period
    data = data %>%
      filter(datum >= startdate & datum <= enddate)
    
    # Change the ordering of cantons to adjust for the ordering of the adjacency matrix: we work with Canton_3 as the
    # cantonal indicator later so we order the data along the list target where we input the ordering of the adjacency matrix
    # used later so the data is ordered.
    target = c(2, 4, 5, 11, 3, 9, 17, 16, 22, 24, 20, 12, 8, 10, 6, 1, 15, 13, 14, 23, 18, 26, 21, 19, 7, 25)
    data = data %>% 
      slice(order(factor(Canton_3, levels = target)))
    
    # Construct X daily
    X = data[,c("median_R_mean.lag", 
                "percentage_age", "Density", "population",
                "schoolClosing", "restGatherings", "cancEvents",
                "testingPolicy", "workClosing2a",
                "grocery_and_pharmacy", "transit_stations", "transactiongrowth", "workplaces",
                "sre000d0", "tre200d0", "ure200d0", "ferien", "Canton_3", "day")]
    
    W = data$mask_treat
    Y = data$median_R_mean
    
  }
  
  if (response == "median_R_mean" && frequency == "weekly") {
    
    # Take the mean over the disjoint weeks for (Y,W,X)
    data = data %>%
      group_by(Canton_3, oneweek) %>%
      summarize(median_R_mean_we      = mean(median_R_mean),
                percentage_age_we = mean(percentage_age),
                Density_we        = mean(Density),
                population_we     = mean(population),
                grocery_and_pharmacy_we   = mean(grocery_and_pharmacy),
                transit_stations_we       = mean(transit_stations),
                newtransactions_we        = sum(Number.of.transactions),
                workplaces_we             = mean(workplaces),
                ferien_we                 = mean(ferien), 
                sre000d0_we               = mean(sre000d0),
                tre200d0_we               = mean(tre200d0),
                ure200d0_we               = mean(ure200d0),
                schoolClosing_we          = mean(schoolClosing), 
                restGatherings_we         = mean(restGatherings), 
                cancEvents_we             = mean(cancEvents), 
                testingPolicy_we          = mean(testingPolicy),
                workClosing2a_we          = mean(workClosing2a), 
                mask_treat_we             = mean(mask_treat)) %>%
      ungroup()
    
    # Delta log transactions 
    data$newtransactions_we[data$newtransactions_we < 0.1] = 0.1
    data$transactiongrowth_we = paneldiff(x = log(data$newtransactions_we), id = data$Canton_3, t = as.Date(data$oneweek), l = 1)
    
    # Lagged target variable
    data$median_R_mean.lag_we = panellag(x = data$median_R_mean_we, i = data$Canton_3, t = data$oneweek, lag = r.infovar/7)
    data$median_R_mean.lag1_we = panellag(x = data$median_R_mean_we, i = data$Canton_3, t = data$oneweek, lag = 1)
    
    # Correct period
    data = data %>%
      filter(oneweek >= start.week & oneweek <= end.week)
    
    # Change the ordering of cantons to adjust for the ordering of the adjacency matrix: we work with Canton_3 as the
    # cantonal indicator later so we order the data along the list target where we input the ordering of the adjacency matrix
    # used later so the data is ordered.
    target = c(2, 4, 5, 11, 3, 9, 17, 16, 22, 24, 20, 12, 8, 10, 6, 1, 15, 13, 14, 23, 18, 26, 21, 19, 7, 25)
    data = data %>% 
      slice(order(factor(Canton_3, levels = target)))
    
    # Construct X weekly
    X = data[,c("median_R_mean.lag_we", "median_R_mean.lag1_we",
                "percentage_age_we", "Density_we", "population_we",
                "schoolClosing_we", "restGatherings_we", "cancEvents_we",
                "testingPolicy_we", "workClosing2a_we",
                "grocery_and_pharmacy_we", "transit_stations_we", "transactiongrowth_we", "workplaces_we",
                "sre000d0_we", "tre200d0_we", "ure200d0_we", "ferien_we", "Canton_3", "oneweek")]
    
    W = data$mask_treat_we
    Y = data$median_R_mean_we
    
  }
  
  # Return  
  return(list(data = data.frame(X = X, Y = Y, W = W)))
  
} # Data preparation

construct.formula = function(response, frequency, model, infovar, type.effect) {
  
  # response: in {casegrowth, median_R_mean}
  # frequency: in {daily, weekly}
  # model: in {pooling, within, random}
  # infovar: include additional info variables for case growth approach
  # type.effect: direct or total to include behavior variables or not
  
  # Output: depending on the three arguments, we generate a different formula, which we output here. This is e.g due
  # to time-invariant covariates, which we want to include for the RE but no for the FE approach.
  
  if (frequency == "daily") {
    
    if (response == "median_R_mean") {
      
      formula = Y ~ X.median_R_mean.lag + X.sre000d0 + X.tre200d0 + X.ure200d0 +
        X.restGatherings + X.cancEvents + X.workClosing2a + W +
        X.ferien +  X.transactiongrowth 
      
      if (model == "pooling" | model == "random") {
        
        formula = update(formula, "~.+ X.population + X.Density + X.percentage_age + X.schoolClosing +
                         X.testingPolicy + X.grocery_and_pharmacy + X.transit_stations + X.workplaces")
        
      } # Pooling and RE
      
    } else {
      
      formula = Y ~ X.casegrowthlag + X.sre000d0 + X.tre200d0 + X.ure200d0 + 
        X.restGatherings  + X.cancEvents + X.workClosing2a + W + 
        X.ferien + X.transactiongrowth 
      
      if (model == "pooling" | model == "random") {
        
        formula = update(formula, "~.+ X.population + X.Density + X.percentage_age + X.schoolClosing +
                         X.testingPolicy + X.grocery_and_pharmacy + X.transit_stations + X.workplaces")
        
      } # Pooling and RE
      
    } 
    
    # Additional information variables
    if (infovar == TRUE && response == "casegrowth") {
      
      formula = update(formula, "~.+ X.logdiffcases + X.casegrowthlag.national + X.logdiffcases.national + X.difflogtests")
      
    } # Info-variables
    
    # Adjust the formula for direct vs. total effect estimation for the pooled and the random effects options
    if ((model == "pooling" | model == "random") && type.effect == "total") {
      
      # Find the positions of the behavior variables in the formula
      remove.vars = which(labels(terms(formula)) %in% c("X.workplaces", "X.transit_stations",
                                                        "X.grocery_and_pharmacy", "X.transactiongrowth"))
      
      # Drop the variables and construct the new formula
      formula = drop.terms(terms(formula), dropx = remove.vars, keep.response = TRUE)
      formula = reformulate(attr(formula, "term.labels"), formula[[2]])
    }
    
    # Adjust the formula for direct vs. total effect estimation for the fixed effects option
    if (model == "within" && type.effect == "total") {
      
      # Find the positions of the behavior variables in the formula
      remove.vars = which(labels(terms(formula)) %in% c("X.transactiongrowth"))
      
      # Drop the variables and construct the new formula
      formula = drop.terms(terms(formula), dropx = remove.vars, keep.response = TRUE)
      formula = reformulate(attr(formula, "term.labels"), formula[[2]])
    }
    
  } else {
    
    if (response == "median_R_mean") {
      
      formula = Y ~ X.median_R_mean.lag_we + X.median_R_mean.lag1_we + X.sre000d0_we + X.tre200d0_we + X.ure200d0_we +
        X.restGatherings_we + X.cancEvents_we + X.workClosing2a_we + W +
        X.ferien_we + X.transactiongrowth_we 
      
      if (model == "pooling" | model == "random") {
        
        formula = update(formula, "~.+ X.population_we + X.Density_we + X.percentage_age_we + X.schoolClosing_we +
                         X.testingPolicy_we + X.grocery_and_pharmacy_we + X.transit_stations_we + X.workplaces_we")
      }
      
    } else {
      
      formula = Y ~ X.casegrowthlag_we + X.casegrowthlag1_we + X.sre000d0_we + X.tre200d0_we + X.ure200d0_we + 
        X.restGatherings_we + X.cancEvents_we + X.workClosing2a_we + W + 
        X.ferien_we + X.transactiongrowth_we 
      
      if (model == "pooling" | model == "random") {
        
        formula = update(formula, "~.+ X.population_we + X.Density_we + X.percentage_age_we + X.schoolClosing_we +
                         X.testingPolicy_we + X.grocery_and_pharmacy_we + X.transit_stations_we + X.workplaces_we")
        
      }
      
    }
    
    if (infovar == TRUE && response == "casegrowth") {
      
      formula = update(formula, "~.+ X.logdiffcases_we + X.casegrowthlag.national_we + X.logdiffcases.national_we + X.difflogtests_we")
      
    } # Info-variables
    
    # Adjust the formula for direct vs. total effect estimation for the pooled and the random effects options
    if ((model == "pooling" | model == "random") && type.effect == "total") {
      
      # Find the positions of the behavior variables in the formula
      remove.vars = which(labels(terms(formula)) %in% c("X.workplaces_we", "X.transit_stations_we",
                                                        "X.grocery_and_pharmacy_we", "X.transactiongrowth_we"))
      
      # Drop the variables and construct the new formula
      formula = drop.terms(terms(formula), dropx = remove.vars, keep.response = TRUE)
      formula = reformulate(attr(formula, "term.labels"), formula[[2]])
    }
    
    # Adjust the formula for direct vs. total effect estimation for the fixed effects option
    if (model == "within" && type.effect == "total") {
      
      # Find the positions of the behavior variables in the formula
      remove.vars = which(labels(terms(formula)) %in% c("X.transactiongrowth_we"))
      
      # Drop the variables and construct the new formula
      formula = drop.terms(terms(formula), dropx = remove.vars, keep.response = TRUE)
      formula = reformulate(attr(formula, "term.labels"), formula[[2]])
    }
    
  } # Frequency
  
  # Return formula
  return(formula)
}

tw.demean = function(data, response, frequency) {
  
  # data: data to two-way demean
  # response: response variable
  # frequency: in {daily, weekly} 
  
  # Output: data but two-way demeaned as described in the paper. We need it for the Hansen covariance matrix estimator
  # as well as for the informal estimator where we need the two-way demeaned y, X. The Frisch-Waugh-Lovell theorem allows us
  # to compute the demeaning in that fashion (see https://stats.stackexchange.com/questions/200358/within-transformation-in-fixed-effect-regression-model).
  
  # Discriminate between frequencies
  if (frequency == "daily") {
    
    # Two-way demean everything apart from index and order such that X.Canton_3 and X.day are at the end of data
    vars.demean = data %>%
      dplyr::select(-c("X.Canton_3", "X.day"))
    data = data %>%
      relocate(X.Canton_3, .after = last_col()) %>%
      relocate(X.day, .after = X.Canton_3)
    
    # Two way demean
    for (k in 1:length(vars.demean)) {
      
      fit = lm(data[,k] ~ factor(data$X.day) + factor(data$X.Canton_3))
      res = fit$residuals 
      
      # Put the residuals into the data frame
      data[,k] = fit$residuals
      
    } # vars demean
    
  } else {
    
    # Two-way demean everything apart from index and order such that X.Canton_3 and X.oneweek are at the end of data
    vars.demean = data %>%
      dplyr::select(-c("X.Canton_3", "X.oneweek"))
    data = data %>%
      relocate(X.Canton_3, .after = last_col()) %>%
      relocate(X.oneweek, .after = X.Canton_3)
    
    # Two way demean
    for (k in 1:length(vars.demean)) {
      
      fit = lm(data[,k] ~ factor(data$X.oneweek) + factor(data$X.Canton_3))
      res = fit$residuals 
      
      # Put the residuals into the data frame
      data[,k] = fit$residuals
      
    } # vars demean
    
  }
  
  # Return demeaned data
  return(data = data)
}

M.Andrews = function(rho, T){
  
  # rho: determined in sd.Thompson.Hansen; corresponds to estimated AR(1) coefficient on S component of covariance matrix
  # T: number of time points in data
  
  # Output: choice of lag truncation parameter as described in Hansen (2022)
  
  nom = sum(rho^2/(1-rho)^4)
  
  denom = sum((1-rho^2)^2/(1-rho)^4)
  
  M = 1.8171*(nom/denom)^(1/3)*T^(1/3)
  
  # Return M
  return(M)
}

sd.Thompson.Hansen = function(y, X, N.cant, beta.hat, M.T=2){
  
  # y: response vector dim dim(y) = NT*1 
  # X: design matrix with dim(X) = NT*p and all time points of a canton are in a row. Further, X contains the intercept.
  # n.cantons: number of unit clusters (26 cantons)
  # beta.hat: vector of estimated coefficients with dim(beta.hat) = p*1 
  # M.T: number of lags in Thompson, default = 2
  
  # Output: covariance matrix of estimated coefficients for Hansen method, Thompson method as well a inferred choice of M
  
  # Time dimension
  T.time = nrow(X)/N.cant
  
  # First index is for time, second index is for canton, third for variable
  X.array = array(X, dim = c(T.time,N.cant,ncol(X))) 
  
  # Residuals
  U = y - X%*%beta.hat 
  
  # Array of residuals along unit and time
  U.array = array(U,dim = c(T.time, N.cant, 1))
  
  # R component from Hansen (2022)
  Ri = sapply(1:N.cant, function(i){
    colSums(X.array[,i,]*U.array[,i,1])
  }) 
  
  # S component from Hansen (2022)
  St = sapply(1:T.time, function(t){
    colSums(X.array[t,,]*U.array[t,,1])
  }) 
  
  # Gm component from Hansen (2022)
  Gm = lapply(1:(T.time-1), function(m){
    Gm.1 = lapply(1:(T.time-m), function(t){
      St[,t]%*%t(St[,t+m])
    })
    
    Reduce("+", Gm.1)
    
  })
  
  # Hm component from Hansen (2022)
  Hm = lapply(1:(T.time-1), function(m){
    Hm.1 = lapply(1:N.cant, function(i){
      
      Hm.2 = lapply(1:(T.time-m), function(t){
        
        cbind(X.array[t,i,]*U.array[t,i,1])%*%X.array[t+m,i,]*U.array[t+m,i,1]
      })
      
      Reduce("+", Hm.2)  
    })
    Reduce("+", Hm.1)  
  })
  
  # Combine the first three parts
  part1 = Reduce("+",lapply(1:N.cant, function(i){
    Ri[,i]%*%t(Ri[,i])  }))
  
  part2 = Reduce("+",lapply(1:T.time, function(t){
    St[,t]%*%t(St[,t])  }))
  
  part3 = Reduce("+",  lapply(1:N.cant, function(i){
    
    p3.1 = lapply(1:(T.time), function(t){
      
      cbind(X.array[t,i,])%*%X.array[t,i,]*U.array[t,i,1]^2
    })
    
    Reduce("+", p3.1)  
  }))
  
  # Part 4 from Thompson (Bruce (2022))
  part4.T = Reduce("+",lapply(1:M.T, function(m){
    
    Gm[[m]] + t(Gm[[m]]) - Hm[[m]]- t(Hm[[m]])
  }))
  
  # Meat from Thompson
  Omega.Thompson = 1/(N.cant*T.time)^2*(part1 + part2 - part3 + part4.T)
  
  Q = 1/(N.cant*T.time)*Reduce("+",  lapply(1:N.cant, function(i){
    
    Q.1 = lapply(1:(T.time), function(t){
      
      cbind(X.array[t,i,])%*%X.array[t,i,]
    })
    
    Reduce("+", Q.1)  
  }))
  
  # Covariance matrix from Thompson
  Var.Thompson = solve(Q)%*%Omega.Thompson%*%solve(Q)
  
  # Compute rho for choice of M as determined in M.Andrews function
  rho = unlist(lapply(1:ncol(X), function(j){
    
    Y.rho = St[j,2:ncol(St)]
    X.rho = St[j,1:(ncol(St)-1)]
    
    fit.rho = lm(Y.rho ~-1+ X.rho)
    rho = coef(fit.rho)
    
  }))
  
  # Determine M
  M = M.Andrews(rho, T = T.time)
  if (M > T.time) {
    M = T.time
  }
  
  # Part 4 from Hansen
  part4.H = Reduce("+",lapply(1:round(M), function(m){
    
    (1-m/(M+1))*(Gm[[m]] + t(Gm[[m]]) - Hm[[m]]- t(Hm[[m]]))
  }))
  
  # Meat from Hansen
  Omega.Hansen = 1/(N.cant*T.time)^2*(part1+ part2 - part3 + part4.H)
  
  # Covariance matrix from Hansen
  Var.Hansen = solve(Q)%*%Omega.Hansen%*%solve(Q)
  
  # Return objects
  return(list(Var.Hansen= Var.Hansen, Var.Thompson=Var.Thompson, M = M))
}

create.A.init = function(){
  
  # Output: adjacency matrix of swiss cantons with the same ordering as our data
  
  # Names
  cantons = c("GE", "VD", "FR", "NE", "VS", "BE", "JU", "SO", "BL", "BS", "AG", 
              "LU", "OW", "NW", "UR", "TI", "ZG", "SZ", "GL", "SG", "ZH", "SH",
              "AR", "AI",  "GR", "TG")
  
  # Initialize
  A.init = matrix(NA, nrow = 26, ncol = 26)
  colnames(A.init) = cantons
  rownames(A.init) = cantons
  
  # GE
  A.init[1,2] = 1
  
  # VD
  A.init[2,c(1,4,3,5,6)] = 1
  
  # FR
  A.init[3, c(2,6, 4)] = 1
  
  # NE
  A.init[4, c(7,6,2,3)] = 1
  
  # VS
  A.init[5, c(2,6,15,16)] = 1
  
  # BE
  A.init[6,c(3,2,5,15,13,14,12,8,7,4)] = 1
  
  # JU
  A.init[7,c(9,10,8,6,4)] = 1
  
  # SO
  A.init[8,c(9,11,6,7)] = 1
  
  # BL
  A.init[9,c(8,11,10,7)] = 1
  
  # BS
  A.init[10,c(7,9)] = 1
  
  # AG
  A.init[11,c(9,8,12, 17,21)] = 1
  
  # LU
  A.init[12,c(11,17,18,14,13,6)] = 1
  
  # OW
  A.init[13,c(12,6,14)] = 1
  
  # NW
  A.init[14,c(15,13,12,18,6)] = 1
  
  # UR
  A.init[15, c(14,6,5,16,25,19,18)] = 1
  
  # TI
  A.init[16, c(5,15, 25)] = 1
  
  # ZG
  A.init[17, c(18,21,11,12)] = 1
  
  # SZ
  A.init[18,c(17,14,15,19,20,21,12)] = 1
  
  # GL
  A.init[19,c(18,15,25,20)] = 1
  
  # SG
  A.init[20, c(25,24,23,26,21,18,19)] = 1
  
  # ZH
  A.init[21,c(20,26,22,11,17,18)] = 1
  
  # SH
  A.init[22,c(21,26)] = 1
  
  # AR
  A.init[23, c(20,24)] = 1
  
  # AI
  A.init[24, c(20,23)] = 1
  
  # GR
  A.init[25,c(20,19,15,16)] = 1
  
  # TG
  A.init[26,c(22,21,20)] = 1
  
  # NA to 0
  A.init[is.na(A.init)] = 0
  
  # Return adjacency matrix
  return(A.init)
}

sd.Informal = function(adjacency.matrix, N.cant, y, X, beta.hat) {
  
  # adjacency.matrix: adjacency matrix as computed in create.A.init
  # N.cant: number of units (26 cantons)
  # y: response vector with dim(y) = NT*1
  # X: design matrix with dim(X) = NT*p
  # beta.hat: vector of estimated coefficients with dim(beta.hat) = p*1 
  
  # Output: covariance matrix of estimated coefficients for informal method
  
  # Time dimension
  T.time = nrow(X)/N.cant
  
  # First index is for time, second index is for canton, third for variable
  X.array = array(X, dim = c(T.time,N.cant,ncol(X))) 
  
  # Residuals
  U = y - X%*%beta.hat 
  
  # Array of residuals along unit and time
  U.array = array(U,dim = c(T.time, N.cant, 1))
  
  # We firstly estimate the coefficient of time dependence as suggested by Hansen (2022)
  # through fitting an AR(1) on the time-sums St. To to this, we compute St and run the AR(1) regression
  # to obtain rho.hat. We do not use St going forward.
  
  # S component from Hansen (2022)
  St = sapply(1:T.time, function(t){
    colSums(X.array[t,,]*U.array[t,,1])
  }) 
  
  # Compute rho.hat
  rho = unlist(lapply(1:ncol(X), function(j){
    
    Y.rho = St[j,2:ncol(St)]
    X.rho = St[j,1:(ncol(St)-1)]
    
    fit.rho = lm(Y.rho ~-1+ X.rho)
    rho = coef(fit.rho)
    
  }))
  
  # Average over variables
  rho.mean = mean(rho)
  
  # Give weight of 0.5 to neighboring cantons as defined in the paper
  lambda = replace(adjacency.matrix, adjacency.matrix == 1, 0.5)
  
  # Make the diagonal values of lambda 1
  diag(lambda) = 1
  
  # We now construct the large weight matrix S with dim(S) = NT*NT. Note that S is a block-matrix with
  # three types of sub-matrices. On the diagonal of the block-matrix is the submatrix of type W.diag. If cantons i,j 
  # are neighbors, we add the sub-matrix W.neighbor to S. If cantons i,j are not neighbors and i is not j, we add a zero
  # sub-matrix. Note that all objects are accordingly ordered.
  
  # Construct the sub-matrix on the diagonal of S [lambda]_ij = 1
  W.diag = matrix(nrow = T.time, ncol = T.time)
  
  for (j in 1:T.time) {
    
    for (k in 1:T.time) {
      
      W.diag[j,k] = 0.5^abs(j-k) 
      
    } 
    
  }
  
  diag(W.diag) = 1
  
  # Construct the sub-matrix for neighboring cantons [lambda]_ij = 0.5
  W.neighbor = matrix(nrow = T.time, ncol = T.time)
  
  for (j in 1:T.time) {
    
    for (k in 1:T.time) {
      
      W.neighbor[j,k] = 0.5*0.5^abs(j-k)
      
    } 
    
  }
  
  # Construct the sub-matrix for cantons that are not neighbors
  W.null = matrix(data = 0, nrow = T.time, ncol = T.time)
  
  # Allocation of intermediate matrix M. M is a matrix with dim(M) = N.cant*N.cant. The elements of M
  # are to corresponding sub-matrices. They are however not unpacked into the large weight matrix S. That is done
  # under S below.
  
  # Allocation
  M = matrix(list(NA), ncol = ncol(lambda), nrow = nrow(lambda))
  
  # Construct the big matrix M
  for(p in 1:nrow(lambda)) {
    
    for (i in 1:ncol(lambda)) {
      
      if (lambda[p,i] == 0.5) {
        
        M[[p,i]] = W.neighbor
        
      } else if (lambda[p,i] == 1) {
        
        M[[p,i]] = W.diag
        
      } else {
        
        M[[p,i]] = W.null
        
      } #else
      
    } # Columns
    
  } # Rows
  
  # Unpack the elements of this matrix and construct the large weight matrix S 
  S = matrix(NA, ncol = ncol(lambda)*T.time, nrow = nrow(lambda)*T.time)
  
  for (s in 1:nrow(M)) {
    
    for (f in 1:ncol(M)) {
      
      S[(T.time*(s-1) + 1):(s*T.time), (T.time*(f-1) + 1):(f*T.time)] = M[[s,f]]
      
    } # Rows
    
  } # Columns
  
  # Residuals
  U = y - X%*%beta.hat
  
  # Meat of informal estimator
  Omega.informal = t(X)%*%(S*(U%*%t(U)))%*%X
  
  # Covariance matrix of informal estimator
  Var.Informal = solve(t(X)%*%X)%*%Omega.informal%*%solve(t(X)%*%X)
  
  # Return the informal covariance matrix
  return(list(Var.Informal = Var.Informal,
              rho.mean = rho.mean))
} 

compute.se = function(fit, model, frequency, response, infovar, type.effect) {
  
  # fit: = fitted model (plm object)
  # model: in {within, random, pooling}
  # frequency: in {daily, weekly}
  # response: in {median_R_mean}
  # model: in {within, random, pooling}
  # infovar: additional information variables
  # type.effect: direct or total to include behavior variables or not
  
  # Output: covariance matrix of estimated coefficients for all the different estimators
  
  # Standard heteroskedasticity consistent covariance matrix 
  hc.standard = vcovHC(fit, method = "white1", type = "HC3")
  
  # Clustered by canton estimator
  hc.group = vcovG(fit, cluster = "group", inner = "cluster", l = 0)
  
  # Clustered by time (day or one week) estimator
  hc.time = vcovG(fit, cluster = "time", inner = "cluster", l = 0)
  
  # Double clustered estimator (canton and time)
  hc.double = vcovDC(fit)
  
  # Newey-West panel estimator
  hc.neweywest = vcovNW(fit)
  
  # Method white1 is needed to match the lm covariance matrix estimator vcovHC. Note that to match the results using lm,
  # you have to use two-way demeaned data and a regression specification without an intercept.
  
  # The self-coded estimators (Hansen and informal) need (y, X) which need to be two-way demeaned for the 
  # fixed effects model. In order to keep the neat structure of code where we only need to input the fit and the
  # model into this function, we "reverse-engineer" which response and which frequency the fitted model is on. We then generate
  # the two-way demeaned data accordingly. we do not do that for RE separately as we only use standard (HC3) estimator for RE.
  
  # Save data from fit
  data.from.fit = fit[["model"]]
  
  # Fixed effects
  if (model == "within") {
    
    # Generate two-way demeaned y,X
    data = tw.demean(data.prep(lag = 7, shift = 14, response = response, frequency = frequency, r.infovar = 21)$data,
                     response = response, frequency = frequency)
    
    # Generate formula for the model 
    formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
    
    # Select y, X and save beta.hat from the fitted model
    y = as.vector(data$Y)
    X = as.matrix(data[,labels(terms(formula))])
    beta.hat = as.numeric(fit[["coefficients"]])
    
    # Covariance matrix Hansen
    hc.hansen = sd.Thompson.Hansen(y = y, X = X, N.cant = 26, beta.hat = beta.hat, M.T = 2)$Var.Hansen
    
    # Covariance matrix informal
    hc.informal = sd.Informal(adjacency.matrix = create.A.init(), N.cant = 26, beta.hat = beta.hat, y = y, X = X)$Var.Informal
    
  } # within
  
  
  # Only return HAC estimator for RE model but 7 times to ensure conformity with functions as the dimension of 
  # the output of the RE model is the same as for the one of the FE model when doing so.
  if (model == "random") {
    
    return(list(
      hc.standard = hc.standard,
      hc.standard = hc.standard,
      hc.standard = hc.standard,
      hc.standard = hc.standard,
      hc.standard = hc.standard,
      hc.standard = hc.standard,
      hc.standard = hc.standard))
    
  } else {
    
    return(list(
      hc.standard = hc.standard,
      hc.group = hc.group,
      hc.time = hc.time,
      hc.double = hc.double,
      hc.neweywest = hc.neweywest,
      hc.hansen = hc.hansen,
      hc.informal = hc.informal))
    
  } # Else
  
} # Compute.se

estimation = function(frequency, response, model, infovar, type.effect) {
  
  # frequency: in {daily, weekly}
  # response: response variable in {median_R_mean, casegrowth}
  # model: in {pooling, within, random}
  # infovar: include additional information variables for case-growth approach
  #type.effect: direct or total to include behavior variables or not
  
  # Output: list of (i) results (regression results depending on covariance matrix estimator)
  #                 (ii) initial fitted model
  #                 (iii) inputted data
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = model, infovar = infovar, type.effect = type.effect)
  
  # Store data; change configuration of data processing here to examine lag, shift and r.infovar
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  
  # Fit model 
  if (frequency == "daily") {
    
    fit = plm(formula, 
              data = data, model = model, random.method = "walhus",
              index = c("X.Canton_3", "X.day"), 
              effect = "twoways")
  } else {
    
    fit = plm(formula, 
              data = data, model = model, random.method = "walhus",
              index = c("X.Canton_3", "X.oneweek"), 
              effect = "twoways")
  }
  
  # Refit the results with the robust covariance matrix estimates
  cov.matrices = compute.se(fit = fit, model = model, frequency = frequency, response = response, 
                            infovar = infovar, type.effect = type.effect)
  results = vector("list", 7)
  
  for (j in 1:length(results)) {
    results[[j]] = coeftest(fit, vcov = cov.matrices[[j]])
  }
  
  # Return the results, the fitted model and the data
  return(list(results = results,
              fit = fit,
              data = data))
} # Estimation

multiple_split = function(response, frequency, infovar, type.effect) {
  
  # response: response variable in {median_R_mean, casegrowth}
  # frequency: frequency in {daily, weekly}
  # infovar: add additional information variables
  # type.effect: direct or total to include behavior variables or not
  
  # Output: bias.corrected estimate of mask policy
  
  # Data
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = "within", infovar = infovar, type.effect = type.effect)
  
  # Number of splits along the cross-section
  s = 20
  
  # Fit non-debiased model daily or weekly and save time and unit for sampling later
  if (frequency == "daily") {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.day"), effect = "twoways", random.method = "walhus")
    time   = as.double(data$X.day)
    unit   = as.double(data$X.Canton_3)
    
  } else {
    
    uc     = plm(formula, data = data, model = "within", index = c("X.Canton_3", "X.oneweek"), effect = "twoways", random.method = "walhus")
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
                   index = c("X.Canton_3", "X.day"), effect = "twoways", random.method = "nerlove")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.day"),effect = "twoways", random.method = "nerlove")
      
      
    } else {
      
      cross1 = plm(formula, data = data[subsample1,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"), effect = "twoways", random.method = "nerlove")
      
      cross2 = plm(formula, data = data[subsample2,], model = "within",
                   index = c("X.Canton_3", "X.oneweek"),effect = "twoways", random.method = "nerlove")
      
    }
    across = across + ((coef(cross1) + coef(cross2))/2)/s
  }
  
  # Average cross over corrected
  acbc = 2 * coef(uc) - across
  return(acbc)
}

data_wb = function(data, mle) {
  
  # data: data used
  # mle: mandatory second argument to ran.gen
  
  # Output: data with added sweight for cantonal clustered multiplier bootstrap. We use exponential weights.
  
  # Two dimensions of the panel data
  if (frequency == "daily") {
    
    n = length(unique(data$X.Canton_3))
    t = length(unique(data$X.day))
    
  } else {
    
    n = length(unique(data$X.Canton_3))
    t = length(unique(data$X.oneweek))
  }
  
  # Exponentially distributed weights (multipliers)
  multipliers = rexp(n)
  
  # For each state, weight is the same for all dates
  weight = rep(multipliers/sum(multipliers), each = t)
  
  # Add to the data and treat it as sampling weight
  data$sweight = weight
  return(data)
}

bootstat_fe = function(data) {
  
  # data: data used
  
  # Output: estimated coefficients; this function is the statistic used in the bootstrap
  
  # Fix the model to within
  formula = construct.formula(response = response, frequency = frequency, model = "within", infovar = infovar, type.effect = type.effect)
  
  if (frequency == "daily") {
    
    fit = plm(formula,
              data = data, model = "within",
              weights = sweight,
              index = c("X.Canton_3", "X.day"),
              effect = "twoways")
  } else {
    
    fit = plm(formula,
              data = data, model = "within",
              weights = sweight,
              index = c("X.Canton_3", "X.oneweek"),
              effect = "twoways")
  }
  
  # Store coefficients
  tot = fit$coefficients
  
  # Return coefficients
  return(unlist(tot))
}

bs.estimation = function(response, frequency, infovar, type.effect) {
  
  # response: response variable
  # frequency: daily or weekly
  # infovar: additional information variables added
  # type.effect: direct or total to include behavior variables or not
  
  # Output: point estimate of W and bootstrapped standard errors
  
  # Load data and add variable sweight as initialization
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  data$sweight = 1
  
  # Cores needed (run sufficiently quick on just 1) 
  ncores = 1
  
  # Run bootstrap (500 replications)
  result_boot = boot(data = data, statistic = bootstat_fe,
                     sim = "parametric", ran.gen = data_wb, mle = 0,
                     parallel = "multicore", ncpus = ncores, R = 1000)
  
  # Obtain the matrix of bootstrap estimates
  result = structure(vapply(result_boot$t, as.double, numeric(1)), 
                     dim = dim(result_boot$t))
  
  # Inference for W via basic bootstrap
  ci = boot.ci(result_boot, index = 8, conf = 0.95, 
               type = "basic")
  
  # Return basic confidence interval if desired. Un-comment the following code if so,
  # which produces bootstrap standard errors to then construct confidence intervals via 
  # the standard normal approximation. 
  # return(ci)
  
  # Inference for W via standard errors and normal approximation
  bse.standard = apply(result, 2, function(x) {
    
    # Standard error
    return(sd(x))
  })
  
  # Inference for W via robust errors and normal approximation
  bse.robust = apply(result, 2, function(x) {
    
    # Standard error
    return(mad(x))
  })
  
  # Combine point estimate and standard errors
  tab = rbind(result_boot$t0, bse.standard, bse.robust)
  
  # Return table
  return(tab)
}

dml.estim = function(response, frequency, infovar, type.effect) {
  
  # frequency: in {daily, weekly}
  # response: response variable
  # infovar: additional information variables 
  # type.effect: direct or total to include behavior variables or not
  
  # Output: list of B, point estimate of W, call and basic confidence interval for W
  
  # Set seed 
  set.seed(42)
  
  # Store data. Add weekly dummies to the double machine learning approach to capture time trends
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  data = dummy_cols(data, select_columns = "X.oneweek")
  
  # Store formula and extract y, X, W
  formula = construct.formula(response = response, frequency = frequency, model = "pooling", infovar = infovar, type.effect = type.effect)
  x_cols = c(labels(terms(formula)), paste("X.oneweek", sep = "_", min(data$X.oneweek):max(data$X.oneweek)))
  x_cols = x_cols[!x_cols %in% c("W")]
  y_col  = "Y"
  d_col  = "W"
  
  if (frequency == "daily") {
    
    dml_data = DoubleMLClusterData$new(data = data,
                                       y_col  = y_col,
                                       d_cols = d_col,
                                       z_cols = NULL,
                                       cluster_cols = c("X.Canton_3", "X.day"), 
                                       x_cols = x_cols)
  } else {
    
    dml_data = DoubleMLClusterData$new(data = data,
                                       y_col  = y_col,
                                       d_cols = d_col,
                                       z_cols = NULL,
                                       cluster_cols = c("X.Canton_3", "X.oneweek"), 
                                       x_cols = x_cols)
  }
  
  # Initialize ML learner; do not use classifiers to ensure that the overlap assumption holds
  lasso  = lrn("regr.cv_glmnet", nfolds = 10, s = "lambda.min")
  lasso_class = lrn("classif.cv_glmnet", nfolds = 10, s = "lambda.min")
  forest = lrn("regr.ranger", num.trees = 300)
  forest_class = lrn("classif.ranger", num.trees = 300)
  
  # Double Machine Learning
  dml_res = DoubleMLPLR$new(
    data = dml_data,
    ml_g = forest, # Outcome model
    ml_m = forest, # Treatment model 
    n_folds = 5,
    n_rep = 1,
    score = "partialling out",
    dml_procedure = "dml1",
    draw_sample_splitting = TRUE,
    apply_cross_fitting = TRUE
  )
  
  # Fit
  dml_res$fit(store_predictions = TRUE)
  dml_res$coef
  dml_res$se
  
  # Nuisance parameters
  #pred.nuisance = dml_res$predictions
  #hist(pred.nuisance[["ml_g"]])
  #hist(pred.nuisance[["ml_m"]], breaks = 60)
  
  # Return the results 
  return(dml_res$summary())
  
}

r.sq = function(fit, frequency, response) {
  
  # fit: fit of plm model
  # frequency: in {daily, weekly}
  # response: response variable
  
  # Output: R^2 for plm model. We need the original data and not the two-way demeaned data
  # so we store them using data.prep
  
  # Store data
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data
  
  # SST and SSE
  sst = sum((data$Y - mean(data$Y))^2)
  sse = t(residuals(fit))%*%residuals(fit)
  
  # Return R^2
  return((sst - sse)/sst)
}

diagnostics = function(fit, frequency, response, infovar, type.effect) {
  
  # fit: fit of plm model
  # frequency: in {daily, weekly}
  # response: response variable
  # infovar: additional information variables
  # type.effect: direct or total to include behavior variables or not
  
  # Output: estimate of W from fixed effects approach of with full data and data restricted on cooks distance
  
  # Normality of residuals
  qqnorm(residuals(fit), ylab = "Residuals")
  qqline(residuals(fit))
  
  # Fitted vs. res plot
  preds = as.numeric(fit[["model"]][["Y"]]) - as.numeric(fit$residuals)
  res = as.numeric(fit$residuals)
  plot(preds, res, ylab = "Residuals", xlab = "Predictions")
  
  # Construct the fitted vs. residuals plot for each canton separately. In this fashion, it 
  # is more like an acf plot for time series objects.
  
  # # Number of observations per canton
  # time.ind = length(fit[["model"]][["Y"]])/26
  # 
  # for (i in 1:26) {
  #   
  #   par(mfrow = c(4,5))
  #   plot(res[((i-1)*time.ind + 1):(i*time.ind)], xlab = "Time", ylab = "Residuals")
  #   
  # } 
  # 
  # # Reset margins
  # par(mfrow = c(1,1))
  
  # Cooks distance; there is no option to compute the cooks distance in the plm package
  # so it`s easiest to run the fixed effects model via lm (which yields the exact same model) and use the designated 
  # function provided in that package. We do it for the fixed effects models.
  
  # Two-way demean data to use in lm 
  data = tw.demean(data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = frequency)$data,
                   response = response, frequency = frequency)
  
  # Store formula
  formula = construct.formula(response = response, frequency = frequency, model = "within", infovar = infovar, type.effect = type.effect)
  
  # Run fixed effects regression via lm
  fit.original = lm(formula, data = data)
  
  # Save orginal R^2
  r.squared.original = r.sq(fit = fit.original, frequency = frequency, response = response)
  
  # Compute cooks distance and determine influential observations
  cooksD = cooks.distance(fit.original)
  influential = cooksD[(cooksD > 4/nrow(data))]
  influential.05 = cooksD[(cooksD > 0.5)]
  #influential = cooksD[(cooksD > (4 * mean(cooksD, na.rm = TRUE)))]
  
  # Remove outlaying observations
  names_of_influential = names(influential)
  outliers = data[names_of_influential,]
  data_without_outliers = data %>% 
    anti_join(outliers)
  
  # Run fixed effects regression via lm on reduced data
  fit.reduced = lm(formula, data = data_without_outliers) 
  
  # Save new R^2
  r.squared.reduced = r.sq(fit = fit.reduced, frequency = frequency, response = response)
  
  # Store point estimate and std err as well as the difference of point estimates
  fit.entire     = summary(fit.original)$coefficients["W",1]
  std.err.entire = summary(fit.original)$coefficients["W",2]
  fit.small      = summary(fit.reduced)$coefficients["W",1]
  std.err.small  = summary(fit.reduced)$coefficients["W",2]
  diff           = fit.entire - fit.small
  
  # Return
  return(list(fit.entire = fit.entire,
              std.err.entire = std.err.entire,
              fit.small = fit.small,
              std.err.small = std.err.small,
              diff = diff,
              data = data,
              influential = influential,
              influential.05 = influential.05,
              preds = preds,
              res = res,
              cooksD = cooksD,
              r.squared.original = r.squared.original,
              r.squared.reduced = r.squared.reduced))
} 


data.prep.halfcantons = function(shift, response, startdate = "2020-07-06",
                                 enddate = "2020-12-21", r.infovar) {
  
  # shift: shifting of response when using casegrowth as defined 
  # response: in {casegrowth, median_R_mean}
  # startdate, enddate: to be adjusted for different time periods
  # r.infovar: lag length of information variables
  
  # Output: data frame used for the estimation (Y, W, X) while merging the half cantons
  
  # Note that the baseline mask policy of mandatory mask wearing was introduced nationally on "2020-07-06" and the 
  # first vaccinations took place on "2020-12-21". This is the reason for our sample selection.
  
  # Read in data
  DataCovid = read.csv(".\\Data\\DataCovid.csv") 
  dat.new = read.csv(".\\Data\\policy_stuff.csv", header = TRUE, sep = ",")
  
  # Merge data on common days
  data = merge(DataCovid, dat.new)
  
  # Delete unnecessary variables to might lead to conflicts
  data = data %>%
    dplyr::select(-c(deaths_1d, recovered_1d, tests_1d, hosp_1d, value.total_1d, sre000d0_1d, tre200d0_1d, ure200d0_1d, O65P, restGatheringsCH)) %>%
    rename(testingPolicy = testing_policyCH)
  
  # Transform percentage to be a real percentage 
  data = data %>%
    mutate(percentage_age = percentage_age*100)
  
  # Get the weekly date as integer. Further construct disjoint weeks
  data$oneweek = cut.Date(as.Date(data$datum), breaks = "1 week", start.on.monday = TRUE, labels = FALSE) 
  
  # Extract which value of oneweek corresponds to startdate, enddate to apply the time selection
  # via the vector oneweek later
  start.day = which(data$datum == startdate)[1]
  start.week = data[start.day, "oneweek"]
  end.day   = which(data$datum == enddate)[1]
  end.week = data[end.day, "oneweek"]
  
  # Make an indicator that treats normal cantons as separate cantons but half cantons as such (same value for the indicator)
  data$geoRegion.2 = ifelse(data$geoRegion %in% c("BS", "BL"), "BSBL",
                            ifelse(data$geoRegion %in% c("AI", "AR"), "AIAR",
                                   ifelse(data$geoRegion %in% c("OW", "NW"), "OWNW", data$geoRegion)))
  
  if (response == "casegrowth") {
    
    # Construct target variable and lagged target variable. We do that by computing the new cases per canton and week. We can then
    # compute the log-growth rate using that series, which constructs 7 times that same growth rate per canton and week. When aggregating, 
    # we can then take the mean over disjoint weeks to obtain the desired variable. We do that for all growth rates.
    
    # Construct national aggregates
    sumTotalCH = aggregate(data$sumTotal, by = list(data$datum), FUN  = sum)
    sumTotalCH = rename(sumTotalCH, datum = Group.1)
    sumTotalCH = rename(sumTotalCH, sumTotalCH = x)
    data = merge(sumTotalCH, data)
    data = data %>%
      arrange(geoRegion)
    
    # Get new cases per day in each canton. Do it for the new cases on a national level as well.
    data = data %>%
      group_by(geoRegion) %>%
      mutate(NewCases = c(0,diff(sumTotal)),
             NewCases.national = c(0, diff(sumTotalCH))) %>%
      ungroup()
    
    data = data %>%
      relocate(sumTotal, .before = NewCases)
    
    # We now have all variables constant across oneweek and canton which need to be constant. The other ones
    # are aggregated to weekly data by simply averaging the values per week and canton. Note that we take the sum
    # for transactions and tests as we need the new number of transactions and tests per canton and week. These variables 
    # already come supplied as new transactions and tests so we do not have to treat them as the cases where we had to explicitly 
    # calculate the new cases as done in the above code.
    data = data %>%
      group_by(geoRegion.2, oneweek) %>%
      summarize(newcases_we           = sum(NewCases),
                newcases.national_we  = sum(NewCases.national),
                newtransactions_we    = sum(Number.of.transactions),
                newtests_we           = sum(tests),
                percentage_age_we = mean(percentage_age),
                Density_we        = mean(Density),
                population_we     = mean(population),
                grocery_and_pharmacy_we   = mean(grocery_and_pharmacy),
                transit_stations_we       = mean(transit_stations),
                workplaces_we             = mean(workplaces),
                ferien_we                 = mean(ferien), 
                sre000d0_we               = mean(sre000d0),
                tre200d0_we               = mean(tre200d0),
                ure200d0_we               = mean(ure200d0),
                schoolClosing_we          = mean(schoolClosing), 
                restGatherings_we         = mean(restGatherings), 
                cancEvents_we             = mean(cancEvents), 
                testingPolicy_we          = mean(testingPolicy), 
                workClosing2a_we          = mean(workClosing2a), 
                mask_treat_we             = mean(mask_treat)) %>%
      ungroup()
    
    # Using newcases_we, newcases.national_we, newtransactions_we and newtests_we, we can now compute the weekly growth rates.
    
    # Create the response variable and lead response variable. 
    # Check if values are below the threshold of 0.1 and set it to 0.1 if so. Does not happen
    data$newcases_we[data$newcases_we < 0.1] = 0.1
    data$casegrowth_we = paneldiff(x = log(data$newcases_we), id = data$geoRegion.2, t = as.Date(data$oneweek), l = 1)
    data$casegrowth_we = panellag(x = data$casegrowth_we, i = data$geoRegion.2, t = as.Date(data$oneweek), lag = -shift/7)
    
    # Lagged target variable
    data$casegrowthlag_we =  panellag(x = data$casegrowth_we, i = data$geoRegion.2, t = data$oneweek, lag = 2)
    data$casegrowthlag1_we = panellag(x = data$casegrowth_we, i = data$geoRegion.2, t = data$oneweek, lag = 1)
    
    # log Delta Cit
    data$logdiffcases_we = lognozero(paneldiff(x = data$newcases_we, id = data$geoRegion.2, t = as.Date(data$oneweek), lag = 1))
    
    # Lagged target variable on national level; works with id argument
    data$newcases.national_we[data$newcases.national_we < 0.1] = 0.1
    data$casegrowthlag.national_we = paneldiff(x = log(data$newcases.national_we), id = data$geoRegion.2, t = as.Date(data$oneweek), l = 1)
    
    # log Delta Cit; same as above but national level
    data$logdiffcases.national_we = lognozero(paneldiff(x = data$newcases.national_we, id = data$geoRegion.2, t = as.Date(data$oneweek), lag = 1))
    
    # Delta log Tit; already new tests so we do not have to apply the difference operator
    data$newtests_we[data$newtests_we < 0.1] = 0.1
    data$difflogtests_we = paneldiff(x = log(data$newtests_we), id = data$geoRegion.2, t = as.Date(data$oneweek), l = 1)
    
    # Delta log transactions 
    data$newtransactions_we[data$newtransactions_we < 0.1] = 0.1
    data$transactiongrowth_we = paneldiff(x = log(data$newtransactions_we), id = data$geoRegion.2, t = as.Date(data$oneweek), l = 1)
    
    # Correct period and name
    data = data %>%
      filter(oneweek >= start.week & oneweek <= end.week)
    data = rename(data, Canton_3 = geoRegion.2)
    
    
    # Construct X weekly
    X = data[,c("casegrowthlag_we", "casegrowthlag1_we", "logdiffcases_we", "casegrowthlag.national_we", "logdiffcases.national_we", "difflogtests_we",
                "percentage_age_we", "Density_we", "population_we",
                "schoolClosing_we", "restGatherings_we", "cancEvents_we",
                "testingPolicy_we", "workClosing2a_we",
                "grocery_and_pharmacy_we", "transit_stations_we", "transactiongrowth_we", "workplaces_we",
                "sre000d0_we", "tre200d0_we", "ure200d0_we", "ferien_we", "Canton_3", "oneweek")]
    W = data$mask_treat_we
    Y = data$casegrowth_we
    
  }
  
  if (response == "median_R_mean") {
    
    # Take the mean over the disjoint weeks for (Y,W,X)
    data = data %>%
      group_by(geoRegion.2, oneweek) %>%
      summarize(median_R_mean_we      = mean(median_R_mean),
                percentage_age_we = mean(percentage_age),
                Density_we        = mean(Density),
                population_we     = mean(population),
                grocery_and_pharmacy_we   = mean(grocery_and_pharmacy),
                transit_stations_we       = mean(transit_stations),
                newtransactions_we        = sum(Number.of.transactions),
                workplaces_we             = mean(workplaces),
                ferien_we                 = mean(ferien), 
                sre000d0_we               = mean(sre000d0),
                tre200d0_we               = mean(tre200d0),
                ure200d0_we               = mean(ure200d0),
                schoolClosing_we          = mean(schoolClosing), 
                restGatherings_we         = mean(restGatherings), 
                cancEvents_we             = mean(cancEvents), 
                testingPolicy_we          = mean(testingPolicy),
                workClosing2a_we          = mean(workClosing2a), 
                mask_treat_we             = mean(mask_treat)) %>%
      ungroup()
    
    # Delta log transactions 
    data$newtransactions_we[data$newtransactions_we < 0.1] = 0.1
    data$transactiongrowth_we = paneldiff(x = log(data$newtransactions_we), id = data$geoRegion.2, t = as.Date(data$oneweek), l = 1)
    
    # Lagged target variable
    data$median_R_mean.lag_we = panellag(x = data$median_R_mean_we, i = data$geoRegion.2, t = data$oneweek, lag = r.infovar/7)
    data$median_R_mean.lag1_we = panellag(x = data$median_R_mean_we, i = data$geoRegion.2, t = data$oneweek, lag = 1)
    
    # Correct period and name
    data = data %>%
      filter(oneweek >= start.week & oneweek <= end.week)
    data = rename(data, Canton_3 = geoRegion.2)
    
    # Construct X weekly
    X = data[,c("median_R_mean.lag_we", "median_R_mean.lag1_we",
                "percentage_age_we", "Density_we", "population_we",
                "schoolClosing_we", "restGatherings_we", "cancEvents_we",
                "testingPolicy_we", "workClosing2a_we",
                "grocery_and_pharmacy_we", "transit_stations_we", "transactiongrowth_we", "workplaces_we",
                "sre000d0_we", "tre200d0_we", "ure200d0_we", "ferien_we", "Canton_3", "oneweek")]
    
    W = data$mask_treat_we
    Y = data$median_R_mean_we
    
  }
  
  # Return  
  return(list(data = data.frame(X = X, Y = Y, W = W)))
  
} # Data preparation half cantons

################################################################################




