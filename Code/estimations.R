# (I) OVERVIEW

# This script runs all the proposed models using the helper function that are sourced
# from helperfunctions.R. Finally, the results are compiled such that they can be further 
# called to create tables and plots.

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

# (III) ESTIMATION

# We will generate all models of interest. We start with the 2^3 models that we run
# through the estimation.R function, which are the fixed effects and random effects models. For case-growth,
# we further look at the results when including the additional information variables.
# We then run the 2^2 models that we run through the 
# multiple_split.R function, which are the de-biased fixed effects models. Lastly, we estimate
# the 2^3 models using the double machine learning framework. We keep frequency in the loop for easy
# access to the results of the daily models if wanted.

# Start with the fixed effects and random effects models
freq.poss = c("weekly")
response.poss = c("median_R_mean", "casegrowth")
model.poss = c("within", "random")
#model.poss = c("within")
type.effect.poss = c("direct", "total")

# Allocation of results
results.list.fe = vector("list", 8)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    for (model in model.poss) {
      
      for (type.effect in type.effect.poss) {
        
        results.list.fe[[i]] = list(results = estimation(frequency = frequency, response = response, model = model, infovar = FALSE, type.effect = type.effect),
                                    frequency = frequency,
                                    response = response,
                                    model = model,
                                    type.effect = type.effect)
        i = i+1
      }
    }
  }
}

# Run the models for response = case-growth and additional information variables. Note that we can only
# run it using the random effects approach as some of the additional information variables are at a national
# level which is not compatible with the fixed effects approach due to time fixed effects.
results.list.fe.infovar = vector("list", 2)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (type.effect in type.effect.poss) {
    
    results.list.fe.infovar[[i]] = list(estimation(frequency = frequency, response = "casegrowth", model = "random", infovar = TRUE, type.effect = type.effect),
                                        frequency = frequency,
                                        response = "casegrowth",
                                        model = "random",
                                        type.effect = type.effect)
    i = i+1
  }
}

# Run the de-biased models via multiple sample splitting
results.list.dfe = vector("list", 4)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    for(type.effect in type.effect.poss) {
      
      results.list.dfe[[i]] = list(multiple_split(response = response, frequency = frequency, infovar = FALSE, type.effect = type.effect),
                                   frequency = frequency,
                                   response = response,
                                   model = "debiased",
                                   type.effect = type.effect)
      i = i+1
    }
  }
}

# Additional, we run the double machine learning models using the lasso to estimate the nuisance functions 
# and cluster-robust standard errors, which are clustered in two-way manner, meaning on time and cantons. For time,
# we use day for the daily model and weeks for the weekly model, just as in the fixed effects models.
results.list.dml = vector("list", 8)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    for (infovar in c(TRUE, FALSE)) {
      
      for (type.effect in type.effect.poss) {
        
        results.list.dml[[i]] = list(dml.estim(response = response, frequency = frequency, infovar = infovar, type.effect = type.effect),
                                     response = response,
                                     frequency = frequency,
                                     infovar = infovar,
                                     type.effect = type.effect)
        i = i+1
      }
    }
  }
}

# Finally, we estimate cantonal-clustered multiplier bootstrapped standard errors for the fixed effects
# models.
results.list.bootstrap = vector("list", 4)
i = 1

# Run the bootstrap with B = 1000
# Run the models
for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    for (type.effect in type.effect.poss) {
      
      results.list.bootstrap[[i]] = list(bs.estimation(frequency = frequency, response = response, infovar = FALSE, type.effect = type.effect),
                                         frequency = frequency,
                                         response = response,
                                         type.effect = type.effect)
      i = i+1
    }
  }
}

################################################################################

# (IV) STORE THE RESULTS (POINT ESTIMATES, STD. ERRORS, P-VALUES)

# Extract the results from the 5 lists that are differently structured. We do this in the order
# that they are run in the above code. Run the code in order (see variable j)

# Empty matrix with the corresponding dimensions
results.dat = matrix(NA, nrow = 102, ncol = 7)
colnames(results.dat) = c("estimate", "std.error", "p.val", "model", "term", "type.effect", "add.infovar")

# Model abbreviations from response and model type
abbrev = function(response, model) {
  
  if (response == "median_R_mean" && model == "within") {
    return("FE R")
  } else if (response == "median_R_mean" && model == "random") {
    return("RE R")
  } else if (response == "casegrowth" && model == "within") {
    return("FE Casegrowth")
  } else if (response == "casegrowth" && model == "random") {
    return("RE Casegrowth")
  } else if (response == "median_R_mean" && model == "debiased") {
    return("DFE R")
  } else if (response == "casegrowth" && model == "debiased") {
    return("DFE Casegrowth")
  }
}

# Names of standard error estimators in correct ordering
cov.names = c("HC3", "Canton", "Time", "Canton-Time", "NW", "CH", "Own")

# Go through standard fixed effects and random effect estimation
j = 1
for (p in 1:8) {
  
  for (k in 1:7) {
    
    results.dat[j,1] = results.list.fe[[p]][["results"]][["results"]][[k]]["W", 1]
    results.dat[j,2] = results.list.fe[[p]][["results"]][["results"]][[k]]["W", 2]
    results.dat[j,3] = results.list.fe[[p]][["results"]][["results"]][[k]]["W", 4]
    results.dat[j,4] = cov.names[k]
    results.dat[j,5] = abbrev(response = results.list.fe[[p]][["response"]],
                              model    = results.list.fe[[p]][["model"]])
    results.dat[j,6] = results.list.fe[[p]][["type.effect"]]
    results.dat[j,7] = "FALSE"
    j = j + 1
  }
}

# Go through random effects models of case-growth with the additional information variables
# as in Chernozhukov (2020)
for (q in 1:2) {
  
  results.dat[j,1] = results.list.fe.infovar[[q]][[1]][["results"]][[1]]["W", 1]
  results.dat[j,2] = results.list.fe.infovar[[q]][[1]][["results"]][[1]]["W", 2]
  results.dat[j,3] = results.list.fe.infovar[[q]][[1]][["results"]][[1]]["W", 4]
  results.dat[j,4] = "HC3"
  results.dat[j,5] = abbrev(response = results.list.fe.infovar[[q]][["response"]],
                            model    = results.list.fe.infovar[[q]][["model"]])
  results.dat[j,6] = results.list.fe.infovar[[q]][["type.effect"]]
  results.dat[j,7] = "TRUE"
  j = j + 1
}

# Extract the results for the de-biased models. Match the specifications from the de-biased models
# to those of the non de-biased models.
for (d in 1:4) {
  
  for (i in 1:7) {
    
    results.dat[j,1] = results.list.dfe[[d]][[1]][["W"]]
    results.dat[j,4] = cov.names[i]
    results.dat[j,5] = abbrev(response = results.list.dfe[[d]][["response"]],
                              model    = results.list.dfe[[d]][["model"]])
    results.dat[j,6] = results.list.dfe[[d]][["type.effect"]]
    results.dat[j,7] = "FALSE"
    j = j + 1
  }
}

results.dat[(59:(59+6)),2] = results.dat[1:7,2]       # Std. errors from non de-biased model for FE R direct
results.dat[(66:(66+6)),2] = results.dat[8:(8+6),2]   # Std. errors from non de-biased model for FE R total
results.dat[(73:(73+6)),2] = results.dat[29:(29+6),2] # Std. errors from non de-biased model for FE Case-growth direct
results.dat[(80:(80+6)),2] = results.dat[36:(36+6),2] # Std. errors from non de-biased model for FE Case-growth direct

# Compute the p-values for the de-biased models via the t-distribution. Depending on the specification, we have a different
# number of degrees of freedom, i.e:
#      - FE Weekly direct: 589
#      - FE Weekly total: 590
#      - RE Weekly direct: 630
#      - RE Weekly total: 634

for (e in 59:86) {
  
  # How many degrees of freedom
  if (grepl("FE", results.dat[e,"term"]) == TRUE && results.dat[e,"type.effect"] == "direct") {
    df = 589
  } else if (grepl("FE", results.dat[e,"term"]) == TRUE && results.dat[e,"type.effect"] == "total") {
    df = 590
  } else if (grepl("RE", results.dat[e,"term"]) == TRUE && results.dat[e,"type.effect"] == "direct") {
    df = 630
  } else if (grepl("RE", results.dat[e,"term"]) == TRUE && results.dat[e,"type.effect"] == "total") {
    df = 634
  }
  
  # Compute p value
  results.dat[e,3] = 2*pt(-abs(as.numeric(results.dat[[e,1]])/as.numeric(results.dat[[e,2]])), df = df)
}

# Get the results from the double machine learning
for (l in 1:8) {
  
  results.dat[j,1] = results.list.dml[[l]][[1]][1]
  results.dat[j,2] = results.list.dml[[l]][[1]][2]
  results.dat[j,3] = results.list.dml[[l]][[1]][4]
  results.dat[j,4] = "DML"
  
  if (results.list.dml[[l]][["response"]] == "median_R_mean") {
    results.dat[j,5] = "DML R"
  } else {
    results.dat[j,5] = "DML Casegrowth"
  }
  
  results.dat[j,6] = results.list.dml[[l]][["type.effect"]]
  results.dat[j,7] = results.list.dml[[l]][["infovar"]]
  j = j + 1
}

# Get the results from the bootstrap
for (n in 1:4) {
  
  results.dat[j,1] = results.list.bootstrap[[n]][[1]][1,"W"]
  results.dat[j,2] = results.list.bootstrap[[n]][[1]][2,"W"]
  
  if (results.list.bootstrap[[n]][["type.effect"]] == "direct") {
    df = 589
  } else {
    df = 590
  }
  results.dat[j,3] = 2*pt(-abs(results.list.bootstrap[[n]][[1]][1,"W"]/results.list.bootstrap[[n]][[1]][2,"W"]), df = df)
  results.dat[j,4] = "Canton-Bootstrap"
  results.dat[j,5] = abbrev(response = results.list.bootstrap[[n]][["response"]], model = "within")
  results.dat[j,6] = results.list.bootstrap[[n]][["type.effect"]]
  results.dat[j,7] = "FALSE"
  j = j + 1
}

# Bootstrapped standard errors for de-biased models. Make an array containing a 
# triple in (row) x (row) x (degree of freedom) x (term) x (type.effect)-space to loop over.
v1 = c(59, 95, 589, "DFE R", "direct") # DFE R direct
v2 = c(66, 96, 590, "DFE R", "total")  # DFE R total
v3 = c(73, 97, 589, "DFE Casegrowth", "direct") # DFE case growth direct
v4 = c(80, 98, 590, "DFE Casegrowth", "total") # DFE case growth total
arr = array(c(v1, v2, v3, v4), dim = c(5,4))

for (g in 1:4) {
  
  results.dat[j,1] = as.numeric(results.dat[as.numeric(arr[1,g]),1])
  results.dat[j,2] = as.numeric(results.dat[as.numeric(arr[2,g]),2])
  results.dat[j,3] = 2*pt(-abs(as.numeric(results.dat[j,1])/as.numeric(results.dat[j,2])), df = as.numeric(arr[3,g]))
  results.dat[j,4] = "Canton-Bootstrap"
  results.dat[j,5] = arr[4,g]
  results.dat[j,6] = arr[5,g]
  results.dat[j,7] = FALSE
  j = j + 1
}

# Save as data frame and export it as csv file
results.dat = as.data.frame(results.dat)

# Remove rows with non-HC3 estimators for the random effects model as they are not actually estimated.
results.dat = results.dat %>%
  filter(!(model != "HC3" & term == "RE R")) %>%
  filter(!(model != "HC3" & term == "RE Casegrowth")) 

# Write csv
write.csv(results.dat,"./Data/results_short_period_weather.csv", row.names = FALSE)


