# (I) OVERVIEW

# This script runs all the proposed models using the helper function that are sourced
# from helperfunctions.R.

################################################################################

# (II) ENVIRONMENT AND PACKAGES

# Empty environment
rm(list = ls())


set.seed(42) 
source("./Code/helperfunctions.R")

################################################################################

# (III) ESTIMATION

# Start with the fixed effects and random effects models
freq.poss = c("weekly")
response.poss = c("median_R_mean", "casegrowth")
model.poss = c("within", "random")
type.effect.poss = c("direct", "total")

# Allocation of results
results.list.fe = vector("list", 8)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    for (model in model.poss) {
      
      for (type.effect in type.effect.poss) {
        
        results.list.fe[[i]] = list(results = estimation(frequency = frequency, response = response, model = model, infovar = FALSE, type.effect = type.effect,
                                                         lag.one = FALSE),
                                    frequency = frequency,
                                    response = response,
                                    model = model,
                                    type.effect = type.effect)
        i = i+1
      }
    }
  }
}



# Additional information variables
results.list.fe.infovar = vector("list", 2)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (type.effect in type.effect.poss) {
    
    results.list.fe.infovar[[i]] = list(estimation(frequency = frequency, response = "casegrowth", model = "random", infovar = TRUE, type.effect = type.effect,
                                                   lag.one = FALSE),
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
      
      results.list.dfe[[i]] = list(multiple_split(response = response, frequency = frequency, infovar = FALSE, type.effect = type.effect,
                                                  lag.one = FALSE),
                                   frequency = frequency,
                                   response = response,
                                   model = "debiased",
                                   type.effect = type.effect)
      i = i+1
    }
  }
}



# Double machine learning
results.list.dml = vector("list", 8)
i = 1

# Run the models
for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    for (infovar in c(TRUE, FALSE)) {
      
      for (type.effect in type.effect.poss) {
        
        results.list.dml[[i]] = list(dml.estim(response = response, frequency = frequency, infovar = infovar, type.effect = type.effect,
                                               lag.one = FALSE),
                                     response = response,
                                     frequency = frequency,
                                     infovar = infovar,
                                     type.effect = type.effect)
        i = i+1
      }
    }
  }
}