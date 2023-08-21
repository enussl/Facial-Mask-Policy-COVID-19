# (I) OVERVIEW

# This script produces the plots from the paper.

################################################################################

# (II) ENVIRONMENT AND PACKAGES

# Empty environment
rm(list = ls())

# Set working directory to the root node of folder structure
#setwd(".\\Mask_Project\\Final")
setwd("C:/Users/eminu/OneDrive/Desktop/Facial-Mask-Policy-COVID-19")

# Read helper functions. Note that they are in the same directory. Add the corresponding path
# otherwise.
#source(".\\Scripts\\helperfunctions.R")
source("./Code/helperfunctions.R")


# # Read in fonts
# font_install("fontcm")
# loadfonts()

################################################################################

# (III) PLOTS

# (i) CI Plots

# Read in data; change col-names to correspond to dw-plot; make ordering consistent; remove those with additional 
# information variables as well as the DML rows.
results = read.csv(".\\Data\\results_short_period_weather.csv", header = T, sep = ",", stringsAsFactors = FALSE)
results = results %>%
  filter(!model == "DML") %>%
  filter(!model == "Canton-Bootstrap") %>%
  filter(add.infovar == FALSE)

results[results == "FE R"] = "FE r"
results[results == "DFE R"] = "DFE r"
results[results == "RE R"] = "RE r"
results[results == "FE Casegrowth"] = "FE growth.new.cases"
results[results == "DFE Casegrowth"] = "DFE growth.new.cases"
results[results == "RE Casegrowth"] = "RE growth.new.cases"
results[results == "Canton-Time"] = "Canton-Week"
results[results == "Time"] = "Week"



results.direct = results %>%
  filter(type.effect == "direct")
results.total = results %>%
  filter(type.effect == "total")



# Direct effect
pdf(".\\Plots\\ci_plot_final_direct_short_weather.pdf", width = 10, height = 10*1.414)
dwplot(results.direct, conf.level = 0.95, dodge_size = 0.6,
       vars_order = c("FE r", "FE growth.new.cases",  
                      "DFE r", "DFE growth.new.cases", 
                      "RE r", "RE growth.new.cases"),
       whisker_args = list(size = 1.5),
       dot_args = list(size = 3.5, shape = 20)) +
  theme_bw(base_line_size = 1.5) +
  scale_color_brewer(palette = "Dark2",
                     name = "Standard Errors",
                     guide = guide_legend(override.aes = list(linetype = 1,
                                                              shape = NA))) +
  theme(text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        panel.grid.major.y = element_blank())+
  xlab("Estimated direct effect of mask policy and corresponding 95%-CI") +
  scale_x_continuous(breaks = c(0.1,0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7), limits = c(-0.7, 0.08),
                     minor_breaks = c(0.05,-0.05,-0.15,-0.25,-0.35,-0.45,-0.55,-0.65)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) 
dev.off()

# Total effect
pdf(".\\Plots\\ci_plot_final_total_short_weather.pdf", width = 10, height = 10*1.414)
dwplot(results.total, conf.level = 0.95, dodge_size = 0.6, 
       vars_order = c("FE r", "FE growth.new.cases",  
                      "DFE r", "DFE growth.new.cases", 
                      "RE r", "RE growth.new.cases"),
       whisker_args = list(size = 1.5),
       dot_args = list(size = 3.5, shape = 20)) +
  theme_bw(base_line_size = 1.5) +
  scale_color_brewer(palette = "Dark2",
                     name = "Standard Errors",
                     guide = guide_legend(override.aes = list(linetype = 1,
                                                              shape = NA))) +
  theme(text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        panel.grid.major.y = element_blank()) +
  xlab("Estimated total effect of mask policy and corresponding 95%-CI") +
  scale_x_continuous(breaks = c(0.1,0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7), limits = c(-0.7, 0.08),
                     minor_breaks = c(0.05,-0.05,-0.15,-0.25,-0.35,-0.45,-0.55,-0.65)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) 
dev.off()

# (ii) Compare R and Case-growth

# Read in data and prepare it
data.r = data.prep(lag = 7, shift = 14, response = "median_R_mean", r.infovar = 21, frequency = "weekly")$data
data.r = data.r %>%
  arrange(X.Canton_3, X.oneweek) %>%
  select(X.Canton_3, X.oneweek, Y, X.population_we)
data.r = rename(data.r, median_R_mean = Y)

# Shift = 11 for maximal correlation
data.case = data.prep(lag = 7, shift = 14, response = "casegrowth", r.infovar = 21, frequency = "weekly")$data
data.case = data.case %>%
  arrange(X.Canton_3, X.oneweek) %>%
  select(X.Canton_3, X.oneweek, Y)
data.case = rename(data.case, casegrowth = Y)

# Merge data and adjust data type and input the abbreviations of the cantons. Also compute correlation
data = merge(data.r, data.case)
dates = rep(seq(as.Date("2020-07-06"), as.Date("2020-10-18"), by = "week"), 26)
data = data %>%
  group_by(X.Canton_3) %>%
  mutate(correl = sprintf("italic(rho) == %.2f", round(cor(median_R_mean, casegrowth), 2)),
         correl.num = round(cor(median_R_mean, casegrowth), 2))
data = cbind(data, dates)
colnames(data) = c("X.Canton_3", "X.oneweek", "median_R_mean", "X.population_we", "casegrowth", "correl", "correl.num", "date")

# Keep first (one) observation per canton to compute the correlation between the correlation of R and case-growth
# and the population size
data.corr = data %>%
  group_by(X.Canton_3) %>% 
  filter(row_number() == 1)

corr.pop = cor(data.corr$correl.num, data.corr$X.population_we) # 0.57
median.correl = median(data.corr$correl.num) # 0.585
mean = mean(data.corr$correl.num) # 0.56

data = data %>%                               
  mutate(X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 2, "GE"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 4, "VD"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 5, "FR"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 11, "NE"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 3, "VS"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 9, "BE"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 17, "JU"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 16, "SO"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 22, "BL"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 12, "LU"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 26, "SH"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 24, "BS"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 20, "AG"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 8, "OW"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 10, "NW"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 6, "UR"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 1, "TI"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 15, "ZG"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 13, "SZ"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 14, "GL"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 23, "SG"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 18, "ZH"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 21, "AR"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 19, "AI"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 7, "GR"),
         X.Canton_3 = replace(X.Canton_3, X.Canton_3 == 25, "TG"))

# Plot
Sys.setlocale("LC_TIME", "English")
pdf(".\\Plots\\r-case-compare_short.pdf", width = 10, height = 10*1.414)
ggplot(mapping = aes(x = date, y = median_R_mean, group = X.Canton_3)) +
  geom_line(data = data, aes(x = date, y = median_R_mean, group = X.Canton_3, colour = "median_R_mean"), size = 0.75) +
  geom_line(data = data, aes(x = date, y = casegrowth, group = X.Canton_3, colour = "casegrowth"), size = 0.75) +
  facet_wrap(~ X.Canton_3, ncol = 3) +
  geom_text(x = as.Date("2020-08-10"), y = 3.6, aes(label = correl), parse = TRUE, data = data, size = 4.5,
            check_overlap = TRUE) +
  xlab("") +
  ylab("") +
  ylim(c(-4,6)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%B") +
  scale_color_manual(name = "", values = c("median_R_mean" = "#1B9E77", 
                                           "casegrowth" = "#E7298A"),
                     labels = c("r", "growth.new.cases")) +
  theme_bw() +
  theme(text = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 14),
        strip.text.y = element_text(color = "black", size = 12),
        legend.position = c(0.8275, 0.04)) 
dev.off()

# (iii) Pairwise correlation plot of variables

# Read in data
data.corrplot = data.prep(lag = 7, shift = 14, response = "median_R_mean", r.infovar = 21, frequency = "daily")$data
data.corrplot = data.corrplot %>%
  select(-c(X.day, X.Canton_3, Y, X.median_R_mean.lag)) %>%
  relocate(W, .before = X.schoolClosing) %>%
  relocate(X.percentage_age, .after = X.transactiongrowth) %>%
  relocate(X.Density, .after = X.percentage_age) %>%
  relocate(X.population, .after = X.Density) %>%
  relocate(X.workplaces, .after = X.transit_stations)


# Compute a measure of correlation between the variables. For two continuous variables, we use the Pearson correlation (in [-1,1]).
# For two discrete variables, we use Cramers V and (in [0,1]). Lastly, we use use a simple regression of the continuous variable
# on the discrete variable to obtain the R^2, which we multiply with the sign of the coefficient (in [-1,1]).

# We first determine which variables are continuous and which are discrete
vars.cont = c("X.grocery_and_pharmacy",  
              "X.transit_stations", "X.transactiongrowth", "X.workplaces",
              "X.sre000d0", "X.tre200d0", "X.ure200d0",
              "X.population", "X.Density", "X.percentage_age")
vars.disc = c("X.restGatherings", "X.schoolClosing", "X.cancEvents", "X.testingPolicy",
              "X.workClosing2a", "W", "X.ferien")

# Allocation
corrmat = matrix(NA, nrow = ncol(data.corrplot), ncol = ncol(data.corrplot))
colnames(corrmat) = rownames(corrmat) = colnames(data.corrplot)

# Build the correlation measure matrix
for (j in colnames(corrmat)) {
  
  for (k in rownames(corrmat)) {
    
    if (j %in% vars.cont && k %in% vars.cont) {
      
      corrmat[j,k] = cor(data.corrplot[,j], data.corrplot[,k])
    }
    
    if (j %in% vars.cont && k %in% vars.disc | k %in% vars.cont && j %in% vars.disc) {
      
      r.squared = summary(lm(data.corrplot[,j] ~ data.corrplot[,k]))$r.squared
      sign.coef = summary(lm(data.corrplot[,j] ~ data.corrplot[,k]))$coefficients[2]/abs(summary(lm(data.corrplot[,j] ~ data.corrplot[,k]))$coefficients[2])
      corrmat[j,k] = sqrt(r.squared)*sign.coef
      
    }
    
    if (j %in% vars.disc && k %in% vars.disc) {
      
      corrmat[j,k] = cramerV(x = data.corrplot[,j], y = data.corrplot[,k])
    }
  }
}  

# Correlation with itself
diag(corrmat) = 1

# Nice column names for the plot
colnames(corrmat) = rownames(corrmat) = c("Mask policy", "School closings", "Restrictions on gatherings",
                                          "Cancellation of events", "Testing policy", "Work closings", "Mobility in grocery and pharmacy",
                                          "Mobility in transit stations", "Mobility in workplaces", "Growth rate of transactions", "Percentage of old", 
                                          "Population density", "Population", "Sunshine", "Temperature", "Humidity", "Holidays")

# Plot
pdf(".\\Plots\\corrmat.pdf", width = 12, height = 12)
corrplot(corrmat, method = "color", addCoef.col = "black", addgrid.col = "gray50", pch.cex = 2,
         tl.col = "black", number.cex = 0.8, number.font = 5, diag = FALSE)
dev.off()

# (iv) Diagnostics for all models in one plot

# We create a plot where have fitted vs. residuals plots for all the fixed effects approaches. We therefore create 4 such plots
# being FE Daily R, FE Daily Case-growth, FE Weekly R and FE Weekly Case-growth. The plots from the de-biased models do not differ and
# do hence not offer more insight. We also do not include the additional information variables and we focus on the total effect.



# Parameters to loop over
freq.poss = c("daily", "weekly")
response.poss = c("median_R_mean", "casegrowth")

# Allocation of results
results.list.fe = vector("list", 4)
i = 1


for (frequency in freq.poss) {
  
  for (response in response.poss) {
    
    results.list.fe[[i]] = estimation(frequency = frequency, response = response, model = "within", infovar = FALSE, type.effect = "total")$fit
    i = i+1
  }
}

# Run the function diagnostics for all the fits
for (k in 1:length(results.list.fe)) {
  
  if (k == 1) {
    
    response = "median_R_mean"
    frequency = "daily"
  } else if (k == 2) {
    
    response = "casegrowth"
    frequency = "daily" 
  } else if (k == 3) {
    
    response = "median_R_mean"
    frequency = "weekly"
  } else if (k == 4) {
    
    response = "casegrowth"
    frequency = "weekly"
  }
  
  # Run diagnostics
  diag.res = diagnostics(fit = results.list.fe[[k]], frequency = frequency, response = response, infovar = FALSE, type.effect = "total")
  data = diag.res[["data"]]
  outliers = data[names(diag.res[["influential"]]),]
  id = rownames(outliers)
  
  data$infl = ifelse(rownames(data) %in% id, 1, 0)
  infl = as.vector(data$infl)
  
  # Data of fitted and residuals with indicator if observation is influential
  data.plot = data.frame(diag.res[["preds"]], diag.res[["res"]], infl)
  colnames(data.plot) = c("Predicted.values", "Residuals", "Influential")
  
  # Plot fitted vs. residuals
  filename = paste("fittedres", k, ".pdf", sep = "")
  pdf(filename, width = 6, height = 6)
  print(ggplot(data = data.plot) +
          geom_point(aes(x = Predicted.values, y = Residuals),
                     position = position_jitter(h = 0.1, w = 0.1), alpha = 0.5, size = 3,
                     colour = "gray70") +
          geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
          geom_smooth(data = data.plot, mapping = aes(x = Predicted.values, y = Residuals), method = "loess", se = FALSE, formula = y ~ x,
                      colour = "#E7298A", fullrange = TRUE) +
          theme_bw() +
          xlab("Predicted values") +
          ylab("Residuals") +
          theme(text = element_text(color = "black", size = 12),
                legend.text = element_text(color = "black", size = 12),
                plot.title = element_text(hjust = 0.5)))
  dev.off()
}

# New plots!

# Total effect
# Parameters to loop over
names = c("RE r total", "RE growth.new.cases total", "FE r total", "FE growth.new.cases total", "DFE r total", "DFE growth.new.cases total")
model.poss = c("random", "within")
response.poss = c("median_R_mean", "casegrowth")

# Allocation of results
results.list.ta = vector("list", 4)
i = 1

# Fixed effects model and random effects model
for (model in model.poss) {
  
  for (response in response.poss) {
    
    results.list.ta[[i]] = estimation(frequency = "weekly", response = response, model = model, infovar = FALSE, type.effect = "total")$fit
    i = i+1
  }
}


# Debiased fixed effects
for (response in response.poss) {
  
  results.list.ta[[i]] = list(multiple_split(response = response, frequency = "weekly", infovar = FALSE, type.effect = "total"),
                              frequency = "weekly",
                              response = response,
                              model = "debiased",
                              type.effect = "total")
  i = i+1
}


# Create the plots
for (k in 1:length(results.list.ta)) {
  
  if (k <= 4) {
    
    residuals = as.numeric(results.list.ta[[k]][["residuals"]])
    predicted = as.numeric(results.list.ta[[k]][["model"]][["Y"]]) - residuals
    
  } else if (k == 5) {
    
    # Reconstruct residuals and fitted values for R.
    data = tw.demean(data.prep(lag = 7, shift = 14, response = "median_R_mean", frequency = "weekly", r.infovar = 21)$data,
                     response = "median_R_mean", frequency = "weekly")
    vars = results.list.ta[[3]][["model"]]
    X = as.matrix(vars[,-1])
    cols.X = colnames(X)
    X.demeaned = as.matrix(data[,(colnames(data) %in% cols.X)])
    X.demeaned = as.matrix(X.demeaned[,c(1,2,6,7,8,3,4,5,10,9)]) # Order
    y.demeaned = data[,"Y"]
    
    coef = as.vector(results.list.ta[[5]][[1]])
    predicted = as.numeric(X.demeaned%*%coef)
    residuals = as.numeric(y.demeaned - predicted)
    
  } else if (k == 6) {
    
    
    # Reconstruct residuals and fitted values for R.
    data = tw.demean(data.prep(lag = 7, shift = 14, response = "casegrowth", frequency = "weekly", r.infovar = 21)$data,
                     response = "casegrowth", frequency = "weekly")
    vars = results.list.ta[[4]][["model"]]
    X = as.matrix(vars[,-1])
    cols.X = colnames(X)
    X.demeaned = as.matrix(data[,(colnames(data) %in% cols.X)])
    X.demeaned = as.matrix(X.demeaned[,c(1,2,6,7,8,3,4,5,10,9)]) # Order
    y.demeaned = data[,"Y"]
    
    coef = as.vector(results.list.ta[[5]][[1]])
    predicted = as.numeric(X.demeaned%*%coef)
    residuals = as.numeric(y.demeaned - predicted)
  }
  
  # create dataframe for ggplot
  data.plot = data.frame(predicted, residuals)
  colnames(data.plot) = c("Predictions", "Residuals")
  
  # Plot fitted vs. residuals
  filename = paste("fittedres_total_short", k, ".pdf", sep = "")
  pdf(filename, width = 6, height = 6)
  print(ggplot(data = data.plot) +
          geom_point(aes(x = predicted, y = residuals),
                     position = position_jitter(h = 0.1, w = 0.1), alpha = 0.5, size = 3,
                     colour = "gray70") +
          geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
          theme_bw() +
          xlab("Predicted values") +
          ylab("Residuals") +
          ggtitle(paste(names[k])) + 
          theme(text = element_text(color = "black", size = 12),
                legend.text = element_text(color = "black", size = 12),
                plot.title = element_text(hjust = 0.5)))
  dev.off()
  
}

# Direct effect
# Parameters to loop over
names = c("RE r direct", "RE growth.new.cases direct", "FE r direct", "FE growth.new.cases direct", "DFE r direct", "DFE growth.new.cases direct")
model.poss = c("random", "within")
response.poss = c("median_R_mean", "casegrowth")

# Allocation of results
results.list.ta = vector("list", 4)
i = 1

# Fixed effects model and random effects model
for (model in model.poss) {
  
  for (response in response.poss) {
    
    results.list.ta[[i]] = estimation(frequency = "weekly", response = response, model = model, infovar = FALSE, type.effect = "direct")$fit
    i = i+1
  }
}


# Debiased fixed effects
for (response in response.poss) {
  
  results.list.ta[[i]] = list(multiple_split(response = response, frequency = "weekly", infovar = FALSE, type.effect = "direct"),
                              frequency = "weekly",
                              response = response,
                              model = "debiased",
                              type.effect = "total")
  i = i+1
}

# Create the plots
for (k in 1:length(results.list.ta)) {
  
  if (k <= 4) {
    
    residuals = as.numeric(results.list.ta[[k]][["residuals"]])
    predicted = as.numeric(results.list.ta[[k]][["model"]][["Y"]]) - residuals
    
  } else if (k == 5) {
    
    # Reconstruct residuals and fitted values for R.
    data = tw.demean(data.prep(lag = 7, shift = 14, response = "median_R_mean", frequency = "weekly", r.infovar = 21)$data,
                     response = "median_R_mean", frequency = "weekly")
    vars = results.list.ta[[3]][["model"]]
    X = as.matrix(vars[,-1])
    cols.X = colnames(X)
    X.demeaned = as.matrix(data[,(colnames(data) %in% cols.X)])
    X.demeaned = as.matrix(X.demeaned[,c(1,2,7,8,9,3,4,5,11,10,6)]) # Order
    y.demeaned = data[,"Y"]
    
    coef = as.vector(results.list.ta[[5]][[1]])
    predicted = as.numeric(X.demeaned%*%coef)
    residuals = as.numeric(y.demeaned - predicted)
    
  } else if (k == 6) {
    
    
    # Reconstruct residuals and fitted values for R.
    data = tw.demean(data.prep(lag = 7, shift = 14, response = "casegrowth", frequency = "weekly", r.infovar = 21)$data,
                     response = "casegrowth", frequency = "weekly")
    vars = results.list.ta[[4]][["model"]]
    X = as.matrix(vars[,-1])
    cols.X = colnames(X)
    X.demeaned = as.matrix(data[,(colnames(data) %in% cols.X)])
    X.demeaned = as.matrix(X.demeaned[,c(1,2,7,8,9,3,4,5,11,10,6)]) # Order
    y.demeaned = data[,"Y"]
    
    coef = as.vector(results.list.ta[[6]][[1]])
    predicted = as.numeric(X.demeaned%*%coef)
    residuals = as.numeric(y.demeaned - predicted)
  }
  
  # create dataframe for ggplot
  data.plot = data.frame(predicted, residuals)
  colnames(data.plot) = c("Predictions", "Residuals")
  
  # Plot fitted vs. residuals
  filename = paste("fittedres_direct_short", k, ".pdf", sep = "")
  pdf(filename, width = 6, height = 6)
  print(ggplot(data = data.plot) +
          geom_point(aes(x = predicted, y = residuals),
                     position = position_jitter(h = 0.1, w = 0.1), alpha = 0.5, size = 3,
                     colour = "gray70") +
          geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
          theme_bw() +
          xlab("Predicted values") +
          ylab("Residuals") +
          ggtitle(paste(names[k])) +
          theme(text = element_text(color = "black", size = 12),
                legend.text = element_text(color = "black", size = 12),
                plot.title = element_text(hjust = 0.5)))
  dev.off()
}


# (v) P-value heat-map

# Read in data; change col-names to correspond to dw-plot; make ordering consistent
results = read.csv(".\\Data\\results-05.csv", header = T, sep = ";", stringsAsFactors = FALSE)
results = results %>%                               
  mutate(model = replace(model, model == "Spatial", "Own"))
results.direct = results[results[,"type.effect"] == "direct",]
results.total = results[results[,"type.effect"] == "total",]
results.direct = results.direct[,-6]
results.total = results.total[,-6]
results.total = results.total[c(1:4,6,5,7:nrow(results.total)),]
colnames(results.direct) = colnames(results.total) = c("estimate", "std.error", "p.val", "model", "term")

# Reshape for heat-map. Do it for the total effect
data.pval = reshape(results.total, direction = "wide", idvar = "term", timevar = "model")
order = c("FE Daily R", "FE Daily Casegrowth", "FE Weekly R", "FE Weekly Casegrowth", 
          "DFE Daily R", "DFE Daily Casegrowth", "DFE Weekly R", "DFE Weekly Casegrowth",
          "RE Daily R", "RE Daily Casegrowth", "RE Weekly R", "RE Weekly Casegrowth")
data.pval = data.pval %>%
  slice(match(order, term))
row.names(data.pval) = data.pval[,1]
data.pval = data.pval[,-1]
data.pval = data.pval %>%
  select(starts_with("p.val"))

# Names 
colnames(data.pval) = c("HC3", "Canton", "Time", "Canton-Time", "NW", "Own", "CH", "Canton-Bootstrap")

# Heat-map
pdf(".\\Plots\\pval.pdf", width = 12, height = 12)
par(mar=c(5.1,10.5,4.1,2.1))
plot(as.pvalue(as.matrix(data.pval)), fmt.cell ="%.6f", reorder = FALSE, na.print = FALSE, main = "",
     col = brewer.pal(5, "Reds"), fmt.key = "%.2f",
     xlab = "", ylab = "", na.col = "white", key = list(side = 3, cex.axis = 0.75),
     polygon.key = NULL, axis.key = NULL, spacing.key = c(0.8,0.4,-0.5),
     border = "black", axis.row = list(side = 2, las = 1))
dev.off()

# (VI) Create overview table of variables and add from which data source them stem from.
# Do that for both the r-value and case-growth. Use all the variables, i.e the set of covariates used in the 
# random effects models.

for (response in c("median_R_mean", "casegrowth")) {
  
  # Store data
  data = data.prep(lag = 7, shift = 14, response = response, r.infovar = 21, frequency = "weekly")$data
  
  # Remove unused variables
  cols.to.drop = c("X.Canton_3", "X.oneweek")
  data.use = data[,!names(data) %in% cols.to.drop]
  
  # Extract the summary statistics for the variables  
  sumstats = matrix(NA, nrow = 4, ncol = ncol(data.use))
  colnames(sumstats) = colnames(data.use)
  rownames(sumstats) = c("median", "mean", "min", "max")
  
  for (j in 1:ncol(data.use)) {
    sumstats[1,j] = median(data.use[,j], na.rm = T)
    sumstats[2,j] = mean(data.use[,j], na.rm = T)
    sumstats[3,j] = min(data.use[,j], na.rm = T)
    sumstats[4,j] = max(data.use[,j], na.rm = T)
  }
  
}


# (VI) Create plot for the adaptation of the strict facial mask policy of the cantons

# Read in data
DataCovid = read.csv(".\\Data\\DataCovid.csv") 
dat.new = read.csv(".\\Data\\policy_stuff.csv", header = TRUE, sep = ",")

# Merge data on common days
data = merge(DataCovid, dat.new)
data = data %>% 
  select(geoRegion, facialCover, datum) %>%
  mutate(datum = as.Date(datum),
         facialCover = facialCover - 2) %>%
  filter(datum >= "2020-07-06" & datum <= "2020-10-25")

# Plot
pdf(".\\Plots\\policy_fraction_short.pdf", width = 20, height = 10)
ggplot(data = data, mapping = aes(x = datum, y = facialCover, group = geoRegion)) +
  geom_line(size = 1.5) +
  geom_line(data = dat.new %>% mutate(datum = as.Date(datum),
                                      facialCover = facialCover - 2) %>% filter(datum >= "2020-07-06" & datum <= "2020-10-25" & geoRegion == "CH"),
            mapping = aes(x = datum, y = facialCover), color = "firebrick", size = 2.5) +
  geom_segment(aes(x = as.Date(c("2020-08-20")), y = 0.79-1, xend  = as.Date(c("2020-08-20")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 0.75) +
  geom_segment(aes(x = as.Date(c("2020-08-23")), y = 0.69-1, xend  = as.Date(c("2020-08-23")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  geom_segment(aes(x = as.Date(c("2020-09-16")), y = 0.89-1, xend  = as.Date(c("2020-09-16")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  geom_segment(aes(x = as.Date(c("2020-10-11")), y = 0.89-1, xend  = as.Date(c("2020-10-11")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  geom_segment(aes(x = as.Date(c("2020-10-13")), y = 0.79-1, xend  = as.Date(c("2020-10-13")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  geom_segment(aes(x = as.Date(c("2020-10-15")), y = 0.69-1, xend  = as.Date(c("2020-10-15")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  geom_segment(aes(x = as.Date(c("2020-10-16")), y = 0.59-1, xend  = as.Date(c("2020-10-16")), yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  geom_segment(aes(x = as.Date(c("2020-10-17"))+0.1, y = 0.49-1, xend  = as.Date(c("2020-10-17"))+0.1, yend = 0.99-1), data = data,
               colour = "grey", linetype = 1, size = 1) +
  annotate("text", x = as.Date(c("2020-08-20")), y = 0.75-1, label = c("NE"), size = 8) +
  annotate("text", x = as.Date(c("2020-08-23")), y = 0.65-1, label = c("BS"), size = 8) +
  annotate("text", x = as.Date(c("2020-09-16")), y = 0.85-1, label = c("VD"), size = 8) +
  annotate("text", x = as.Date(c("2020-10-11")), y = 0.85-1, label = c("BE"), size = 8) +
  annotate("text", x = as.Date(c("2020-10-13")), y = 0.75-1, label = c("GE"), size = 8) +
  annotate("text", x = as.Date(c("2020-10-15")), y = 0.65-1, label = c("SZ"), size = 8) +
  annotate("text", x = as.Date(c("2020-10-16")), y = 0.55-1, label = c("LU, GR, FR"), size = 8) +
  annotate("text", x = as.Date(c("2020-10-17")), y = 0.45-1, label = c("VS"), size = 8) +
  scale_x_date(breaks = "1 month", date_labels = "%B") +
  scale_y_continuous(breaks= c(-1,0,1)) +
  labs(x = "", y = "Facial-Mask Policy", color = "Cantons") +
  theme_minimal() +
  theme(text = element_text(size = 30), 
        panel.grid = element_line(linetype = 2, size = 1.25),
        panel.grid.minor = element_line(linetype = 2, size = 1.25),
        panel.border = element_blank()) 
dev.off()

################################################################################




