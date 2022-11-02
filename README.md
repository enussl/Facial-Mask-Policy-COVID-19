# Facial-Mask-Policy-COVID-19
Repository containing all the code used in the paper "The effect of a strict facial-mask policy on the spread of COVID-19 in Switzerland during the early phase of the pandemic" by Nussli* et al.

## Table of Contents
* [General Info](#general-info)
* [Scripts](#scripts)

## General Info
The project is structured in 3 folders being:
* Data
* Code
* Plots

Using the data provided under Data, we run all the estimations using the scripts provided under Code, which generates the results saved in Data and the
plots outputted to Plots.

## Scripts
We shortly describe the scripts and the implemented functionalities therein.
### helperfunctions.R
* panellag: compute lagged instances of data with panel structure
* paneldiff: compute differenced instances of data with panel stucture
* panelma: compute moving averaged instances of data with panel structure
* lognozero: avoid taking $\log$ of zero by using 0.1 for values smaller than 0.1
* dlogd: compute log-differenced growth rates
* data.prep: prepare data for estimation
* data.prep.halfcantons: prepare data while merging halfcantons
* construct.formula: construct regression formula
* tw.demean: perform two-way demeaning
* M.Andrews: determine cut-off lag for Chiang-Hansen covariance matrix estimator
* sd.Thompson.Hansen: compute Chiang-Hansen covariance estimator (also computes the covariance matrix introduced in https://www.sciencedirect.com/science/article/pii/S0304405X10001923)
* create.A.init: create adjacency matrix for cantons in Switzerland using the same ordering as the data frame produced by data.prep
* sd.Informal: compute own (informal) covariance matrix estimator
* compute.se: apply all the covariance matrix estimators on the fitted model
* estimation: FE and RE estimation procedure
* multiple_split: DFE estimation procedure
* dml.estim: Double Machine Learning estimation procedure

Currently not used:
* data.wb: data preparation for bootstrapped standard errors
* bootstat_fe: bootstrapped model estimation
* bs_estimation: bootstrap standard error estimation
* r.sq: $R^2$ for plm models
* diagnostics: visual diagnostics for plm models

### estimations.R
We run all the models using the functions defined in helperfunctions.R while looping over the parameters of interest.
We write the results to the file results.csv.

### plots.R
We construct the plots used in the paper, which are saved in the folder Plots.

### robustness.R
We run the robustness checks discussed in the appendix using the functions defined in helperfunctions.R


