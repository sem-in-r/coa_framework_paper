# Load packages and libraries
# installed.packages("seminr")
# installed.packages("parallel") 
# install.packages("TeachingDemos")

####   NOTE:  Please do install pls-predict from the github repo as the code requires the latest version.

library(seminr)
library(parallel)
library(TeachingDemos)
source(file = "library.R")

# Variables to be specified for model
## Constructs of interest
constructs <- c("BI")
## Number of iterations
iter <- 10000
## Seed
seed <- 1234
## Folds (Null = observations)
folds <- NULL

# Load the project data ----
data <- read.csv(file = "utaut2_216.csv")
data <- data[,-66]
# Structural Model ----
structural_model <- relationships(
  paths(from = c("PE","EE","SI","FC","HM","PV","HAB","Exp","Age","Gender"), to = "BI")
)

measurement_model <- constructs(
  composite("PE", multi_items("PERF", 1:4)),
  composite("EE", c("PEOU1","PEOU3","PEOU5","PEOU6")),
  composite("SI", c(multi_items("NORM", 1:2),"INFL3")),
  composite("FC", multi_items("FACL", 1:4)),
  composite("HM", multi_items("MOTIV", 1:3)),
  composite("PV", multi_items("VALUE", 1:3)),
  composite("HAB", multi_items("HAB", 1:4)),
  composite("BI", multi_items("INT", 1:3)),
  composite("Exp", single_item("Experience")),
  composite("Age", single_item("age")),
  composite("Gender", single_item("gender"))
)

# Estimate the PLS Model
# Estimating the full model
pls_model <- estimate_pls(data = data,
                          measurement_model = measurement_model,
                          structural_model = structural_model)

# Summarize the model ----
pls_summary <- summary(pls_model)

# Perform the predictions
plspredict_model <- predict_pls_construct(pls_model)



# Automated from here ----
set.seed(seed)

# Report the cross validated in-sample and out-sample RMSE ----
IS_RMSE <- mean((plspredict_model$composites$actuals_star[,"BI"] - plspredict_model$composites$composite_in_sample[,"BI"])^2)
OOS_RMSE <- mean((plspredict_model$composites$actuals_star[,"BI"] - plspredict_model$composites$composite_out_of_sample[,"BI"])^2)

# Calculate the influential cases ----
# Calculate the list of Delta values for constructs
BI_delta <- plspredict_model$composites$composite_in_sample[,"BI"] - plspredict_model$composites$composite_out_of_sample[,"BI"]

# BI: List the top 1%, 5% most influential using delta and perBIntile 
BI_delta_5_percent <- BI_delta[BI_delta < quantile(BI_delta, probs = c(0.025)) | BI_delta > quantile(BI_delta, probs = c(0.975))]
BI_delta < quantile(BI_delta, probs = c(0.025)) | BI_delta > quantile(BI_delta, probs = c(0.975))
sort(round(BI_delta_5_percent, digits = 3))

# Identify the poorly fit cases:
poorly_fit <- abs(plspredict_model$composites$actuals_star[,"BI"] - plspredict_model$composites$composite_in_sample[,"BI"])
poorly_fit_5_percent <- poorly_fit[poorly_fit > quantile(poorly_fit, probs = c(0.95))]

sort(round(poorly_fit_5_percent, digits = 3))


# Identify the poorly predicted cases:
poorly_predicted <- abs(plspredict_model$composites$actuals_star[,"BI"] - plspredict_model$composites$composite_out_of_sample[,"BI"])
poorly_predicted_5_percent <- poorly_predicted[poorly_predicted > quantile(poorly_predicted, probs = c(0.95))]
sort(round(poorly_predicted_5_percent, digits = 3))


# # Why are they deviant?
deviant <- data[names(BI_delta_5_percent),"CODE"]

# 
# fast_and_low_sd <- too_fast[too_fast %in% low_sd]
# deviant[deviant %in% fast_and_low_sd]

# Wobble: ----

# What paths are significant?
boot_pls <- bootstrap_model(pls_model, nboot = 1000, seed = 123)
boot_sum <- summary(boot_pls)
# Wobble from 81
pls_model_no81 <- estimate_pls(data = data[-81,],
                          measurement_model = measurement_model,
                          structural_model = structural_model)

pls_model$path_coef - pls_model_no81$path_coef
boot_pls_no81 <- bootstrap_model(pls_model_no81, nboot = 1000, seed = 123)
boot_sum_no81 <- summary(boot_pls_no81)

# pls_model$outer_weights - pls_model_no81$outer_weights
# pls_model$outer_loadings - pls_model_no81$outer_loadings
# Wobble from 99
pls_model_no99 <- estimate_pls(data = data[-99,],
                               measurement_model = measurement_model,
                               structural_model = structural_model)

pls_model$path_coef - pls_model_no99$path_coef
as.vector(pls_model$path_coef[1:10,11] - pls_model_no99$path_coef[1:10,11]) / boot_sum$bootstrapped_paths[,3]
as.vector(pls_model$outer_weights - pls_model_no99$outer_weights)
boot_pls_no99 <- bootstrap_model(pls_model_no99, nboot = 2000, seed = 123)
boot_sum_no99 <- summary(boot_pls_no99)
as.vector(pls_model$path_coef[1:10,11] - pls_model_no99$path_coef[1:10,11]) / boot_sum$bootstrapped_paths[,3]
(boot_sum$bootstrapped_weights[,1] - boot_sum_no99$bootstrapped_weights[,1])/boot_sum$bootstrapped_weights[,3]

# Wobble from 106
pls_model_no106 <- estimate_pls(data = data[-106,],
                                measurement_model = measurement_model,
                                structural_model = structural_model)

pls_model$path_coef - pls_model_no106$path_coef
boot_pls_no106 <- bootstrap_model(pls_model_no106, nboot = 1000, seed = 123)
boot_sum_no106 <- summary(boot_pls_no106)

# Wobble from removal of top 5% 
new_data <- data[-as.numeric(names(BI_delta_10_percent)),]
rownames(new_data)  <- 1:nrow(new_data)      
         
pls_model_no5perc <- estimate_pls(data = new_data,
                                measurement_model = measurement_model,
                                structural_model = structural_model)

summary(pls_model_no5perc)
pls_model$path_coef - pls_model_no5perc$path_coef
# boot_pls_no5perc <- bootstrap_model(pls_model_no5perc, nboot = 1000, seed = 123)
# boot_sum_no5perc <- summary(boot_pls_no5perc)

plspredict_model_5perc <- predict_pls(pls_model_no5perc,
                                technique = predict_DA)

# Report the cross validated in-sample and out-sample RMSE ----
sum_plspredict_model_5perc <- summary(plspredict_model_5perc)
sum_plspredict_model_5perc$composite_accuracy[c("IS_RMSE","OOS_RMSE"),constructs]

# Wobble from 12
pls_model_no12 <- estimate_pls(data = data[-12,],
                                measurement_model = measurement_model,
                                structural_model = structural_model)

pls_model_no12$path_coef
pls_model$outer_weights - pls_model_no12$outer_weights
boot_pls_no12 <- bootstrap_model(pls_model_no12, nboot = 1000, seed = 123)
boot_sum_no12 <- summary(boot_pls_no12)
as.vector(pls_model$path_coef[1:10,11] - pls_model_no12$path_coef[1:10,11]) / (boot_sum$bootstrapped_paths[,3])
(boot_sum$bootstrapped_weights[,1] - boot_sum_no12$bootstrapped_weights[,1])/boot_sum$bootstrapped_weights[,3]


# Wobble from 187
pls_model_no187 <- estimate_pls(data = data[-187,],
                               measurement_model = measurement_model,
                               structural_model = structural_model)

pls_model_no187$path_coef
boot_pls_no187 <- bootstrap_model(pls_model_no187, nboot = 2000, seed = 123)
boot_sum_no187 <- summary(boot_pls_no187)

# Wobble from 180
pls_model_no180 <- estimate_pls(data = data[-180,],
                                measurement_model = measurement_model,
                                structural_model = structural_model)

pls_model_no180$path_coef
boot_pls_no180 <- bootstrap_model(pls_model_no180, nboot = 2000, seed = 123)
boot_sum_no180 <- summary(boot_pls_no180)

# Wobble from 37
pls_model_no37 <- estimate_pls(data = data[-37,],
                                measurement_model = measurement_model,
                                structural_model = structural_model)

pls_model_no37$path_coef
boot_pls_no37 <- bootstrap_model(pls_model_no37, nboot = 1000, seed = 123)
boot_sum_no37 <- summary(boot_pls_no37)


diff_vec <- abs(plspredict_model$composites$composite_in_sample[,"BI"] - plspredict_model$composites$composite_out_of_sample[,"BI"])
infl_cases <- data[diff_vec > quantile(diff_vec, probs = 0.9),"CODE"]

# CE: Calculate Cooks distance for this PLS model
BI_cooks_dist_pls <- sapply(X = 1:nrow(data),FUN = cook_dist_pls, data = data, 
                            mm = measurement_model,
                            sm = structural_model,
                            construct = "BI")

#CE: Calculate the top 1%, 5% most influential using Cooks and percentile 
names(BI_cooks_dist_pls) <- 1:nrow(data)
BI_cook_5_percent <- BI_cooks_dist_pls[BI_cooks_dist_pls > quantile(BI_cooks_dist_pls, probs = c(0.95))]

# 4/n rule 
round(BI_cooks_dist_pls[BI_cooks_dist_pls > 4/216], digits = 3)
round(BI_cook_5_percent, digits = 3)

# BI: 5% most extreme observations in Cooks Distance not in influential cases - "33"  "99"  "112"
setdiff(names(BI_cook_5_percent), names(BI_delta_5_percent))

# BI: 5% most extreme observations in influential cases not in  Cooks Distance - "12"  "185" "192"
setdiff(names(BI_delta_5_percent),names(BI_cook_5_percent))

