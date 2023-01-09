#' This is the main script to reproduce the results from the COA paper and its
#' appendices. Run the script from top to bottom to reproduce all results. 
#' Code for specific tables and figures in the paper can be found by searching 
#' the comments in this script (e.g., look for "Table 2" or "Figure 6").
#'
#' This script may be modified in future to improve readibility and 
#' organization -- care will be taken to ensure outputs do not change. 
#' Please see CHANGELOG file for update history.

#' Load packages and libraries (uncomment and run the first time)
# install.packages("seminr")
# install.packages("rpart")
# install.packages("rpart.plot")
# install.packages("party")
# install.packages("tree")
# install.packages("devtools")
# devtools::install_github("https://github.com/cran/gesca")
#
#' SEMCOA Package installation
#' 
#' The COA framework code is packaged separately so it can be used for 
#' reproducibility or reuse.
#'
#' Reproducing MS Paper Results: install version of COA framework code from 
#' `ms-coa-paper` branch of SEMCOA repo (uncomment and run)
# devtools::install_github("https://github.com/sem-in-r/semcoa", ref = "ms-coa-paper", force = TRUE)

#' Reusing COA framework: See SEMCOA package website for general installation 
#' instructions for reuse on new models or data: https://github.com/sem-in-r/semcoa

#' Load required libraries
library(party)
library(tree)
library(rpart)
library(seminr)
library(semcoa)
library(gesca)
library(rpart.plot)

# Load the survey data used in the empirical demonstration ----
utaut_216 <- read.csv(file = "utaut2_216.csv")

# Structural Model for UTAUT demonstration ----
structural_model <- relationships(
  paths(from = c("PE","EE","SI","FC","HM","PV","HAB","Exp","Age","Gender"), to = "BI")
)

measurement_model <- constructs(
  composite("PE", multi_items("PE", 1:4)),
  composite("EE", multi_items("EE", 1:4)),
  composite("SI", multi_items("SI", 1:3)),
  composite("FC", multi_items("FC", 1:4)),
  composite("HM", multi_items("HM", 1:3)),
  composite("PV", multi_items("PV", 1:3)),
  composite("HAB", multi_items("HAB", 1:4)),
  composite("BI", multi_items("INT", 1:3)),
  composite("Exp", single_item("experience")),
  composite("Age", single_item("age")),
  composite("Gender", single_item("gender"))
)

# Estimate the PLS Model
pls_model <- estimate_pls(data = utaut_216[,-66],
                          measurement_model = measurement_model,
                          structural_model = structural_model)

# Summarize the model ----
pls_summary <- summary(pls_model)

# Bootstrap model
pls_bootstrap <- bootstrap_model(pls_model, nboot = 2000, seed = 1234)
pls_boot_summary <- summary(pls_bootstrap, alpha = 0.01)

# Output for Figure 6 ----
pls_boot_summary$bootstrapped_paths

# Conduct the coa framework analysis
coa_output <- coa(pls_model, focal_construct = "BI", deviance_bounds = )

# output for Figure 7
plot_pd(coa_output)

# Identify deviant groups and members
names(coa_output$dtree$deviant_groups)
coa_output$dtree$deviant_groups

# output for Figure 8 and Table 1 ----
# group A
group_rules("A", coa_output$dtree)
# group B
group_rules("B", coa_output$dtree)
# group C
group_rules("C", coa_output$dtree)
# group D
group_rules("D", coa_output$dtree)
# group E
group_rules("E", coa_output$dtree)

# Collect the unstable parameters
unstable <- unstable_params(analysis = coa_output)


# output of Table 2 ----
original <- pls_model$path_coef[,"BI"]
A <- unstable$group_diffs$A$param_diffs$path_coef[,"BI"]
B <- unstable$group_diffs$B$param_diffs$path_coef[,"BI"]
C <- unstable$group_diffs$C$param_diffs$path_coef[,"BI"]
D <- unstable$group_diffs$D$param_diffs$path_coef[,"BI"]
E <- unstable$group_diffs$E$param_diffs$path_coef[,"BI"]
round((cbind(original,
             A, (A)/original, 
             B, (B)/original,
             C, (C)/original,
             D, (D)/original,  
             E, (E)/original)[1:7,]),3)

# Online Appendices results: ----
# sort the predictive deviants by absolute size and select largest (absolute) 22 cases
order <- names(sort(abs(coa_output$predictions$PD), decreasing = TRUE))[1:22]
# assign to deviants variable
predictive_deviants <- order
# assign the predictive deviance for each of the deviants
predictive_deviance <- round(coa_output$predictions$PD[order],3)

# subset the indicators of the antecedent constructs
lev_data <- as.matrix(pls_model$rawdata[,pls_model$mmVariables[-c(26,27,28)]])
# calcualte leverage
leverage <- diag(((lev_data %*% solve(t(lev_data) %*% lev_data)) %*% t(lev_data)))
# name the vector
names(leverage) <- 1:216
# assign the output to a variable
leverage_x <- as.vector(leverage[order])

# subset the indicators of the outcome constructs
mahal_data_items <- utaut_216[,pls_model$mmVariables[c(26,27,28)]]
# calculate nahalanobis distance for the indicators
mahal_items <- mahalanobis(mahal_data_items, colMeans(mahal_data_items), cov(mahal_data_items))
# name the vector
names(mahal_items) <- 1:216
# assign the output to a storage variable
mahal_y <- as.vector(mahal_items[order])

# calculate residuals
residuals <- pls_model$construct_scores[,"BI"]  - coa_output$predictions$fitted_score
# assign the output to a storage variable
residuals_BI <- residuals[order]

# online appendix Table A1 ----
round(cbind(as.numeric(order), predictive_deviance, leverage_x, mahal_y, residuals_BI),3)

# online Appendix Table C1 ----

myModel <- "
		# Measurement model 
		PE =~ PE1 + PE2 + PE3 + PE4
		EE =~ EE1 + EE2 + EE3 + EE4
		SI =~ SI1 + SI2 + SI3
		FC =~ FC1 + FC2 + FC3 + FC4
		HM =~ HM1 + HM2 + HM3
		PV =~ PV1 + PV2 + PV3
		HAB =~ HAB1 + HAB2 + HAB3 + HAB4
		BI =~ INT1 + INT2 + INT3
		Exp =~ experience
		Age =~ age
		Gender =~ gender
		
		# Structural model 
		BI ~ PE + EE + SI + FC + HM + PV + HAB + Exp + Age + Gender
"
utaut_GSCA <- gesca.run(myModel, utaut_216)
summary(utaut_GSCA)

scaled_data <- scale(utaut_216[,-32])
GSCA_construct_scores <- (scaled_data[,utaut_GSCA$wname] %*% utaut_GSCA$WR) 
GSC_actual_star <- GSCA_construct_scores[,8]
GSCA_fit <- (GSCA_construct_scores %*% t(utaut_GSCA$BR))[,8]
GSCA_prediction <- c()

for (i in 1:216) {
  utaut_GSCA_train <- gesca.run(myModel, utaut_216[-i,], nbt = 1)
  scaled_data_train <- scale(utaut_216[-i,-32])
  means_vector <- attr(scaled_data_train, "scaled:center")
  sd_vector <- attr(scaled_data_train, "scaled:scale")
  standard_data <- (utaut_216[i,-32] - means_vector) / sd_vector
  GSCA_prediction[i] <- ((as.matrix(standard_data[,utaut_GSCA_train$wname]) %*% utaut_GSCA_train$WR) %*% t(utaut_GSCA_train$BR))[,8]
  
}

GSCA_predict_MSE <- mean((GSC_actual_star - GSCA_prediction)^2)
GSCA_fit_MSE <-mean((GSC_actual_star - GSCA_fit)^2)

GSCA_overfit_ratio <- (GSCA_predict_MSE - GSCA_fit_MSE)/GSCA_fit_MSE

GSCA_PD <- GSCA_fit - GSCA_prediction
names(GSCA_PD) <- 1:216

PD <- coa_output$predictions$PD 
GSCA_deviants <- (GSCA_PD)[(GSCA_PD > quantile(GSCA_PD, probs = c(0.975))) | (GSCA_PD < quantile(GSCA_PD, probs = c(0.025)))]
PLS_deviants <- (PD)[(PD > quantile(PD, probs = c(0.975))) | (PD < quantile(PD, probs = c(0.025)))]

# Table C1 ----
round(matrix(c(coa_output$predictions$IS_MSE,
coa_output$predictions$OOS_MSE,
coa_output$predictions$overfit_ratio,
GSCA_fit_MSE, 
GSCA_predict_MSE, 
GSCA_overfit_ratio,
coa_output$predictions$IS_MSE - GSCA_fit_MSE,
coa_output$predictions$OOS_MSE - GSCA_predict_MSE,
coa_output$predictions$overfit_ratio - GSCA_overfit_ratio),ncol = 3, dimnames = list(c("MSEin", "MSEout", "overfit ratio"), c("PLS", "GSCA", "PLS - GSCA"))),3)

PD_order <- c(12,32,37,40,81,93,99,106,109,134,151,180,187)
round(cbind(PD_order, 
      coa_output$predictions$PD[PD_order],
      GSCA_PD[PD_order],
      coa_output$predictions$PD[PD_order] - GSCA_PD[PD_order]),3)

gsca_comp_data <- as.data.frame(cbind(GSCA_construct_scores, GSCA_PD))
colnames(gsca_comp_data) <- c(utaut_GSCA$lname, "PD")
GSCA_tree <- rpart(PD ~ ., 
                   data = gsca_comp_data, 
                   minsplit = 2, 
                   minbucket = 1, cp = 0.001)
# Table C2 can be derived from the following plots
# Table C2 ----

rpart.plot(GSCA_tree, type = 2)
rpart.plot(coa_output$dtree$tree , type =2 )

# Compare nodes to PLS ----
quant <- quantile(GSCA_PD, probs = c(0.025,0.975))
# Node 4
GSCA_tree$frame[40,]
(1:216)[GSCA_tree$where == 40]
# Cases 81 and 109
GSCA_tree$frame[GSCA_tree$where[32],]


# Checking whether the stable paths are affected
GSCA_param_diffs <- function(remove_cases, gsca_model, gsca_data, params="BR") {
  no_dgroup_data <- gsca_data[-remove_cases,]
  suppressMessages(
    no_dgroup_model <- gesca.run(myModel,
                                 data=no_dgroup_data,
                                 nbt = 1
    )
  )
  diffs <-    gsca_model[["BR"]][8,] - no_dgroup_model[["BR"]][8,]
  names(diffs) <- gsca_model$lname[1:11]
  diffs
}
GSCA_no81_109 <- GSCA_param_diffs(c(81, 109), utaut_GSCA, utaut_216)
GSCA_no93_96_134_151 <- GSCA_param_diffs(c(93, 96, 134, 151), utaut_GSCA, utaut_216)
GSCA_no27_37_180 <- GSCA_param_diffs(c(27, 37, 180), utaut_GSCA, utaut_216)

# Difference between GSCA LDO and PLS LDO
# Table C3 ----
round(matrix(c(coa_output$pls_model$path_coef[1:10,11],
  utaut_GSCA$BR[8,-8], 
  -coa_output$pls_model$path_coef[1:10,11] - utaut_GSCA$BR[8,-8], 
  -coa_output$unstable$group_diffs$E$param_diffs$path_coef[1:10,11],
  GSCA_no81_109[-8],
  -coa_output$unstable$group_diffs$E$param_diffs$path_coef[1:10,11] - GSCA_no81_109[-8],
  -coa_output$unstable$group_diffs$D$param_diffs$path_coef[1:10,11],
  GSCA_no93_96_134_151[-8],
  -coa_output$unstable$group_diffs$D$param_diffs$path_coef[1:10,11] - GSCA_no93_96_134_151[-8],
  -coa_output$unstable$group_diffs$C$param_diffs$path_coef[1:10,11],
  GSCA_no27_37_180[-8],
  -coa_output$unstable$group_diffs$C$param_diffs$path_coef[1:10,11] - GSCA_no27_37_180[-8]), ncol = 12, dimnames = list(c("PE", "EE", "SI", "FC", "HM", "PV", "Hab", "Exp", "Age", "Gender"), c("Orig PLS", "Orig GSCA", "Orig PLS - GSCA", "A: PLS", "A: GSCA", "A: PLS - GSCA","B: PLS", "B: GSCA", "B: PLS - GSCA","E: PLS", "E: GSCA", "E: PLS - GSCA" ))),3)
  
# code for Figure D1 ----
# You can bypass the commented out code below and load the simulation results:
load(file = "returnlist2_11082021.rda")
# orig_cases <- 1:216
# utaut_data <- cbind(utaut_216, orig_cases)
# 
# utaut_model <- estimate_pls(data = utaut_data,
#                             measurement_model = measurement_model,
#                             structural_model = structural_model)
# 
# utaut_overfit <- coa(pls_model = utaut_model, 
#                      focal_construct = "BI",
#                      params = c("path_coef", "outer_weights", "rSquared"))
# 
# 
# # b. Sample size and stability -----
deviant_all <- (1:nrow(utaut_216))[coa_output$predictions$PD > quantile(coa_output$predictions$PD, probs = coa_output$deviance_bounds)[2] | coa_output$predictions$PD < quantile(coa_output$predictions$PD, probs = coa_output$deviance_bounds)[1]]
# 
# ### Downsample ----
# bootstrap_results2 <- matrix(NA, nrow = 15, ncol = 500)
# return_list2 <- vector(mode = "list", length = 500)
# set.seed(123)
# utaut_data_no_dev <- utaut_data[-deviant_all,]
# for (i in 1:500) {
#   new_sample <- sample(1:nrow(utaut_data_no_dev), (nrow(utaut_data_no_dev) - 11), replace = FALSE)
#   
#   # Add the deviants to bootstrapped non-deviants
#   new_data <- rbind(utaut_data_no_dev[new_sample,], utaut_data[deviant_all,])
#   # new_data <- utaut_data[new_sample,]
#   rownames(new_data) <- 1:nrow(new_data)
#   utaut_model_resample <- estimate_pls(data = new_data,
#                                        measurement_model = utaut_mm,
#                                        structural_model = utaut_sm)
#   
#   new_utaut_overfit <- coa(pls_model = utaut_model_resample, 
#                            focal_construct = "BI",
#                            params = c("path_coef", "outer_weights", "rSquared"))
#   focal_dev <- (1:nrow(new_data))[new_utaut_overfit$predictions$PD > quantile(new_utaut_overfit$predictions$PD, probs = new_utaut_overfit$deviance_bounds)[2] | new_utaut_overfit$predictions$PD < quantile(new_utaut_overfit$predictions$PD, probs = new_utaut_overfit$deviance_bounds)[1]]
#   bootstrap_results2[1:length(focal_dev),i] <- new_data[focal_dev, "orig_cases"]
#   
#   deviants_in <-  deviant_all %in% new_sample
#   
#   new_pd <- new_utaut_overfit$predictions$PD[focal_dev] 
#   names(new_pd) <- new_data[focal_dev, "orig_cases"]
#   
#   groups <- new_utaut_overfit$dtree$deviant_groups
#   for(int in names(new_utaut_overfit$dtree$deviant_groups)) {
#     groups[[int]] <- new_data[,"orig_cases"][groups[[int]]]
#   }
#   
#   unique <- new_data[,"orig_cases"][new_utaut_overfit$dtree$unique_deviants]
#   
#   return_list2[[i]] <- list(IS_MSE = new_utaut_overfit$predictions$IS_MSE,
#                             OOS_MSE = new_utaut_overfit$predictions$OOS_MSE,
#                             overfit_ratio = new_utaut_overfit$predictions$overfit_ratio,
#                             new_pds = new_pd,
#                             deviants_in = deviants_in,
#                             deviant_groups = groups,
#                             group_differences = new_utaut_overfit$unstable$group_diffs,
#                             deviant_unique = unique,
#                             unique_differences = new_utaut_overfit$unstable$unique_diffs,
#                             new_sample = new_sample) 
# }

col_vec <- (as.numeric(names(table(bootstrap_results2)) %in% deviant_all)*0.7) +0.3
# Figure D1 ----
barplot(table(bootstrap_results2[1:12,1:500]), col= col_vec, ylim = c(0,500))

# Table D1 and D2 are derived from the return_list2 object

# Code for Table E1 ----
# c. Removing controls and stability -----
### No Gender ----
all_deviants <- c(12, 32, 37, 40, 76, 81, 93, 99, 106, 109, 134, 151, 153, 180, 187)
utaut_sm_no_gender <- relationships(
  paths(from = c("PE","EE","SI","FC","HM","PV","HAB","Exp","Age"), to = "BI")
)

utaut_model_no_gender <- estimate_pls(data = utaut_216,
                                measurement_model = measurement_model,
                                structural_model = utaut_sm_no_gender)

utaut_overfit_no_gender <- coa(pls_model = utaut_model_no_gender, 
                         focal_construct = "BI",
                         params = c("path_coef", "outer_weights", "rSquared"))

# Prediction metrics
coa_output$predictions$IS_MSE - utaut_overfit_no_gender$predictions$IS_MSE
coa_output$predictions$OOS_MSE - utaut_overfit_no_gender$predictions$OOS_MSE
coa_output$predictions$overfit_ratio - utaut_overfit_no_gender$predictions$overfit_ratio

# Focal deviants
(1:nrow(utaut_216))[utaut_overfit_no_gender$predictions$PD > quantile(utaut_overfit_no_gender$predictions$PD, probs = utaut_overfit_no_gender$deviance_bounds)[2] | utaut_overfit_no_gender$predictions$PD < quantile(utaut_overfit_no_gender$predictions$PD, probs = utaut_overfit_no_gender$deviance_bounds)[1]]
coa_output$predictions$PD[all_deviants] - utaut_overfit_no_gender$predictions$PD[all_deviants]

# Cluster
utaut_overfit_no_gender$dtree$deviant_groups
utaut_overfit_no_gender$dtree$dev_group_rules
utaut_overfit_no_gender$dtree$unique_deviants

# Leave deviants out
utaut_overfit_no_gender$unstable$group_diffs

### No Age ----
utaut_sm_no_age <- relationships(
  paths(from = c("PE","EE","SI","FC","HM","PV","HAB","Exp","Gender"), to = "BI")
)

# Estimate model and run deviance trees

utaut_model_no_age <- estimate_pls(data = utaut_216,
                                measurement_model = measurement_model,
                                structural_model = utaut_sm_no_age)

utaut_overfit_no_age <- coa(pls_model = utaut_model_no_age, 
                         focal_construct = "BI",
                         params = c("path_coef", "outer_weights", "rSquared"))

# Prediction metrics
coa_output$predictions$IS_MSE - utaut_overfit_no_age$predictions$IS_MSE
coa_output$predictions$OOS_MSE - utaut_overfit_no_age$predictions$OOS_MSE
coa_output$predictions$overfit_ratio - utaut_overfit_no_age$predictions$overfit_ratio

# Focal deviants
(1:nrow(utaut_216))[utaut_overfit_no_age$predictions$PD > quantile(utaut_overfit_no_age$predictions$PD, probs = utaut_overfit_no_age$deviance_bounds)[2] | utaut_overfit_no_age$predictions$PD < quantile(utaut_overfit_no_age$predictions$PD, probs = utaut_overfit_no_age$deviance_bounds)[1]]
round(coa_output$predictions$PD[all_deviants] - utaut_overfit_no_age$predictions$PD[all_deviants],3)

# Cluster
coa_output$dtree$deviant_groups
utaut_overfit_no_age$dtree$deviant_groups
utaut_overfit_no_age$dtree$dev_group_rules
utaut_overfit_no_age$dtree$unique_deviants

# Leave deviants out
utaut_overfit_no_age$unstable$group_diffs

### No Experience ----
utaut_sm_no_exp <- relationships(
  paths(from = c("PE","EE","SI","FC","HM","PV","HAB","Age","Gender"), to = "BI")
)

utaut_model_no_exp <- estimate_pls(data = utaut_216,
                                measurement_model = measurement_model,
                                structural_model = utaut_sm_no_exp)

utaut_overfit_no_exp <- coa(pls_model = utaut_model_no_exp, 
                         focal_construct = "BI",
                         params = c("path_coef", "outer_weights", "rSquared"))

# Prediction metrics
coa_output$predictions$IS_MSE - utaut_overfit_no_exp$predictions$IS_MSE
coa_output$predictions$OOS_MSE - utaut_overfit_no_exp$predictions$OOS_MSE
coa_output$predictions$overfit_ratio - utaut_overfit_no_exp$predictions$overfit_ratio

# Focal deviants
(1:nrow(utaut_216))[utaut_overfit_no_exp$predictions$PD > quantile(utaut_overfit_no_exp$predictions$PD, probs = utaut_overfit_no_exp$deviance_bounds)[2] | utaut_overfit_no_exp$predictions$PD < quantile(utaut_overfit_no_exp$predictions$PD, probs = utaut_overfit_no_exp$deviance_bounds)[1]]
round(coa_output$predictions$PD[all_deviants] - utaut_overfit_no_exp$predictions$PD[all_deviants],3)

# Cluster
utaut_overfit_no_exp$dtree$deviant_groups
utaut_overfit_no_exp$dtree$dev_group_rules
utaut_overfit_no_exp$dtree$unique_deviants

# Leave deviants out
utaut_overfit_no_exp$unstable$group_diffs

# Table E1 ----
round(matrix(c(coa_output$predictions$IS_MSE,
  coa_output$predictions$OOS_MSE,
  coa_output$predictions$overfit_ratio,
  utaut_overfit_no_gender$predictions$IS_MSE,
  utaut_overfit_no_gender$predictions$OOS_MSE,
  utaut_overfit_no_gender$predictions$overfit_ratio,
  coa_output$predictions$IS_MSE - utaut_overfit_no_gender$predictions$IS_MSE,
  coa_output$predictions$OOS_MSE - utaut_overfit_no_gender$predictions$OOS_MSE,
  coa_output$predictions$overfit_ratio - utaut_overfit_no_gender$predictions$overfit_ratio,
  utaut_overfit_no_exp$predictions$IS_MSE,
  utaut_overfit_no_exp$predictions$OOS_MSE,
  utaut_overfit_no_exp$predictions$overfit_ratio,
  coa_output$predictions$IS_MSE - utaut_overfit_no_exp$predictions$IS_MSE,
  coa_output$predictions$OOS_MSE - utaut_overfit_no_exp$predictions$OOS_MSE,
  coa_output$predictions$overfit_ratio - utaut_overfit_no_exp$predictions$overfit_ratio,
  utaut_overfit_no_age$predictions$IS_MSE,
  utaut_overfit_no_age$predictions$OOS_MSE,
  utaut_overfit_no_age$predictions$overfit_ratio,
  coa_output$predictions$IS_MSE - utaut_overfit_no_age$predictions$IS_MSE,
  coa_output$predictions$OOS_MSE - utaut_overfit_no_age$predictions$OOS_MSE,
  coa_output$predictions$overfit_ratio - utaut_overfit_no_age$predictions$overfit_ratio), ncol = 7, dimnames = list(c("MSEin", "MSEout", "overfit ratio"),c("Orig PLS", "Excl Gender", "Orig - Excl Gender", "Excl Exp", "Orig - Excl Exp","Excl Age", "Orig - Excl Age"))),3)

round(matrix(c(all_deviants, 
  coa_output$predictions$PD[all_deviants],
  utaut_overfit_no_gender$predictions$PD[all_deviants],
  coa_output$predictions$PD[all_deviants] -  utaut_overfit_no_gender$predictions$PD[all_deviants],
  utaut_overfit_no_exp$predictions$PD[all_deviants],
  coa_output$predictions$PD[all_deviants] -  utaut_overfit_no_exp$predictions$PD[all_deviants],
  utaut_overfit_no_age$predictions$PD[all_deviants],
  coa_output$predictions$PD[all_deviants] -  utaut_overfit_no_age$predictions$PD[all_deviants]), ncol = 8, dimnames = list(all_deviants, c("Case","Orig PLS", "Excl gender", "Diff","Excl Exp", "Diff","Excl Age", "Diff"))),3)

# Table E2 can be derived from the following plots:
rpart.plot(coa_output$dtree$tree)
rpart.plot(utaut_overfit_no_gender$dtree$tree)
rpart.plot(utaut_overfit_no_exp$dtree$tree)
rpart.plot(utaut_overfit_no_age$dtree$tree)

# Table F1 ----
PD <- coa_output$predictions$PD
cs_data <- as.data.frame(cbind(pls_model$construct_scores, PD))
# party
ctreecontrol <- ctree_control(minsplit = 2,
                              maxdepth = 0,
                              minbucket = 1,
                              mincriterion = 0.5
)
partytree <- ctree(
  PD ~ ., 
  data = cs_data,
  controls = ctreecontrol
)

# rpart comparative tree
rparttree <- coa(pls_model = pls_model, focal_construct = "BI")
rpart.plot(rparttree$dtree$tree)
rparttree$deviance_tree

# tree comparative tree
treecontrol <- tree.control(216, mincut = 1, minsize = 2, mindev = 0)
treetree <- tree(PD ~ ., data = cs_data, control = treecontrol, split = "gini")
plot(treetree, type = "uniform")

# Compare deviant clusters across algorithms----
# Cluster A
## party
(1:216)[partytree@get_where() == partytree@get_where()[187]]
## rpart
rparttree$dtree$deviant_groups[1]$A
## tree
treetree$frame[treetree$where[[106]],]
treetree$frame[treetree$where[[187]],]
treetree$frame["22",]
# Cases 106 and 187 belong to leaves 44 and 45 which stem from node 22.

# Cluster B
## party
(1:216)[partytree@get_where() == partytree@get_where()[99]]
## rpart
rparttree$dtree$deviant_groups[2]$B
## tree
treetree$frame[treetree$where[[12]],]
treetree$frame[treetree$where[[71]],]
treetree$frame[treetree$where[[99]],]
treetree$frame["82",]
# Cases 12, 71, and 99 belong to leaves 330, 164, and 331 which stem from 
# node 82.

# Cluster C
## party
(1:216)[partytree@get_where() == partytree@get_where()[27]]
## rpart
rparttree$dtree$deviant_groups[3]$C
## tree
treetree$frame[treetree$where[[27]],]
treetree$frame[treetree$where[[37]],]
treetree$frame[treetree$where[[180]],]
treetree$frame["8",]
# Cases 27, 37, and 180 belong to leaves 17, 33, and 32 which stem from 
# node 8.

# Cluster D
## party
(1:216)[partytree@get_where() == partytree@get_where()[96]]
## rpart
rparttree$dtree$deviant_groups[4]$D
## tree
treetree$frame[treetree$where[[93]],]
treetree$frame[treetree$where[[134]],]
treetree$frame[treetree$where[[151]],]
treetree$frame["59",]
# Cases 93, 134, and 151 belong to leaves 118, 238, and 239 which stem from 
# node 59.

# Cluster E
## party
(1:216)[partytree@get_where() == partytree@get_where()[81]]
## rpart
rparttree$dtree$deviant_groups[5]$E
## tree
treetree$frame[treetree$where[[81]],]
treetree$frame[treetree$where[[109]],]
treetree$frame["6",]
# Cases 81 and 109 belong to leaves 13 and 12 which stem from 
# node 6.
