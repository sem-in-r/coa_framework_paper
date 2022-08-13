# Load packages and libraries
# install.packages("seminr")
# install.packages("rpart")
# install.packages("gesca")
# install.packages("fancyRpartPlot")
# install.packages("party")
# install.packages("tree")
# devtools::install_github("https://github.com/sem-in-r/semcoa", ref = "ms-coa-paper", force = TRUE)
library(party)
library(tree)
library(rpart)
library(seminr)
library(semcoa)
library(gesca)

# Load the project data ----
utaut_216 <- read.csv(file = "utaut2_216.csv")

# Structural Model ----
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

# Output for figure 6 ----
pls_boot_summary$bootstrapped_paths

# Conduct the coa framework analysis
coa_output <- coa(pls_model, focal_construct = "BI", deviance_bounds = )

# output for Figure 7
plot_pd(coa_output)

# Identify deviant groups and members
names(coa_output$dtree$deviant_groups)
coa_output$dtree$deviant_groups

# output for figure 8 and Table 1 ----
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
fancyRpartPlot(GSCA_tree)
fancyRpartPlot(coa_output$deviance_tree)

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
  
# code for figure D1 ----
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
# figure D1 ----
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

# Table E2 can be dervied from the following plots:
fancyRpartPlot(coa_output$deviance_tree)
fancyRpartPlot(utaut_overfit_no_gender$deviance_tree)
fancyRpartPlot(utaut_overfit_no_exp$deviance_tree)
fancyRpartPlot(utaut_overfit_no_age$deviance_tree)

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
plot(partytree)
partytree@where
partytree@tree

# rpart comparative tree
rparttree <- coa(pls_model = pls_model, focal_construct = "BI")
fancyRpartPlot(rparttree$deviance_tree)
rparttree$deviance_tree


# tree comparative tree
treecontrol <- tree.control(216, mincut = 1, minsize = 2, mindev = 0)
treetree <- tree(PD ~ ., data = cs_data, control = treecontrol, split = "gini")
treetree$frame[287,1:2] 
plot(treetree, type = "uniform")
cutoff <- quantile(PD, probs = c(0.025, 0.975))
treetree$frame[treetree$frame[,"yval"] < cutoff[[1]] | treetree$frame[,"yval"] > cutoff[[2]],]
treetree$frame[treetree$frame[,"yval"] < cutoff[[1]] | treetree$frame[,"yval"] > cutoff[[2]],]
treetreeframe <- cbind(treetree$frame, 1:nrow(treetree$frame))
whichcase <- function(node) {(1:216)[treetree$where == treetreeframe[rownames(treetreeframe) == node,6]]}
whichcase(117)
whichrule <- function(node) {treetreeframe[rownames(treetreeframe) == node,]}
whichrule(3)

# Compare deviant clusters ----
# 99
## party
(1:216)[partytree@get_where() == partytree@get_where()[99]]
### Cases 99, 106, 166 cluster together in node 19
# rules:
# 27
## party
(1:216)[partytree@get_where() == partytree@get_where()[32]]
### Cases 99, 106, 166 cluster together in node 19
# rules:

# 12
## party
(1:216)[partytree@get_where() == partytree@get_where()[12]]
### Cases 12 and 127 cluster together in node 18
# rules:
# BI < 0.372
# PE > -0.668
# BI < 0.651
# HM < 0.591
# PE < 0.508
# EE < -2.2

## rpart
cstreecomp$frame[cstreecomp$where[12],]
cstreecomp$where == 34
### Cases 12 and 99 cluster together in node 239
# SI < 0.037
# FC < -0.93
# FC > -2.4
# PE > -1.1
# HAB < 1
# BI < 0.46

#####################################################
# 187
## party
(1:216)[partytree@get_where() == partytree@get_where()[187]]
### Cases 63 and 187 cluster together in node 45
# rules:
# EXP < -1.323
# HAB > 0.997
# BI > -0.651
# PE > -0.668
# BI < 0.372

## rpart
cstreecomp$frame[cstreecomp$where[187],]
(1:216)[cstreecomp$where == 37]
### Cases 106 and 187 cluster together in node 31
# FC < -1.3
# HAB > 1
# FC > -2.4
# BI < -0.46

#####################################################
# 81
## party
(1:216)[partytree@get_where() == partytree@get_where()[81]]
### Cases 81 and 109 cluster together in node 49
# rules:
# SI < -1.996
# BI < 1.572
# BI > 0.372

## rpart
cstreecomp$frame[cstreecomp$where[81],]
(1:216)[cstreecomp$where == 3]
### Cases 81 and 109 cluster together in node 4
# SI < -1.6
# BI > 0.46

#####################################################
# 134
## party
(1:216)[partytree@get_where() == partytree@get_where()[134]]
### Cases 93  96 134 151 cluster together in node 68
# rules:
# HAB < 0.194
# BI > 1.572
# BI > 0.372

## rpart
cstreecomp$frame[cstreecomp$where[134],]
(1:216)[cstreecomp$where == 6]
### Cases 93  96 134 151 cluster together in node 20
# HAB < -0.22
# BI > 1.8
# SI > -1.6
# BI > 0.46

#####################################################
# 1
## party
(1:216)[partytree@get_where() == partytree@get_where()[1]]
### Cases  1  31  38  44  72  95 108 145 202 205 cluster together in node 66
# rules:
# FC > 0.172
# SI > 0.7
# HAB > 0.15
# SI > -1.996
# BI < 1.572
# BI > 0.372

## rpart
cstreecomp$frame[cstreecomp$where[1],]
(1:216)[cstreecomp$where == 12]
### Cases 1   3  10  14  22  25  31  38  41  44  45  47  48  50  57  60  61  64  70  72  77  80
#        86  90  95  97 103 104 108 118 119 126 131 138 144 145 148 149 150 177 178 179 202 205
#       207 215 cluster together in node 23
# HM > -0.11
# HAB > 0.22
# SI > -1.66
# BI > 0.46
