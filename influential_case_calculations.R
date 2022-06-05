
### Residuals and Influence in Regression - Cook and Weisberg ----

## P33 Sec 2.2.3
# Equation 2.2.22
# e_oos_i <- actual_i - t(x)_i %*% B_(i)
# e_oos_i <- actual_i - y_hat_oos_i

# Equation 2.2.23
# e_oos_i <- e_is_i / (1 - V_ii)
lev_data <- pls_model$construct_scores[,-11]
leverage <- diag((lev_data %*% solve(t(lev_data) %*% lev_data)) %*% t(lev_data))

e_oos_2 <- fit_error[2]/(1 - leverage[2])
predict_error[2]
### Thus we don't need to perform LOOCV!!
mean((fit_error/(1-leverage))^2)
mean((predict_error^2))

new_matrix <- cbind(fit_error, predict_error, plspredict_model$composites$composite_out_of_sample[,"BI"], plspredict_model$composites$composite_in_sample[,"BI"], abs(plspredict_model$composites$composite_in_sample[,"BI"]-plspredict_model$composites$composite_out_of_sample[,"BI"]), abs(plspredict_model$composites$composite_in_sample[,"BI"]-plspredict_model$composites$composite_out_of_sample[,"BI"])/fit_error)

1 - (fit_error/predict_error)

### Calculate the leverage for the structural model
construct_lev_data <- pls_model$construct_scores[,-11]
construct_leverage <- diag((pls_model$construct_scores[,-11] %*% solve(t(pls_model$construct_scores[,-11]) %*% pls_model$construct_scores[,-11])) %*% t(pls_model$construct_scores[,-11]))

oobs_3_lev <- t(construct_lev_data[3,]) %*% solve(t(construct_lev_data) %*% construct_lev_data) %*% construct_lev_data[3,]


scaled_data[,pls_model$mmVariables] %*% (pls_model$outer_weights %*% pls_model$path_coef[,"BI"])

pls_model$construct_scores[,"BI"]


###### removie this---- 
measurement_model <- constructs(
  composite("PE", multi_items("PERF", 1:4)),
  composite("EE", c("PEOU1","PEOU3","PEOU5","PEOU6")),
  composite("SI", c(multi_items("NORM", 1:2),"INFL3")),
  composite("FC", multi_items("FACL", 1:4)),
  composite("HM", multi_items("MOTIV", 1:3)),
  composite("PV", multi_items("VALUE", 1:3)),
  composite("HAB", multi_items("HAB", 1:4)),
  composite("Exp", single_item("frequency")),
  composite("Age", single_item("age")),
  composite("Gender", single_item("gender"))
)
######


### Calculate the leverage for the measurement model
pls_model$mmMatrix
PE_lev_data <- as.matrix(pls_model$rawdata[,multi_items("PERF", 1:4)])
PE_leverage <- ((PE_lev_data %*% solve(t(PE_lev_data) %*% PE_lev_data)) %*% t(PE_lev_data))

EE_lev_data <- as.matrix(pls_model$rawdata[,c("PEOU1","PEOU3","PEOU5","PEOU6")])
EE_leverage <- ((EE_lev_data %*% solve(t(EE_lev_data) %*% EE_lev_data)) %*% t(EE_lev_data))

SI_lev_data <- as.matrix(pls_model$rawdata[,c(multi_items("NORM", 1:2),"INFL3")])
SI_leverage <- ((SI_lev_data %*% solve(t(SI_lev_data) %*% SI_lev_data)) %*% t(SI_lev_data))

FC_lev_data <- as.matrix(pls_model$rawdata[,multi_items("FACL", 1:4)])
FC_leverage <- ((FC_lev_data %*% solve(t(FC_lev_data) %*% FC_lev_data)) %*% t(FC_lev_data))

HM_lev_data <- as.matrix(pls_model$rawdata[,multi_items("MOTIV", 1:3)])
HM_leverage <- ((HM_lev_data %*% solve(t(HM_lev_data) %*% HM_lev_data)) %*% t(HM_lev_data))

PV_lev_data <- as.matrix(pls_model$rawdata[,multi_items("VALUE", 1:3)])
PV_leverage <- ((PV_lev_data %*% solve(t(PV_lev_data) %*% PV_lev_data)) %*% t(PV_lev_data))

HAB_lev_data <- as.matrix(pls_model$rawdata[,multi_items("HAB", 1:4)])
HAB_leverage <- ((HAB_lev_data %*% solve(t(HAB_lev_data) %*% HAB_lev_data)) %*% t(HAB_lev_data))

EXP_lev_data <- as.matrix(pls_model$rawdata[,single_item("frequency")])
EXP_leverage <- ((EXP_lev_data %*% solve(t(EXP_lev_data) %*% EXP_lev_data)) %*% t(EXP_lev_data))

AGE_lev_data <- as.matrix(pls_model$rawdata[,single_item("age")])
AGE_leverage <- ((AGE_lev_data %*% solve(t(AGE_lev_data) %*% AGE_lev_data)) %*% t(AGE_lev_data))

GEN_lev_data <- as.matrix(pls_model$rawdata[,single_item("gender")])
GEN_leverage <- ((GEN_lev_data %*% solve(t(GEN_lev_data) %*% GEN_lev_data)) %*% t(GEN_lev_data))

TOT_leverage <- PE_leverage + SI_leverage + EE_leverage + FC_leverage + HM_leverage + PV_leverage + HAB_leverage + EXP_leverage + AGE_leverage + GEN_leverage  


TOT_leverage %*% construct_lev_data
scaled_data <- as.matrix(scale(pls_model$data))
pls_model$construct_scores[,"PE"] - (scaled_data[,multi_items("PERF", 1:4)] %*% pls_model$outer_weights[1:4,"PE"])
  
((FC_lev_data %*% solve(t(FC_lev_data) %*% FC_lev_data)) %*% t(FC_lev_data)) %*% construct_lev_data


### Remove case and check impact on measurement model
# Estimating the full model
new_data <- data[-1,]
rownames(new_data)<- 1:215
pls_model_no1 <- estimate_pls(data = new_data,
                          measurement_model = measurement_model,
                          structural_model = structural_model)
plspredict_model_no1 <- predict_pls(pls_model,
                                technique = predict_DA)

pls_model$outer_weights - pls_model_no1$outer_weights



### PERF example estimating 
PE_lev_data <- as.matrix(scaled_data[,multi_items("PERF", 1:4)])
PE_lev_data2 <- as.matrix(data[,multi_items("PERF", 1:4)])
PE_leverage <- diag((PE_lev_data %*% solve(t(PE_lev_data) %*% PE_lev_data)) %*% t(PE_lev_data))
PE_leverage2 <- diag((PE_lev_data2 %*% solve(t(PE_lev_data2) %*% PE_lev_data2)) %*% t(PE_lev_data2))

PE_mm_lm <- lm(pls_model$construct_scores[,"PE"] ~ scaled_data[,"PERF1"] + scaled_data[,"PERF2"] + scaled_data[,"PERF3"] + scaled_data[,"PERF4"])
PE_mm_lm$residuals


W_hat <- pls_model$outer_weights[1:4,1]
error_hat <- pls_model$construct_scores[,"BI"] - PE_lev_data %*% pls_model$outer_weights[1:4,1]

W_oos <- W_hat - (solve(t(PE_lev_data) %*% PE_lev_data) %*% PE_lev_data[2,] %*% error_hat[2])/(1 - PE_leverage[2])
W_oos <- PE_mm_lm1$coefficients[2:5] - (solve(t(PE_lev_data) %*% PE_lev_data) %*% PE_lev_data[2,] %*% error_hat[2])/(1 - PE_leverage[2])

#### Check if it works without the scaling
construct_score1 <- scale(pls_model$construct_scores[,"BI"] * cor(pls_model$construct_scores[,"PE"], pls_model$construct_scores[,"BI"]))
construct_score2 <- scale(pls_model_no1$construct_scores[,"BI"] * cor(pls_model_no1$construct_scores[,"PE"], pls_model_no1$construct_scores[,"BI"]))
PE_leverage1 <- diag((PE_lev_data %*% solve(t(PE_lev_data) %*% PE_lev_data)) %*% t(PE_lev_data))
PE_leverage2 <- diag((PE_lev_data[-1,] %*% solve(t(PE_lev_data[-1,]) %*% PE_lev_data[-1,])) %*% t(PE_lev_data[-1,]))

PE_leverage1[2]

PE_mm_lm1 <- lm(construct_score1 ~ PE_lev_data[,"PERF1"] + PE_lev_data[,"PERF2"] + PE_lev_data[,"PERF3"] + PE_lev_data[,"PERF4"])
PE_mm_lm2 <- lm(construct_score2 ~ PE_lev_data[-1,"PERF1"] + PE_lev_data[-1,"PERF2"] + PE_lev_data[-1,"PERF3"] + PE_lev_data[-1,"PERF4"])
PE_mm_lm1$residuals

new_construct_score <- PE_lev_data %*% PE_mm_lm1$coefficients[2:5]
lm(new_construct_score ~ PE_lev_data[,"PERF1"] + PE_lev_data[,"PERF2"] + PE_lev_data[,"PERF3"] + PE_lev_data[,"PERF4"])
pls_model$construct_scores[,"PE"]

outer_weights_model <- seminr:::standardize_outer_weights(PE_lev_data, mmVariables = c("PERF1","PERF2","PERF3","PERF4"),outer_weights = PE_mm_lm1$coefficients[2:5])
PE_mm_lm2$coefficients[2:5]


