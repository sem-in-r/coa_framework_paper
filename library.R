
# Function to randomly sample and return index
getRandomIndex <- function(d) {return(sample.int(nrow(d),replace = TRUE))}

# Method to calculate cooks distance for a pls model
cook_dist_pls <- function(index, data, mm, sm, construct) {
  full_model <- estimate_pls(data = data,
                             measurement_model = mm,
                             structural_model = sm)
  cook_model <- estimate_pls(data = data[-index,],
                             measurement_model = mm,
                             structural_model = sm)
  numerator <- sum(((full_model$construct_scores %*% full_model$path_coef)[-index,construct] - (cook_model$construct_scores %*% cook_model$path_coef)[,construct])^2)
  denominator <- (ncol(full_model$construct_scores)+1)* ( mean((full_model$construct_scores[,construct] - (full_model$construct_scores %*% full_model$path_coef)[,construct])^2) )
  return(numerator/denominator)
}

source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}


rmse <- function(vect) {
  return(sqrt(mean(vect^2)))
}

### predict_pls ----
predict_pls_construct <- function (model, technique = predict_DA, noFolds = NULL, cores = NULL)
{
  stopifnot(inherits(model, "seminr_model"))
  order <- sample(nrow(model$data), nrow(model$data), replace = FALSE)
  ordered_data <- model$data[order, ]
  pred_matrices <- prediction_matrices_composite(noFolds, ordered_data,
                                       model, technique, cores)
  PLS_predicted_outsample_construct <- pred_matrices$out_of_sample_construct[as.character(c(1:nrow(model$data))), ]
  PLS_predicted_insample_construct <- pred_matrices$in_sample_construct[as.character(c(1:nrow(model$data))), ]
  PLS_predicted_outsample_item <- pred_matrices$out_of_sample_item[as.character(c(1:nrow(model$data))),]
  PLS_predicted_insample_item <- pred_matrices$in_sample_item[as.character(c(1:nrow(model$data))), ]
  LM_predicted_outsample_item <- pred_matrices$out_of_sample_lm_item[as.character(c(1:nrow(model$data))),]
  LM_predicted_insample_item <- pred_matrices$in_sample_lm_item[as.character(c(1:nrow(model$data))), ]
  results <- list(composites = list(composite_out_of_sample = PLS_predicted_outsample_construct,
                                    composite_in_sample = PLS_predicted_insample_construct,
                                    actuals_star = model$construct_scores[as.character(c(1:nrow(model$data))),
                                    ]), items = list(item_out_of_sample = PLS_predicted_outsample_item,
                                                     item_in_sample = PLS_predicted_insample_item, lm_out_of_sample = LM_predicted_outsample_item,
                                                     lm_in_sample = LM_predicted_insample_item, item_actuals = ordered_data[as.character(c(1:nrow(model$data))),
                                                     ], lm_in_sample_residuals = pred_matrices$lm_in_sample_item_residuals[as.character(c(1:nrow(model$data))),
                                                     ], pls_in_sample_residuals = pred_matrices$pls_in_sample_item_residuals[as.character(c(1:nrow(model$data))),
                                                     ]))
  class(results) <- "pls_prediction_kfold"
  return(results)
}

prediction_matrices_composite <- function (noFolds, ordered_data, model, technique, cores)
{
  out <- tryCatch({
    if (is.null(noFolds)) {
      noFolds = nrow(model$data)
      folds <- cut(seq(1, nrow(ordered_data)), breaks = noFolds,
                   labels = FALSE)
      suppressWarnings(ifelse(is.null(cores), cl <- parallel::makeCluster(parallel::detectCores()),
                              cl <- parallel::makeCluster(cores)))
      parallel::clusterExport(cl = cl, varlist = c("generate_lm_predictions",
                                                   "predict_lm_matrices",
                                                   "predict_composite"), envir = environment())
      utils::capture.output(matrices <- parallel::parSapply(cl,
                                                            1:noFolds, in_and_out_sample_predictions_composite, folds = folds,
                                                            ordered_data = ordered_data, model = model, technique = technique))
      parallel::stopCluster(cl)
    }
    else {
      folds <- cut(seq(1, nrow(ordered_data)), breaks = noFolds,
                   labels = FALSE)
      matrices <- sapply(1:noFolds, in_and_out_sample_predictions_composite,
                         folds = folds, ordered_data = ordered_data, model = model,
                         technique = technique)
    }
    in_sample_construct_matrix <- do.call(cbind, matrices[(1:(noFolds *
                                                                8))[1:(noFolds * 8)%%8 == 1]])
    out_sample_construct_matrix <- do.call(cbind, matrices[(1:(noFolds *
                                                                 8))[1:(noFolds * 8)%%8 == 2]])
    in_sample_item_matrix <- do.call(cbind, matrices[(1:(noFolds *
                                                           8))[1:(noFolds * 8)%%8 == 3]])
    out_sample_item_matrix <- do.call(cbind, matrices[(1:(noFolds *
                                                            8))[1:(noFolds * 8)%%8 == 4]])
    in_sample_lm_matrix <- do.call(cbind, matrices[(1:(noFolds *
                                                         8))[1:(noFolds * 8)%%8 == 5]])
    out_sample_lm_matrix <- do.call(cbind, matrices[(1:(noFolds *
                                                          8))[1:(noFolds * 8)%%8 == 6]])
    PLS_in_sample_item_residuals <- do.call(cbind, matrices[(1:(noFolds *
                                                                  8))[1:(noFolds * 8)%%8 == 7]])
    LM_in_sample_item_residuals <- do.call(cbind, matrices[(1:(noFolds *
                                                                 8))[1:(noFolds * 8)%%8 == 0]])
    average_insample_construct <- sapply(1:length(model$constructs),
                                         seminr:::mean_rows, matrix = in_sample_construct_matrix, noFolds = noFolds,
                                         constructs = model$constructs)
    average_insample_item <- sapply(1:length(model$mmVariables),
                                    seminr:::mean_rows, matrix = in_sample_item_matrix, noFolds = noFolds,
                                    constructs = model$mmVariables)
    average_outsample_construct <- sapply(1:length(model$constructs),
                                          seminr:::sum_rows, matrix = out_sample_construct_matrix, noFolds = noFolds,
                                          constructs = model$constructs)
    average_outsample_item <- sapply(1:length(model$mmVariables),
                                     seminr:::sum_rows, matrix = out_sample_item_matrix, noFolds = noFolds,
                                     constructs = model$mmVariables)
    average_insample_pls_item_residuals <- sqrt(sapply(1:length(model$mmVariables),
                                                       seminr:::mean_rows, matrix = PLS_in_sample_item_residuals^2,
                                                       noFolds = noFolds, constructs = model$mmVariables))
    endogenous_items <- unlist(sapply(unique(model$smMatrix[,
                                                            2]), function(x) model$mmMatrix[model$mmMatrix[,
                                                                                                           "construct"] == x, "measurement"]), use.names = FALSE)
    average_insample_lm <- sapply(1:length(endogenous_items),
                                  seminr:::mean_rows, matrix = in_sample_lm_matrix, noFolds = noFolds,
                                  constructs = endogenous_items)
    average_outsample_lm <- sapply(1:length(endogenous_items),
                                   seminr:::sum_rows, matrix = out_sample_lm_matrix, noFolds = noFolds,
                                   constructs = endogenous_items)
    average_insample_lm_item_residuals <- sqrt(sapply(1:length(endogenous_items),
                                                      seminr:::mean_rows, matrix = LM_in_sample_item_residuals^2,
                                                      noFolds = noFolds, constructs = endogenous_items))
    colnames(average_insample_construct) <- colnames(average_outsample_construct) <- model$constructs
    colnames(average_insample_item) <- colnames(average_insample_pls_item_residuals) <- colnames(average_outsample_item) <- model$mmVariables
    colnames(average_insample_lm) <- colnames(average_outsample_lm) <- colnames(average_insample_lm_item_residuals) <- endogenous_items
    return(list(out_of_sample_construct = average_outsample_construct,
                in_sample_construct = average_insample_construct,
                out_of_sample_item = average_outsample_item, in_sample_item = average_insample_item,
                out_of_sample_lm_item = average_outsample_lm, in_sample_lm_item = average_insample_lm,
                pls_in_sample_item_residuals = average_insample_pls_item_residuals,
                lm_in_sample_item_residuals = average_insample_lm_item_residuals))
  }, error = function(cond) {
    message("Parallel encountered this ERROR: ")
    message(cond)
    parallel::stopCluster(cl)
    return(NULL)
  }, warning = function(cond) {
    message("Parallel encountered this WARNING:")
    message(cond)
    parallel::stopCluster(cl)
    return(NULL)
  }, finally = {
  })
}

generate_lm_predictions <- function (x, model, ordered_data, testIndexes, endogenous_items,
                                     trainIndexes)
{
  dependant_items <- model$mmMatrix[model$mmMatrix[, 1] ==
                                      x, 2]
  in_sample_matrix <- matrix(0, nrow = nrow(ordered_data),
                             ncol = length(dependant_items), dimnames = list(rownames(ordered_data),
                                                                             dependant_items))
  out_sample_matrix <- matrix(0, nrow = nrow(ordered_data),
                              ncol = length(dependant_items), dimnames = list(rownames(ordered_data),
                                                                              dependant_items))
  independant_matrix <- ordered_data[, -which(names(ordered_data) %in%
                                                dependant_items)]
  dependant_matrix <- as.matrix(ordered_data[, dependant_items])
  indepTestData <- independant_matrix[testIndexes, ]
  indepTrainData <- independant_matrix[-testIndexes, ]
  if (length(testIndexes) == 1) {
    depTestData <- t(as.matrix(dependant_matrix[testIndexes,
    ]))
  }
  else {
    depTestData <- as.matrix(dependant_matrix[testIndexes,
    ])
  }
  depTrainData <- as.matrix(dependant_matrix[-testIndexes,
  ])
  colnames(depTrainData) <- colnames(depTestData) <- dependant_items
  lm_prediction_list <- sapply(dependant_items, predict_lm_matrices,
                               depTrainData = depTrainData, indepTrainData = indepTrainData,
                               indepTestData = indepTestData, endogenous_items = endogenous_items)
  in_sample_matrix[trainIndexes, ] <- matrix(unlist(lm_prediction_list[(1:length(lm_prediction_list))[1:length(lm_prediction_list)%%2 ==
                                                                                                        1]]), ncol = length(dependant_items), nrow = nrow(depTrainData),
                                             dimnames = list(rownames(depTrainData), dependant_items))
  out_sample_matrix[testIndexes, ] <- matrix(unlist(lm_prediction_list[(1:length(lm_prediction_list))[1:length(lm_prediction_list)%%2 ==
                                                                                                        0]]), ncol = length(dependant_items), nrow = nrow(depTestData),
                                             dimnames = list(rownames(depTestData), dependant_items))
  return(list(in_sample_matrix, out_sample_matrix))
}


predict_lm_matrices <- function (x, depTrainData, indepTrainData, indepTestData, endogenous_items)
{
  trainLM <- stats::lm(depTrainData[, x] ~ ., indepTrainData)
  lmprediction_out_sample <- stats::predict(trainLM, newdata = indepTestData)
  lmprediction_in_sample <- stats::predict(trainLM, newdata = indepTrainData)
  return(list(lm_prediction_in_sample = lmprediction_in_sample,
              lm_prediction_out_sample = lmprediction_out_sample))
}

in_and_out_sample_predictions_composite <- function (x, folds, ordered_data, model, technique)
{
  testIndexes <- which(folds == x, arr.ind = TRUE)
  trainIndexes <- which(folds != x, arr.ind = TRUE)
  testingData <- ordered_data[testIndexes, ,drop = F]
  trainingData <- ordered_data[-testIndexes, ,drop = F]
  PLS_predicted_outsample_construct <- matrix(0, nrow = nrow(ordered_data),
                                              ncol = length(model$constructs), dimnames = list(rownames(ordered_data),
                                                                                               model$constructs))
  PLS_predicted_insample_construct <- matrix(0, nrow = nrow(ordered_data),
                                             ncol = length(model$constructs), dimnames = list(rownames(ordered_data),
                                                                                              model$constructs))
  PLS_predicted_outsample_item <- matrix(0, nrow = nrow(ordered_data),
                                         ncol = length(model$mmVariables), dimnames = list(rownames(ordered_data),
                                                                                           model$mmVariables))
  PLS_predicted_insample_item <- matrix(0, nrow = nrow(ordered_data),
                                        ncol = length(model$mmVariables), dimnames = list(rownames(ordered_data),
                                                                                          model$mmVariables))
  PLS_predicted_insample_item_residuals <- matrix(0, nrow = nrow(ordered_data),
                                                  ncol = length(model$mmVariables), dimnames = list(rownames(ordered_data),
                                                                                                    model$mmVariables))
  utils::capture.output(train_model <- seminr::estimate_pls(data = trainingData,
                                                            measurement_model = model$measurement_model, structural_model = model$smMatrix,
                                                            inner_weights = model$inner_weights))
  test_predictions <- predict_composite(object = train_model,
                                     testData = testingData, technique = technique)
  PLS_predicted_outsample_construct[testIndexes, ] <- test_predictions$predicted_constructs
  PLS_predicted_outsample_item[testIndexes, ] <- test_predictions$predicted_items
  train_predictions <- predict_composite(object = train_model,
                                      testData = trainingData, technique = technique)
  PLS_predicted_insample_construct[trainIndexes, ] <- train_predictions$predicted_constructs
  PLS_predicted_insample_item[trainIndexes, ] <- train_predictions$predicted_items
  PLS_predicted_insample_item_residuals[trainIndexes, ] <- as.matrix(train_predictions$item_residuals)
  endogenous_items <- unlist(sapply(unique(model$smMatrix[,
                                                          2]), function(x) model$mmMatrix[model$mmMatrix[, "construct"] ==
                                                                                            x, "measurement"]), use.names = FALSE)
  lm_holder <- sapply(unique(model$smMatrix[, 2]), generate_lm_predictions,
                      model = model, ordered_data = ordered_data[, model$mmVariables],
                      testIndexes = testIndexes, endogenous_items = endogenous_items,
                      trainIndexes = trainIndexes)
  lmprediction_in_sample <- matrix(0, ncol = 0, nrow = length(trainIndexes))
  lmprediction_out_sample <- matrix(0, ncol = 0, nrow = length(testIndexes))
  lmprediction_in_sample_residuals <- matrix(0, nrow = nrow(ordered_data),
                                             ncol = length(endogenous_items), byrow = TRUE, dimnames = list(rownames(ordered_data),
                                                                                                            endogenous_items))
  lmprediction_in_sample <- do.call(cbind, lm_holder[((1:(length(unique(model$smMatrix[,
                                                                                       2])) * 2))[1:(length(unique(model$smMatrix[, 2])) * 2)%%2 ==
                                                                                                    1])])
  lmprediction_out_sample <- do.call(cbind, lm_holder[((1:(length(unique(model$smMatrix[,
                                                                                        2])) * 2))[1:(length(unique(model$smMatrix[, 2])) * 2)%%2 ==
                                                                                                     0])])
  lmprediction_in_sample_residuals[trainIndexes, ] <- as.matrix(ordered_data[trainIndexes,
                                                                             endogenous_items]) - lmprediction_in_sample[trainIndexes,
                                                                                                                         endogenous_items]
  return(list(PLS_predicted_insample = PLS_predicted_insample_construct,
              PLS_predicted_outsample = PLS_predicted_outsample_construct,
              PLS_predicted_insample_item = PLS_predicted_insample_item,
              PLS_predicted_outsample_item = PLS_predicted_outsample_item,
              LM_predicted_insample_item = lmprediction_in_sample,
              LM_predicted_outsample_item = lmprediction_out_sample,
              PLS_predicted_insample_item_residuals = PLS_predicted_insample_item_residuals,
              LM_predicted_insample_item_residuals = lmprediction_in_sample_residuals))
}

predict_composite <- function(object, testData, technique = predict_DA, na.print=".", digits=3, ...){
  stopifnot(inherits(object, "seminr_model"))
  
  # Abort if received a higher-order-model or moderated model
  if (!is.null(object$hoc)) {
    message("There is no published solution for applying PLSpredict to higher-order-models")
    return()
  }
  if (!is.null(object$interaction)) {
    message("There is no published solution for applying PLSpredict to moderated models")
    return()
  }
  
  #Extract Measurements needed for Predictions
  normData <- testData[,object$mmVariables]
  
  # Standardize data
  normData[,object$mmVariables] <- seminr:::standardize_data(normData[,object$mmVariables],object$meanData[object$mmVariables],object$sdData[object$mmVariables])
  
  #Convert dataset to matrix
  normData<-data.matrix(normData)
  
  #Estimate Factor Scores from Outter Path
  predicted_construct_scores <- normData%*%object$outer_weights
  
  #Estimate Factor Scores from Inner Path and complete Matrix
  predicted_construct_scores <- technique(object$smMatrix, object$path_coef, predicted_construct_scores)
  
  #Predict Measurements with loadings
  predictedMeasurements<-predicted_construct_scores%*% t(object$outer_loadings)
  
  # Unstandardize data
  predictedMeasurements[,object$mmVariables] <- seminr:::unstandardize_data(predictedMeasurements[,object$mmVariables],object$meanData[object$mmVariables],object$sdData[object$mmVariables])
  
  #Calculating the residuals
  residuals <- testData[,object$mmVariables] - predictedMeasurements[,object$mmVariables]
  
  #Prepare return Object
  predictResults <- list(testData = testData[,object$mmVariables],
                         predicted_items = predictedMeasurements[,object$mmVariables],
                         item_residuals = residuals,
                         predicted_constructs = predicted_construct_scores)
  
  class(predictResults) <- "predicted_seminr_model"
  return(predictResults)
}
