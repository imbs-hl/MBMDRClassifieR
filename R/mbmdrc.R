
#' MB-MDR based classifier
#'
#' Use the MB-MDR HLO classification of top models to predict the disease risk.
#'
#' @param formula       		    	[\code{formula}]\cr
#'                      		    	A formula of the form LHS ~ RHS, where LHS is the dependent
#'                      		    	variable and RHS are the features to enter the MB-MDR
#'                      		     	genotype classification.
#' @param data          			    [\code{data.frame}]\cr
#'                      			    The data object.
#' @param order         			    [\code{int} or \code{integer}]\cr
#'                      			    Single integer specifying the interaction depth used
#'                      			    in the MB-MDR algorithm. Possible options are
#'                      			    1 and 2.
#' @param min_cell_size           [\code{int}]\cr
#'                                Single integer specifying the minimum number of
#'                                observation in a genotype combination to be statistically
#'                                relevant.
#' @param alpha         			    [\code{double}]\cr
#'                      			    Single numeric as significance level used during HLO
#'                      			    classification of genotype combinations.#'
#' @param adjustment              [\code{string}]\cr
#'                                Adjust method to be used.
#'                                "CODOMINANT", "ADDITIVE" or "NONE" (default).
#' @param max_results   			    [\code{int}]\cr
#'                      			    Single integer specifying the number of top models to
#'                      			    report.
#' @param top_results   			    [\code{int}]\cr
#'                      			    Single integer specifying how many models shall be used
#'                      			    for prediction. If \code{folds} and \code{cv.loss} are
#'                      			    set, this value specifies the upper limit of top results
#'                      			    to assess in CV.
#' @param folds         			    [\code{int}]\cr
#'                      			    Single interger specifying how many folds for internal
#'                      			    cross validation to find optimal \code{top_results}
#'                      			    should be used.
#' @param cv_loss       			    [\code{character}]\cr
#'                      			    One of \code{auc} or \code{bac}, specifying which loss
#'                      			    should be used to find optimal \code{top_results}.
#' @param dependent_variable_name [\code{string}]\cr
#' 									              Name of dependent variable, needed if no formula given.
#' @param verbose                 [\code{int}]\cr
#' 									              Level of verbosity. Default is level 1 giving some basic
#' 									              information about progress. Level 0 will switch off any
#' 									              output. Levels 2 and 3 are for debugging purposes.
#' @param log_file					      [\code{string}]\cr
#' 									              File path for a log file. If empty (default) logging is
#' 									              sent to console.
#' @param ...                     Further arguments passed to or from other methods.
#'
#' @return A S3 object of class \code{mbmdrc}.
#'
#' @details If the data type of the dependent variable is a factor, classification mode is started automatically. Otherwise MB-MDR will run in regression mode.
#' In classification mode, the first factor level is assumed to code for the negative class.
#'
#' @export
#' @import data.table
mbmdrc <- function(formula, data,
                   order = 2,
                   min_cell_size = 10,
                   alpha = 0.1,
                   adjustment = "NONE",
                   max_results = 1000,
                   top_results = 1000,
                   folds, cv_loss,
                   dependent_variable_name,
                   verbose,
                   ...) {

  # Input checks
  assertions <- checkmate::makeAssertCollection()

  checkmate::assertDataFrame(data, min.cols = 2, min.rows = 2,
                             add = assertions)
  checkmate::assertInteger(order, lower = 1, upper = 2,
                           max.len = 2, min.len = 1,
                           any.missing = FALSE, unique = TRUE, null.ok = FALSE,
                           add = assertions)
  checkmate::assertInt(min_cell_size, lower = 0,
                       add = assertions)
  checkmate::assertNumber(alpha, lower = 0, upper = 1,
                          add = assertions)
  checkmate::assertChoice(adjustment, c("CODOMINANT", "ADDITIVE", "NONE"),
                          add = assertions)
  checkmate::assertInt(max_results, lower = 1, upper = 1e4)
  checkmate::assertInt(top_results, lower = 1, upper = max_results)

  if(missing(verbose)) {
    verbose <- 1
  } else {
    checkmate::assertInt(verbose, lower = 0, upper = 3,
                         add = assertions)
  }

  # Formula interface ----
  if(missing(formula)) {
    if(missing(dependent_variable_name)) {
      assertions$push("Please give formula or dependent variable name.")
    }
    checkmate::assertString(dependent_variable_name,
                            add = assertions)
    checkmate::assertChoice(dependent_variable_name, choices = colnames(data),
                            add = assertions)

    response <- data[, dependent_variable_name]
    data_selected <- data[, -which(colnames(data)==dependent_variable_name)]
  } else {
    formula <- stats::as.formula(formula)
    checkmate::assertClass(formula, classes = "formula",
                           add = assertions)
    data_selected <- stats::model.frame(formula, data, na.action = stats::na.pass)
    response <- data_selected[[1]]
  }

  # Prediction type ----
  if(is.factor(response)) {
    checkmate::assertFactor(response, n.levels = 2,
                            add = assertions)
    model_type <- "binary" # classification
  } else {
    model_type <- "continuous" # regression
  }

  # Interaction depth ----
  dim <- sapply(order, function(o) {
    switch(o,
           "1" = "1D",
           "2" = "2D",
           "3" = "3D")
  })

  checkmate::reportAssertions(assertions)

  # Dependent variable name ----
  if(!missing(formula)) {
    dependent_variable_name <- names(data_selected)[1]
    independent_variable_names <- names(data_selected)[-1]
  } else {
    indepentend_variable_names <- colnames(data_selected)[colnames(data_selected) != dependent_variable_name]
  }

  # Input data and variable names, create final data matrix
  if(is.matrix(data_selected)) {
    data_final <- data_selected
  } else {
    data_final <- data.matrix(data_selected)[, -1]
  }
  variable_names <- colnames(data_final)

  # Write MB-MDR file ----
  file <- tempfile()
  data.table::fwrite(data.frame("y" = response, data_final),
                     file = file,
                     sep = " ",
                     append = FALSE,
                     na = "-9")

  # Clean up ----
  rm("data_selected")

  # Initialize output ----
  result <- list()

  # Data dependent max_results and top_results
  max_results <- min(sum(sapply(order, function(o) choose(ncol(data_final), o))), max_results)
  top_results <- min(max_results, top_results)

  # Internal cross validation ----
  if(!missing(folds) & !missing(cv_loss)) {
    checkmate::assertInt(folds, na.ok = TRUE, null.ok = TRUE, lower = 2, upper = 10,
                         add = assertions)
    checkmate::assertChoice(cv_loss, choices = c("auc", "bac"),
                            add = assertions)

    # Select prediction type
    pred_type <- switch(cv_loss,
                        "auc" = "scoreprob",
                        "bac" = "response")

    # Select loss function
    measure <- switch(cv_loss,
                      "auc" = auc,
                      "bac" = bac)

    # Set search space for optimal top_results value
    if(length(top_results) == 1) {
      top_results <- 1:max(top_results)
    }

    fold_idx <- sample(1:folds, nrow(data_final), replace = TRUE)

    # Calculate the MB-MDR for each fold and assess current top_results value
    cv_performance <- rbindlist(lapply(1:folds, function(f) {
      # Generate MB-MDR models on CV training data
      data_cv_train <- data_final[fold_idx != f,]
      response_cv_train <- response[fold_idx != f]
      cv_file <- tempfile()
      data.table::fwrite(data.frame("y" = response_cv_train, data_cv_train),
                         file = cv_file,
                         sep = " ",
                         append = FALSE,
                         na = "-9")

      mbmdr <- do.call('c', sapply(dim, function(d) {
        mbmdR::mbmdr(file = cv_file, trait = model_type,
                     cpus.topfiles = 1,
                     cpus.permutations = 1,
                     work.dir = tempdir(),
                     n.pvalues = max_results,
                     permutations = 0,
                     group.size = min_cell_size,
                     alpha = alpha,
                     dim = d,
                     multi.test.corr = "NONE",
                     adjustment = adjustment,
                     verbose = "MEDIUM")
      }))

      # Predict on CV testing data
      data_cv_test <- data_final[fold_idx == f,]
      pred <- predict.mbmdr(object = mbmdr, newdata = data_cv_test,
                            all = TRUE,
                            o_as_na = FALSE,
                            type = pred_type)

      # Prepare CV loss for all top_result values
      response_cv_test <- response[fold_idx == f]
      pred <- melt(pred, id.vars = "ID",
                   variable.name = "top_results",
                   value.name = pred_type)
      pred[, list(cv_loss = measure(get(pred_type),
                                    response_cv_test,
                                    positive = 1,
                                    negative = 0)), by = "top_results"]
    }), idcol = "fold")

    # Select top_results value with optimal loss
    mean_cv_loss <- NULL # hack to circumvent R CMD CHECK notes
    top_results <- cv_performance[, list(mean_cv_loss = mean(cv_loss)),
                                  by = "top_results"][, which.max(mean_cv_loss)]

    # Save CV results
    result$cv_performance <- cv_performance
  }

  # Call MB-MDR ----
  mbmdr <- do.call('c', sapply(dim, function(d) {
               mbmdR::mbmdr(file = file, trait = model_type,
                            work.dir = tempdir(),
                            cpus.topfiles = 1, cpus.permutations = 1,
                            n.pvalues = top_results,
                            permutations = 0,
                            group.size = min_cell_size, alpha = alpha,
                            dim = d,
                            multi.test.corr = "NONE", adjustment = adjustment,
                            verbose = "MEDIUM")
             }))
  result$mbmdr <- mbmdr

  result$call <- sys.call()
  result$num_samples <- nrow(data_final)
  result$num_features <- ncol(data_final)
  result$num_combinations <- choose(ncol(data_final), order)
  result$model_type <- model_type
  result$top_results <- top_results
  if (is.factor(response)) {
    result$levels <- levels(droplevels(response))
  }

  class(result) <- "mbmdrc"

  return(result)

}

#' MB-MDR classifier prediction
#'
#' Prediction with new data and saved MB-MDR classifier object.
#'
#' @param object        [\code{mbmdr}]\cr
#'                      MB-MDR models and HLO tables as output from \code{\link{mbmdrc}}.
#' @param newdata       [\code{newdata}]\cr
#'                      New data to predict class status for.
#' @param type          [\code{string}]\cr
#'                      Type of prediction. One of \code{response}, \code{prob},
#'                      \code{score} or \code{scoreprob}. See details.
#'                      Default is \code{response}.
#' @param top_results   [\code{int}]\cr
#'                      How many models are used for prediction.
#' @param all           [\code{bool}]\cr
#' 						          Output predictions for all possible top results.
#' @param o_as_na       [\code{bool}]\cr
#'                      Encode non informative cells with NA or with 0.5.
#' @param ...           Further arguments passed to or from other methods.
#'
#' @return A \code{data.table} object with an ID column and the prediction value.
#'
#' @details For \code{type='response'} (the default), the predicted classes, for
#' \code{type='prob'} the predicted case probabilities, for \code{type='score'}
#' the predicted risk scores and for \code{type='scoreprob'} the predicted risk
#' scores transformed to the [0, 1] interval are returned.
#'
#' For \code{type='score'} and \code{type='scoreprob'} genotypes classified as
#' H contribute +1, as L contribute -1 and as O contribute 0 to the score.
#'
#' If a genotype combination is classified as O by MB-MDR, the case probability
#' is not significantly different from 0.5. On the other hand, there might have
#' been just too few observations in the training data so that \code{NA} might
#' be more reasonable as contribution to \code{response} and \code{prob} type
#' predictions.
#'
#' @import data.table
predict.mbmdr <- function(object, newdata, type = "response", top_results, all = FALSE, o_as_na = TRUE, ...) {

  # Input checks ----
  assertions <- checkmate::makeAssertCollection()
  checkmate::assertClass(object, "mbmdr",
                         add = assertions)
  checkmate::assert(checkmate::checkDataFrame(newdata),
                    checkmate::checkMatrix(newdata),
                    combine = "or")
  features <- unique(unlist(sapply(object, function(model) model$features)))
  checkmate::assertSubset(features, colnames(newdata),
                          add = assertions)
  checkmate::assertChoice(type, c("response", "prob", "score", "scoreprob"),
                          add = assertions)
  checkmate::assertFlag(all,
                        add = assertions)
  if(missing(top_results) & !all) {
    checkmate::reportAssertions(assertions)
    stop("Please specify the number of top results to enter predictions or set 'all=TRUE'")
  } else if(!missing(top_results)) {
    checkmate::assertInt(top_results,
                         lower = 0, upper = length(object),
                         na.ok = FALSE, null.ok = FALSE,
                         add = assertions)
  }

  checkmate::reportAssertions(assertions)

  # Get number of models and number of samples
  num_models <- length(object)
  num_samples <- nrow(newdata)

  # Initialize predictions
  predictions <- data.table::data.table(MODEL = rep(1:num_models, each = num_samples),
                                        ID = rep(1:num_samples, times = num_models),
                                        PROB = numeric(num_models * num_samples))

  # Iterate through MB-MDR models
  for(m in 1:num_models) {
    # Get genotypes
    genotypes <- as.matrix(newdata[, object[[m]]$features])

    # Construct bases for indexing
    num_rows <- attr(object[[m]]$cell_labels, "num_rows")
    bases <- if(num_rows == 1) {
      length(object[[m]]$cell_labels)
    } else {
      c(length(object[[m]]$cell_labels)/num_rows, num_rows)
    }

    # Construct index as linear combination of bases and feature combination
    idx <- genotypes %*% bases^(0:(length(bases)-1)) + 1

    # Get cell predictions
    prob <- object[[m]]$cell_predictions[idx]

    if(!o_as_na) {
      # Set cell predictions to 0.5 for non-informative cells or feature combinations not present in training data
      prob[object[[m]]$cell_labels[idx]=="O" | is.na(prob)] <- 0.5
    } else {
      # Set cell predictions to NA for non-informative cells or feature combinations not present in training data
      prob[object[[m]]$cell_labels[idx]=="O" | is.na(prob)] <- NA
    }

    predictions[MODEL == m, PROB := prob]
  }

  if(all) {
    switch(type,
           # Round the mean case probability to 0 or 1 to return hard classification
           "response" = return(dcast(predictions[, list(RESPONSE = round(cumsum(PROB)/1:max(MODEL)),
                                                        TOPRESULTS = 1:max(MODEL)), by = c("ID")],
                                     ID~TOPRESULTS, value.var = "RESPONSE")),
           # Return mean case probability over all models for all top.results
           "prob" = return(dcast(predictions[, list(PROB = cumsum(PROB)/1:max(MODEL),
                                                    TOPRESULTS = 1:max(MODEL)), by = c("ID") ],
                                 ID~TOPRESULTS, value.var = "PROB")),
           # Return a risk score. Genotype combinations classified as H contribute +1,
           # genotype combinations classified as L contribute -1 and genotype combinations
           # classified as O contribute 0
           "score" = return(dcast(predictions[, list(SCORE = cumsum(sign(PROB-0.5)),
                                                     TOPRESULTS = 1:max(MODEL)), by = c("ID") ],
                                  ID~TOPRESULTS, value.var = "SCORE")),
           # Return the score transformed to a [0, 1] interval
           "scoreprob"= return(dcast(predictions[, list(SCORE = cumsum(sign(PROB-0.5)),
                                                        TOPRESULTS = 1:max(MODEL)), by = c("ID")][, SCOREPROB := range01(SCORE),
                                                                                                  by = c("TOPRESULTS")],
                                     ID~TOPRESULTS,
                                     value.var = "SCOREPROB")))
  } else {
    switch(type,
           # Round the mean case probability to 0 or 1 to return hard classification
           "response" = return(predictions[MODEL<=top_results, list(predictions = round(mean(PROB))),
                                           by = c("ID")]$predictions),
           # Return mean case probability over all models
           "prob" = return(predictions[MODEL<=top_results, list(predictions = mean(PROB)),
                                       by = c("ID")]$predictions),
           # Return a risk score. Genotype combinations classified as H contribute +1,
           # genotype combinations classified as L contribute -1 and genotype combinations
           # classified as O contribute 0
           "score" = return(predictions[MODEL<=top_results, list(predictions = sum(sign(PROB-0.5))),
                                        by = c("ID")]$predictions),
           # Return the score transformed to a [0, 1] interval
           "scoreprob" = return(predictions[MODEL<=top_results, list(predictions = sum(sign(PROB-0.5))),
                                            by = c("ID")][, list(ID, predictions = range01(predictions))]$predictions))
  }

}

#' @rdname predict.mbmdr
#'
#' @export
predict.mbmdrc <- function(object, newdata, type = "response", top_results, o_as_na = TRUE, ...) {

  # Input checks ----
  assertions <- checkmate::makeAssertCollection()

  checkmate::assertClass(object, "mbmdrc",
                         add = assertions)
  checkmate::assertDataFrame(newdata,
                             add = assertions)
  checkmate::assertSubset(object$mbmdr$feature_names, colnames(newdata))
  checkmate::assertChoice(type, choices = c("response", "prob", "score", "scoreprob"),
                          add = assertions)
  if(!missing(top_results)) {
    checkmate::assertInt(top_results, lower = 1, upper = length(object$mbmdr),
                         add = assertions)
    top_results <- min(length(object$mbmdr), top_results)
  } else {
    top_results <- object$top_results
  }
  checkmate::assertChoice(o_as_na,
                          add = assertions)

  checkmate::reportAssertions(assertions)

  predictions <- stats::predict(object$mbmdr, newdata = newdata,
                                model_type = model_type, type = type,
                                top_results = top_results,
                                o_as_na = o_as_na, ...)

  if(type %in% c("prob", "scoreprob")) {
    predictions <- cbind(1-predictions, predictions)
    colnames(predictions) <- object$levels
  }

  return(predictions)

}

range01 <- function(x, ...) {
  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}
