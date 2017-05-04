
#' MB-MDR based classifier
#'
#' Use the MB-MDR HLO classification of top models to predict the disease risk.
#'
#' @param formula       		    	[\code{formula}]\cr
#'                      		    	A formula of the form LHS ~ RHS, where LHS is the disease
#'                      			    status and RHS are the features to enter the MB-MDR
#'                      		     	genotype classification.
#' @param data          			    [\code{data.frame}]\cr
#'                      			    The data object.
#' @param order         			    [\code{int}]\cr
#'                      			    Single integer specifying the interaction depth used
#'                      			    in the MB-MDR algorithm.
#' @param alpha         			    [\code{double}]\cr
#'                      			    Single numeric as significance level used during HLO
#'                      			    classification of genotype combinations.
#' @param max_results   			    [\code{int}]\cr
#'                      			    Single integer specifying the number of top models to
#'                      			    report.
#' @param top_results   			    [\code{int} or \code{integer}]\cr
#'                      			    Single integer specifying how many models shall be used
#'                      			    for prediction. If \code{folds} and \code{cv.loss} are
#'                      			    set, this value specifies the upper limit of top results
#'                      			    to assess in CV, or, if a \code{integer} vector the top
#'                      			    result values to assess in CV.
#' @param folds         			    [\code{int}]\cr
#'                      			    Single interger specifying how many folds for internal
#'                      			    cross validation to find optimal \code{top_results}
#'                      			    should be used.
#' @param cv_loss       			    [\code{character}]\cr
#'                      			    One of \code{auc} or \code{bac}, specifying which loss
#'                      			    should be used to find optimal \code{top_results}.
#' @param dependent_variable_name [\code{string}]\cr
#' 									              Name of dependent variable, needed if no formula given.
#' @param num_threads			      	[\code{int}]\cr
#' 									              Number of threads. Default is number of CPUs available.
#' @param verbose                 [\code{int}]\cr
#' 									              Level of verbosity. Default is level 1 giving some basic
#' 									              information about progress. Level 0 will switch off any
#' 									              output. Levels 2 and 3 are for debugging purposes.
#' @param log_file					      [\code{string}]\cr
#' 									              File path for a log file. If empty (default) logging is
#' 									              sent to console.
#'
#' @return A S3 object of class \code{mbmdrc}.
#'
#' @export
#' @import data.table
#' @useDynLib MBMDRClassifieR
#' @importFrom Rcpp sourceCpp
mbmdrc <- function(formula, data,
                   order = 2,
                   alpha = 0.1,
                   max_results = 100,
                   top_results = 10,
                   folds, cv_loss,
                   dependent_variable_name,
                   num_threads,
                   verbose,
                   log_file) {

  # Input checks
  assertions <- checkmate::makeAssertCollection()

  checkmate::assertDataFrame(data, min.cols = 2, min.rows = 2,
                             add = assertions)
  checkmate::assertInt(order, lower = 2, upper = 5,
                       add = assertions)
  checkmate::assertNumber(alpha, lower = 0, upper = 1,
                          add = assertions)
  checkmate::assertInt(max_results, lower = 1, upper = 1e4)
  checkmate::assertInt(top_results, lower = 1, upper = max_results)

  if(missing(verbose)) {
    verbose <- 1
  } else {
    checkmate::assertInt(verbose, lower = 0, upper = 3,
                         add = assertions)
  }

  if(missing(log_file)) {
    log_file <- ""
  } else {
    checkmate::assertPathForOutput(log_file, overwrite = TRUE,
                                   add = assertions)
    file.create(log_file)
  }

  if(missing(num_threads)) {
    # Determine number of CPUs within C++
    num_threads <- 0
  } else {
    checkmate::assertInt(num_threads, lower = 1,
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
    model_type <- 1 # classification
  } else {
    model_type <- 2 # regression
  }

  checkmate::reportAssertions(assertions)

  # Dependent variable name ----
  if(!missing(formula)) {
    dependent_variable_name <- names(data_selected)[1]
    independent_variable_names <- names(data_selected)[-1]
  } else {
    indepentend.variable.names <- colnames(data_selected)[colnames(data_selected) != dependent_variable_name]
  }

  # Input data and variable names, create final data matrix
  if(is.matrix(data_selected)) {
    data_final <- data_selected
  } else {
    data_final <- data.matrix(data_selected)[, -1]
  }
  variable_names <- colnames(data_final)

  # Clean up ----
  rm("data_selected")

  # Initialize output ----
  result <- list()

  # Internal cross validation ----
  if(!missing(folds) & !missing(cv_loss)) {
    checkmate::assertInt(folds, na.ok = TRUE, null.ok = TRUE, lower = 2, upper = 10,
                         add = assertions)
    checkmate::assertChoice(cv_loss, choices = c("auc", "bac"),
                            add = assertions)

    # Select prediction type
    pred_type <- switch(cv_loss,
                        "auc" = "prob",
                        "bac" = "response")

    # Select loss function
    measure <- switch(cv_loss,
                      "auc" = mlr::measureAUC,
                      "bac" = mlr::measureBAC)

    # Set search space for optimal top_results value
    if(length(top_results) == 1) {
      top_results <- 1:max(top_results)
    }

    fold_idx <- sample(1:folds, nrow(data_final), replace = TRUE)
    saved_mbmdr <- list()

    # Calculate the MB-MDR for each fold and assess current top_results value
    cv_performance <- rbindlist(lapply(1:folds, function(f) {
      # Generate MB-MDR models on CV training data
      data_cv_train <- data_final[fold_idx != f,]
      response_cv_train <- response[fold_idx != f]
      mbmdr <- mbmdrCpp(model_type = model_type, input_data = data_cv_train, response = response_cv_train,
                        order = order, alpha = alpha,
                        max_results = max_results,
                        num_threads = num_threads,
                        verbose = verbose,
                        log_file = log_file,
                        saved_mbmdr = saved_mbmdr)
      mbmdr <- mbmdr$mbmdr

      # Predict on CV testing data
      data_cv_test <- data_final[fold_idx == f,]
      pred <- predict.mbmdr(object = mbmdr, newdata = data_cv_test,
                            model_type = model_type,
                            top_results = max(top_results),
                            all = TRUE,
                            type = pred_type,
                            o_as_na = FALSE,
                            num_threads = num_threads,
                            verbose = verbose,
                            log_file = log_file)

      # Prepare CV loss for all top_result values
      pred <- melt(pred, id.vars = "ID",
                   variable.name = "top_results",
                   value.name = pred_type)
      pred[, list(cv_loss = measure(get(pred_type),
                                    response[fold_idx == f],
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

  # Call C++ ----
  result <- c(result, mbmdrCpp(model_type = model_type,
                               input_data = data_final, response = response,
                               order = order, alpha = alpha,
                               max_results = max_results,
                               num_threads = num_threads,
                               verbose = verbose,
                               log_file = log_file,
                               saved_mbmdr = saved_mbmdr))

  result$call <- sys.call()
  result$num_samples <- nrow(data_final)
  result$num_features <- ncol(data_final)
  result$num_combinations <- choose(ncol(data_final), order)
  result$model_type <- model_type
  result$top_results <- top_results

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
#' @param model_type    [\code{int}]\cr
#' 						          Code for type of MB-MDR models. 1 is for classification,
#' 						          2 is for regression.
#' @param top_results   [\code{int}]\cr
#'                      How many models are used for prediction.
#' @param o_as_na       [\code{bool}]\cr
#'                      Indicates how to handle samples in O classified genotype
#'                      combinations. See details.
#'                      Default is \code{TRUE}.
#' @param all           [\code{bool}]\cr
#' 						          Output predictions for all possible top results.
#' 						          Default is \code{FALSE}.
#' @param num_threads		[\code{int}]\cr
#' 									    Number of threads. Default is number of CPUs available.
#' @param verbose       [\code{int} or \code{string}]\cr
#' 						          Level of verbosity. Default is level 1 giving some basic
#' 						          information about progress. Level 0 will switch off any
#' 						          output. Levels 2 and 3 are for debugging purposes.
#' @param log_file		  [\code{string}]\cr
#'                      File path for a log file. If empty (default) logging is
#' 						          sent to console.
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
predict.mbmdr <- function(object, newdata, type = "response", model_type, top_results, o_as_na, all, num_threads, verbose, log_file, ...) {

  # hack to circumvent R CMD CHECK notes
  PROB <- NULL
  MODEL <- NULL
  ID <- NULL
  SCOREPROB <- NULL
  SCORE <- NULL

  if(missing(all)) {
    all <- FALSE
  } else {
    checkmate::assertFlag(all)
  }

  # Extract information ----
  response <- factor(sample(0:1, nrow(newdata), TRUE))
  order <- object$order
  alpha <- object$alpha
  max_results <- object$max_models
  top_results <- min(object$num_models, top_results)

  # Prepare data ----
  data_final <- data.matrix(newdata[, object$feature_names])

  # Call C++ ----
  predictions <- data.table::as.data.table(mbmdrCpp(model_type = model_type,
                                                    input_data = data_final, response = response,
                                                    order = order, alpha = alpha,
                                                    max_results = max_results,
                                                    num_threads = num_threads,
                                                    verbose = verbose,
                                                    log_file = log_file,
                                                    saved_mbmdr = object))
  setnames(predictions, "predictions", "PROB")

  if(!o_as_na) {
    # Convert to null hypothesis value
    predictions[is.na(PROB), PROB := 0.5]
  }

  # Add model and sample identifier
  predictions[, MODEL := rep(1:max_results, each = nrow(newdata))]
  predictions[, ID := rep(1:nrow(newdata), times = max_results)]

  if(all) {
    switch(type,
           # Round the mean case probability to 0 or 1 to return hard classification
           "response" = return(dcast(predictions[MODEL<=top_results, list(RESPONSE = round(cumsum(PROB)/1:max(MODEL)),
                                                                          TOPRESULTS = 1:max(MODEL)), by = c("ID")],
                                     ID~TOPRESULTS, value.var = "RESPONSE")),
           # Return mean case probability over all models for all top.results
           "prob" = return(dcast(predictions[MODEL<=top_results, list(PROB = cumsum(PROB)/1:max(MODEL),
                                                                      TOPRESULTS = 1:max(MODEL)), by = c("ID") ],
                                 ID~TOPRESULTS, value.var = "PROB")),
           # Return a risk score. Genotype combinations classified as H contribute +1,
           # genotype combinations classified as L contribute -1 and genotype combinations
           # classified as O contribute 0
           "score" = return(dcast(predictions[MODEL<=top_results, list(SCORE = cumsum(sign(PROB-0.5)),
                                                                       TOPRESULTS = 1:max(MODEL)), by = c("ID") ],
                                  ID~TOPRESULTS, value.var = "SCORE")),
           # Return the score transformed to a [0, 1] interval
           "scoreprob"= return(dcast(predictions[MODEL<=top_results, list(SCORE = cumsum(sign(PROB-0.5)),
                                                                          TOPRESULTS = 1:max(MODEL)), by = c("ID")][, SCOREPROB := range01(SCORE),
                                                                                                                    by = c("TOPRESULTS")],
                                     ID~TOPRESULTS,
                                     value.var = "SCOREPROB")))
  } else {
    switch(type,
           # Round the mean case probability to 0 or 1 to return hard classification
           "response" = return(predictions[MODEL<=top_results, list(predictions = round(mean(PROB, na.rm = TRUE))),
                                           by = c("ID")]$predictions),
           # Return mean case probability over all models
           "prob" = return(predictions[MODEL<=top_results, list(predictions = mean(PROB, na.rm = TRUE)),
                                       by = c("ID")]$predictions),
           # Return a risk score. Genotype combinations classified as H contribute +1,
           # genotype combinations classified as L contribute -1 and genotype combinations
           # classified as O contribute 0
           "score" = return(predictions[MODEL<=top_results, list(predictions = sum(sign(PROB-0.5), na.rm = TRUE)),
                                        by = c("ID")]$predictions),
           # Return the score transformed to a [0, 1] interval
           "scoreprob" = return(predictions[MODEL<=top_results, list(predictions = sum(sign(PROB-0.5), na.rm = TRUE)),
                                            by = c("ID")][, list(ID, predictions = range01(predictions))]$predictions))
  }

}

#' @rdname predict.mbmdr
#'
#' @export
predict.mbmdrc <- function(object, newdata, type = "response", top_results, o_as_na, num_threads, verbose, log_file, ...) {

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
    checkmate::assertInt(top_results, lower = 1, upper = object$mbmdr$num_models,
                         add = assertions)
    top_results <- min(object$mbmdr$num_models, top_results)
  }
  if(!missing(o_as_na)) {
    checkmate::assertFlag(o_as_na, add = assertions)
  } else {
    o_as_na = TRUE
  }
  if(!missing(num_threads)) {
    checkmate::assertInt(num_threads, lower = 0,
                         add = assertions)
  }
  if(missing(verbose)) {
    verbose <- 1
  } else {
    checkmate::assertInt(verbose, lower = 0, upper = 3,
                         add = assertions)
  }
  if(missing(log_file)) {
    log_file <- ""
  } else {
    checkmate::assertPathForOutput(log_file, overwrite = TRUE,
                                   add = assertions)
  }

  checkmate::reportAssertions(assertions)

  model_type <- object$model_type

  stats::predict(object$mbmdr, newdata = newdata,
                 model_type = model_type, type = type,
                 top_results = top_results,
                 o_as_na = o_as_na,
                 num_threads = num_threads,
                 verbose = verbose, log_file = log_file, ...)

}

range01 <- function(x, ...) {
  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

.onUnload <- function (libpath) {
  library.dynam.unload("MBMDRClassifieR", libpath)
}