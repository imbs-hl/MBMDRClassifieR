
#' MB-MDR based classifier
#'
#' Use the MB-MDR HLO classification of top models to predict the disease risk.
#'
#' @param formula       [\code{formula}]\cr
#'                      A formula of the form LHS ~ RHS, where LHS is the disease
#'                      status and RHS are the features to enter the MB-MDR
#'                      genotype classification.
#' @param data          [\code{data.frame}]\cr
#'                      The data object.
#' @param order         [\code{int}]\cr
#'                      Single integer specifying the interaction depth used
#'                      in the MB-MDR algorithm.
#' @param alpha         [\code{double}]\cr
#'                      Single numeric as significance level used during HLO
#'                      classification of genotype combinations.
#' @param max.results   [\code{int}]\cr
#'                      Single integer specifying the number of top models to
#'                      report.
#' @param top.results   [\code{int} or \code{integer}]\cr
#'                      Single integer specifying how many models shall be used
#'                      for prediction. If \code{folds} and \code{cv.loss} are
#'                      set, this value specifies the upper limit of top results
#'                      to assess in CV, or, if a \code{integer} vector the top
#'                      result values to assess in CV.
#' @param folds         [\code{int}]\cr
#'                      Single interger specifying how many folds for internal
#'                      cross validation to find optimal \code{top.results}
#'                      should be used.
#' @param cv.loss       [\code{character}]\cr
#'                      One of \code{auc} or \code{bac}, specifying which loss
#'                      should be used to find optimal \code{top.results}.
#'
#' @return A S3 object of class \code{mbmdrc}.
#'
#' @export
#' @import data.table
mbmdrc <- function(formula, data, order = 2, alpha = 0.1, max.results = 100,
                   top.results = 10, folds, cv.loss) {

  # Initialize output
  mbmdrc <- list()

  # Create model frame
  mf <- stats::model.frame(formula, data)

  # Extract model terms
  mt <- attr(mf, "terms")

  # Replace data object with the RHS columns
  data <- mf[, attr(mt, "term.labels")]

  # Recode the outcome as factor with 0/1 labelling
  y <- factor(mf[, 1], labels = c(0, 1))

  # Bind data and outcome
  data <- cbind(pheno = y, data)

  # Replace top.results with highest value possible based on order and number of
  # columns
  top.results <- min(choose(ncol(data) -1 , order), top.results)

  if(!missing(folds) & !missing(cv.loss)) {
    # Enter CV to determine optimal top.results

    # Select prediction type
    pred_type <- switch(cv.loss,
                        "auc" = "prob",
                        "bac" = "response")

    # Select loss function
    measure <- switch(cv.loss,
                      "auc" = mlr::measureAUC,
                      "bac" = mlr::measureBAC)

    # Set search space for optimal top.results value
    if(length(top.results) == 1) {
      top_results <- 1:top.results
    } else {
      top_results <- top.results
    }

    # Perform CV for each possible value for top.results
    cv_performance <- sapply(top_results, function(tr) {

      # Split the data
      data$fold <- sample(1:folds, nrow(data), replace = TRUE)

      # Calculate the MB-MDR for each fold and assess current top.results value
      sapply(1:folds, function(f) {

        # Generate MB-MDR models on training data
        mbmdr <- mbmdr(data[data$fold != f, -ncol(data)], order, alpha, max.results)

        # Predict the status for current top.results value on leave out data
        pred <- predict(mbmdr, data[data$fold == f, -ncol(data)], top.results = tr, type = pred_type)

        # Calculate CV loss
        measure(pred$predictions, truth = data[data$fold == f, "pheno"],
                positive = 1, negative = 0)

      })

    })

    # Select top.results value with optimal loss
    top.results <- top_results[which.max(mean_cv_performance)]

    # Save CV results
    mbmdrc$cv_performance <- cv_performance
    mbmdrc$mean_cv_performance <- mean(cv_performance)

  }

  # Train MB-MDR on the full dataset
  mbmdr <- mbmdr(data = data,
                 order = order,
                 alpha = alpha,
                 max.results = max.results)

  # Save mbmdr object
  mbmdrc$mbmdr <- mbmdr

  # Save top.results value
  mbmdrc$top.results <- top.results

  # Save function call
  mbmdrc$call <- sys.call()


  class(mbmdrc) <- "mbmdrc"

  return(mbmdrc)

}

#' MB-MDR classifier prediction
#'
#' Prediction with new data and saved MB-MDR classifier object.
#'
#' @param object        [\code{mbmdr}]\cr
#'                      MB-MDR models and HLO tables as output from \code{\link{mbmdr}}.
#' @param newdata       [\code{newdata}]\cr
#'                      New data to predict class status for.
#' @param type          [\code{string}]\cr
#'                      Type of prediction. One of \code{response}, \code{prob},
#'                      \code{score} or \code{scoreprob}. See details.
#'                      Default is \code{response}.
#' @param o.as.na       [\code{bool}]\cr
#'                      Indicates how to handle samples in O classified genotype
#'                      combinations. See details.
#'                      Default is \code{TRUE}.
#' @param top.results   [\code{int}]\cr
#'                      How many models are used for prediction.
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
predict.mbmdr <- function(object, newdata, type = "response", o.as.na = TRUE, top.results, ...) {

  # Convert data to a data.table object, keep rownames
  data <- data.table::as.data.table(newdata, keep.rownames = "ID")

  # Ensure correct order of MB-MDR models
  setorderv(object$result, "STATISTIC", order = -1L)

  # Get case probability for all sample genotypes for the first top.results models
  model_names <- sapply(object$result[1:top.results, MODEL],
                        function(model) {
                          # Extract HLO table
                          hlo_table <- object$hlo_tables[[paste(model, collapse = ",")]]

                          # Create model name from feature names
                          model_name <- paste(names(hlo_table)[1:object$order], collapse = ".")

                          # Merge the data with the HLO table
                          data <<- merge(data, hlo_table,
                                         by = names(hlo_table)[1:object$order],
                                         sort = FALSE, all.x = TRUE)

                          # Add a new column containing the case probabilities
                          # for all samples
                          data[, c(model_name) := list(y1)]

                          # Set case probability to NA if genotype combination
                          # is classified as O (not significant)
                          data[label == "O", c(model_name) := list(NA)]

                          # Remove obsolete columns
                          data[, names(hlo_table)[-(1:object$order)] := list(NULL)]

                          return(model_name)
                        })

  # Transform data to long format, drop original genotypes, keep case probabilities
  data <- melt(data, id.vars = c("ID", "pheno"),
               measure.vars = model_names,
               variable.name = "MODEL",
               value.name = "PROB")

  # If requested by user, replace case probabilities in O genotype combinations
  # with 0.5, which is the null hypothesis
  if(!o.as.na) {
    data[is.na(PROB), PROB := 0.5]
  }

  if(type == "response") {
    # Round the mean case probability to 0 or 1 to return hard classification
    return(data[, .(predictions = round(mean(PROB, na.rm = TRUE))),
                by = c("ID")])
  }
  if(type == "prob") {
    # Return mean case probability over all models
    return(data[, .(predictions = mean(PROB, na.rm = TRUE)),
                by = c("ID")])
  }
  if(type == "score") {
    # Return a risk score. Genotype combinations classified as H contribute +1,
    # genotype combinations classified as L contribute -1 and genotype combinations
    # classified as O contribute 0
    return(data[, .(predictions = sum(round(PROB)*2-1, na.rm = TRUE)),
                by = c("ID")])
  }
  if(type == "scoreprob") {
    # Return the score transformed to a [0, 1] interval
    return(data[, .(predictions = range01(sum(round(PROB)*2-1, na.rm = TRUE))),
                by = c("ID")])
  }

}

#' @rdname predict.mbmdr
#'
#' @export
predict.mbmdrc <- function(object, newdata, type = "response", o.as.na = TRUE, top.results = object$top.results, ...) {

  top.results <- min(nrow(object$result), top.results)

  predict(object$mbmdr, newdata = newdata, type = type, na.O = na.O, top.results = top.results)

}

range01 <- function(x, ...) {
  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}