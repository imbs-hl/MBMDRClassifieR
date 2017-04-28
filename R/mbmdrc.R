
#' MB-MDR based classifier
#'
#' Use the MB-MDR HLO classification of top models to predict the disease risk.
#'
#' @param formula       			[\code{formula}]\cr
#'                      			A formula of the form LHS ~ RHS, where LHS is the disease
#'                      			status and RHS are the features to enter the MB-MDR
#'                      			genotype classification.
#' @param data          			[\code{data.frame}]\cr
#'                      			The data object.
#' @param order         			[\code{int}]\cr
#'                      			Single integer specifying the interaction depth used
#'                      			in the MB-MDR algorithm.
#' @param alpha         			[\code{double}]\cr
#'                      			Single numeric as significance level used during HLO
#'                      			classification of genotype combinations.
#' @param max_results   			[\code{int}]\cr
#'                      			Single integer specifying the number of top models to
#'                      			report.
#' @param top_results   			[\code{int} or \code{integer}]\cr
#'                      			Single integer specifying how many models shall be used
#'                      			for prediction. If \code{folds} and \code{cv.loss} are
#'                      			set, this value specifies the upper limit of top results
#'                      			to assess in CV, or, if a \code{integer} vector the top
#'                      			result values to assess in CV.
#' @param folds         			[\code{int}]\cr
#'                      			Single interger specifying how many folds for internal
#'                      			cross validation to find optimal \code{top_results}
#'                      			should be used.
#' @param cv_loss       			[\code{character}]\cr
#'                      			One of \code{auc} or \code{bac}, specifying which loss
#'                      			should be used to find optimal \code{top_results}.
#' @param dependent.variable.name	[\code{string}]\cr
#' 									Name of dependent variable, needed if no formula given.
#' @param num_threads				[\code{int}]\cr
#' 									Number of threads. Default is number of CPUs available.
#' @param verbose                   [\code{int}]\cr
#' 									Level of verbosity. Default is level 1 giving some basic
#' 									information about progress. Level 0 will switch off any
#' 									output. Levels 2 and 3 are for debugging purposes.
#' 									Level 3 disables multithreading.
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
		verbose) {
	
	# Input checks
	assertions <- checkmate::makeAssertCollection()
	
	checkmate::assertDataFrame(data, min.cols = 2, min.rows = 2,
			add = assertions)
	checkmate::assertInt(order, lower = 2, upper = 5,
			add = assertions)
	checkmate::assertNumber(alpha, lower = 0, upper = 1,
			add = assertions)
	checkmate::assertInt(max_results, lower = 1, upper = 1e6)
	checkmate::assertInt(top_results, lower = 1, upper = max_results)
	checkmate::assertInt(folds, na.ok = TRUE, null.ok = TRUE, lower = 2, upper = 10,
			add = assertions)
	checkmate::assertChoice(cv_loss, choices = c("auc", "bac"),
			add = assertions)
	
	if(missing(verbose)) {
		verbose <- 1
	} else {
		checkmate::assertInt(verbose, lower = 0, upper = 3,
				add = assertions)
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
		formula <- as.formula(formula)
		checkmate::assertClass(formula, classes = "formula",
				add = assertions)
		data_selected <- model.frame(formula, data, na.action = na.pass)
		response <- data_selected[[1]]
	}
	
	# Prediction type ----
	if(is.factor(response)) {
		checkmate::assertFactor(response, n.levels = 2,
				add = assertions)
		pred_type <- 1 # classification
	} else {
		pred_type <- 2 # regression
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
	
	# Call C++ ----
	saved_mbmdr <- list()
	result <- mbmdrCpp(pred_type,
			data_final, response,
			order, alpha,
			max_results, top_results,
			num_threads,
			verbose,
			saved_mbmdr)
	
	result$call <- sys.call()
	result$num_samples <- nrow(data_final)
	result$num_features <- ncol(data_final)
	result$num_combinations <- choose(ncol(data_final), order)
	result$pred_type <- pred_type
	
	class(result) <- "mbmdrc"
	
	return(result)
	
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
#' @param top_results   [\code{int}]\cr
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
#'
#' @export
predict.mbmdrc <- function(object, newdata, type = "response", top_results, num_threads, verbose, ...) {
	
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
		checkmate::assertInt(top_results, lower = 1,
				add = assertions)
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
	
	checkmate::reportAssertions(assertions)
	
	# Extract information ----
	mbmdr <- object$mbmdr
	pred_type <- object$pred_type
	response <- factor(sample(0:1, n, TRUE))
	order <- mbmdr$order
	alpha <- mbmdr$alpha
	max_results <- mbmdr$max_models
	
	# Prepare data ----
	data_final <- data.matrix(newdata[, mbmdr$feature_names])
	
	# Call C++ ----
	predictions <- mbmdrCpp(pred_type,
			data_final, response,
			order, alpha,
			max_results, top_results,
			num_threads,
			verbose,
			mbmdr)
	
	return(predictions)
	
}

range01 <- function(x, ...) {
	(x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

.onUnload <- function (libpath) {
	library.dynam.unload("MBMDRClassifieR", libpath)
}