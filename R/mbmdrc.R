
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
		checkmate::assertChoice(verbose, choices = 0:3,
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
			variable_names,
			order, alpha,
			max_results, top_results,
			num_threads,
			verbose,
			saved_mbmdr)
	
	result$call <- sys.call()
	result$num_samples <- nrow(data_final)
	result$num_features <- ncol(data_final)
	result$num_combinations <- choose(ncol(data_final), order)
	
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
predict.mbmdr <- function(object, newdata, type = "response", top_results, all = FALSE, ...) {
	
	# Convert data to a data.table object, keep rownames
	data <- data.table::as.data.table(newdata, keep.rownames = "ID")
	
	# Ensure correct order of MB-MDR models
	setorderv(object$result, "STATISTIC", order = -1L, na.last = TRUE)
	
	# Get case probability for all sample genotypes for the first top_results models
	model_names <- unlist(sapply(1:top_results,
					function(idx) {
						if(object$result[idx, is.na(STATISTIC)]) {
							return()
						}
						
						# Extract model
						model <- object$result[idx, MODEL][[1]]
						
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
					}))
	
	# Transform data to long format, drop original genotypes, keep case probabilities
	data <- melt(data, id.vars = c("ID"),
			measure.vars = model_names,
			variable.name = "MODEL",
			value.name = "PROB")
	
	# Replace case probabilities in O genotype combinations with 0.5, which is the
	# null hypothesis
	data[is.na(PROB), PROB := 0.5]
	
	if(all) {
		
		model_ranking <- data.table(MODEL = model_names, RANK = seq_along(model_names))
		data[model_ranking, RANK := i.RANK, on = "MODEL"]
		
		if(type == "response") {
			# Round the mean case probability to 0 or 1 to return hard classification
			return(dcast(data[, .(RESPONSE = round(cumsum(PROB)/1:max(RANK)),
											TOPRESULTS = 1:max(RANK)), by = c("ID") ],
							ID~TOPRESULTS, value.var = "RESPONSE"))
		}
		if(type == "prob") {
			# Return mean case probability over all models for all top_results
			return(dcast(data[, .(PROB = cumsum(PROB)/1:max(RANK),
											TOPRESULTS = 1:max(RANK)), by = c("ID") ],
							ID~TOPRESULTS, value.var = "PROB"))
		}
		if(type == "score") {
			# Return a risk score. Genotype combinations classified as H contribute +1,
			# genotype combinations classified as L contribute -1 and genotype combinations
			# classified as O contribute 0
			return(dcast(data[, .(SCORE = cumsum(round(PROB)*2-1),
											TOPRESULTS = 1:max(RANK)), by = c("ID") ],
							ID~TOPRESULTS, value.var = "SCORE"))
		}
		if(type == "scoreprob") {
			# Return the score transformed to a [0, 1] interval
			return(dcast(data[, .(SCORE = cumsum(round(PROB)*2-1),
											TOPRESULTS = 1:max(RANK)), by = c("ID")][, SCOREPROB := range01(SCORE),
									by = c("TOPRESULTS")],
							ID~TOPRESULTS,
							value.var = "SCOREPROB"))
		}
		
	} else {
		
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
			return(data[, .(predictions = sum(round(PROB)*2-1, na.rm = TRUE)),
							by = c("ID")][, .(ID, predictions = range01(predictions))])
		}
	}
	
}

#' @rdname predict.mbmdr
#'
#' @export
predict.mbmdrc <- function(object, newdata, type = "response", top_results = object$top_results, all = FALSE, ...) {
	
	top_results <- min(nrow(object$result), top_results)
	
	stats::predict(object$mbmdr, newdata = newdata, type = type, top_results = top_results, all = all)
	
}

range01 <- function(x, ...) {
	(x - min(x, ...)) / (max(x, ...) - min(x, ...))
}