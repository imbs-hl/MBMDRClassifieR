#' Calculate AUC
#'
#' @description Computes the Area under the Receiver Operating Curve.
#'
#' @param probabilities [\code{numeric}]\cr
#'                      The predicted probabilities for the positive class as a
#'                      numeric vector.
#' @param truth         [\code{factor} or \code{numeric}]\cr
#'                      Vector of the true class.
#' @param negative      [\code{character}]\cr
#'                      The name of the negative class.
#' @param positive      [\code{character}]\cr
#'                      The name of the positive class
#'
#' @return The Area under the Receiver Operating Curve as single \code{numeric}.
#'
auc <- function(probabilities, truth, negative, positive) {
  if (is.factor(truth)) {
    i = as.integer(truth) == which(levels(truth) == positive)
  }
  else {
    i = truth == positive
  }
  if (length(unique(i)) < 2L) {
    stop("truth vector must have at least two classes")
  }
  if (length(i) > 5000L) {
    r = data.table::frankv(probabilities)
  }
  else {
    r = rank(probabilities)
  }
  n_pos = as.numeric(sum(i))
  n_neg = length(i) - n_pos
  (sum(r[i]) - n_pos * (n_pos + 1)/2)/(n_pos * n_neg)
}

#' Calculate BAC
#'
#' @description Computes the balanced accuracy.
#'
#' @param response      [\code{factor}]\cr
#'                      Vector of the predicted class.
#' @param truth         [\code{factor} or \code{numeric}]\cr
#'                      Vector of the true class.
#' @param negative      [\code{character}]\cr
#'                      The name of the negative class.
#' @param positive      [\code{character}]\cr
#'                      The name of the positive class
#'
#' @return The balanced accuracy as single \code{numeric}.
#'
bac <- function(response, truth, negative, positive) {
  mean(c(tp(response, truth, positive)/sum(truth == positive),
         tn(response, truth, negative)/sum(truth == negative)))
}

#' Calculate TP
#'
#' @description Computes the number of true positives.
#'
#' @param response      [\code{factor}]\cr
#'                      Vector of the predicted class.
#' @param truth         [\code{factor} or \code{numeric}]\cr
#'                      Vector of the true class.
#' @param positive      [\code{character}]\cr
#'                      The name of the positive class.
#'
#' @return The number of true positives as single \code{numeric}.
#'
tp <- function(response, truth, positive) {
  sum(truth == response & response == positive)
}

#' Calculate TN
#'
#' @description Computes the number of true negatives
#'
#' @param response      [\code{factor}]\cr
#'                      Vector of the predicted class.
#' @param truth         [\code{factor} or \code{numeric}]\cr
#'                      Vector of the true class.
#' @param negative      [\code{character}]\cr
#'                      The name of the negative class.
#'
#' @return The number of true negatives as single \code{numeric}.
#'
tn <- function(response, truth, negative) {
  sum(truth == response & response == negative)
}