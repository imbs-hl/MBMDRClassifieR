
#' Run the MB-MDR algorithm on the data
#'
#' Successively calculate the test statistics of all possible models, i.e.
#' feature combinations in the data.
#'
#' @param data          [\code{data.frame}]\cr
#'                      The data.
#' @param order         [\code{int}]\cr
#'                      The interaction depth, i.e. how many features are in one
#'                      model.
#' @param alpha         [\code{double}]\cr
#'                      Significance threshold for classifying a genotype combination
#'                      as (H)igh or (L)ow risk.
#' @param max.results   [\code{int}]\cr
#'                      Report only this amount of top models.
#'
#' @return A \code{list} object of class \code{mbmdr} containing the top models
#' and their respective test statistic, the HLO tables of the top models, the
#' model tables of the top models and the function arguments \code{order},
#' \code{alpha} and \code{max.results}.
mbmdr <- function(data, order = 2L, alpha = 0.1, max.results = 100) {

  # Count the cases
  cases <- sum(data$pheno == 1)

  # Total sample size
  N <- nrow(data)

  # Number of features
  k <- ncol(data) - 1

  # Case ratio under null hypothesis
  mu <- sum(cases)/N

  # Ensure integer type of interaction order
  order <- as.integer(order)

  # The initial model
  model <- 1L:order

  # Set max.results to the highest possible combination count
  max.results <- min(max.results, choose(k, order))

  # Initialize output
  result <- data.table::data.table(MODEL = list(vector("integer", order)),
                                   STATISTIC = rep(-Inf, max.results + 1))
  hlo_tables <- vector("list", max.results + 1)
  model_tables <- vector("list", max.results + 1)

  # Iterate through all possible models
  while (sum(model)) {

    # Get HLO table
    counts <- Step1(data = data, model = model, mu = mu, N = N, alpha = alpha)

    # Get model table
    model_counts <- Step2(counts = counts, mu = mu, N = N)

    # Create model name
    model_name <- colnames(data)[model + 1]

    # Insert model and statistic at the bottom
    # Use -99 to enable sorting if we are in the degenerative case of only NAs
    # as test statistics (which results in NaN and destroyes sorting)
    result[max.results + 1, c("MODEL", "STATISTIC") := list(list(model),
                                                            max(-99, model_counts[, max(S, na.rm = TRUE)]))]

    # Sort the results
    setorderv(result, "STATISTIC", order = -1L)

    # Save HLO table
    hlo_tables[[max.results + 1]] <- counts

    # Set model name
    names(hlo_tables)[max.results + 1] <- paste(model, collapse = ",")

    # Save model table
    model_tables[[max.results + 1]] <- model_counts

    # Set model name
    names(model_tables)[max.results + 1] <- paste(model, collapse = ",")

    # Reorder the lists as in the results list
    hlo_tables <- hlo_tables[sapply(result$MODEL, paste, collapse = ",")]
    model_tables <- model_tables[sapply(result$MODEL, paste, collapse = ",")]

    # Get the next model
    model <- NextModel(model, k = k, order = order)

  }

  # Replace strange value dummy with NA
  result[STATISTIC == -99, STATISTIC := NA]

  # Prepare output
  mbmdr <- list(result = result[1:max.results],
                hlo_tables = hlo_tables[1:max.results],
                model_tables = model_tables[1:max.results],
                order = order,
                alpa = alpha,
                max.results = max.results)

  # Set class
  class(mbmdr) <- "mbmdr"

  return(mbmdr)

}

#' Calculate HLO matrix
#'
#' For each model calculate the case ratio for each genotype combination
#' and test if this ratio is significantly different to the case ratio
#' in all other genotype combinations.
#'
#' @param data   [\code{data.frame}]\cr
#'               The data.
#' @param model  [\code{integer}]\cr
#'               The feature combination to consider.
#' @param mu     [\code{double}]\cr
#'               The overall case to control ratio.
#' @param N      [\code{int}]\cr
#'               The total sample size.
#' @param alpha  [\code{double}]\cr
#'               The significane level.
#'
#' @return A \code{data.table} object. Each row is a genotype combination of the
#' model. The columns contain the counts of cases and controls, the case ratio in
#' the genotype combination, the case ratio in all other genotype combinations,
#' the chi-square statistic, the corresponding p value and a class label.
Step1 <- function(data, model, mu, N, alpha) {

  tables <- data.table::melt(table(data[, c(model + 1, 1)]))

  counts <- data.table::as.data.table(data.table::dcast(tables, ...~pheno))
  data.table::setnames(counts, c("0", "1"), c("controls", "cases"))

  counts[, n1 := cases + controls]
  counts[, n0 := N - n1]

  counts[, y1 := cases/n1]
  counts[, y0 := (sum(cases) - cases)/n0]

  counts[, S := 1/(mu*(1-mu))*(n1*(y1-mu)^2+n0*(y0-mu)^2)]

  counts[, p := stats::pchisq(S, df = 1, lower.tail = FALSE)]

  counts[, label := "O"]
  counts[p<alpha & y1>y0, label := "H"]
  counts[p<alpha & y1<y0, label := "L"]

  return(counts)

}

#' Assess significance of model
#'
#' After classifying each genotype combination of a model as (H)igh, (L)ow or
#' uninformative (O), assess the significance of the model.
#'
#' @param counts  [\code{data.table}]\cr
#'                Output of \code{\link{Step1}}.
#' @param mu      [\code{double}]\cr
#'               The overall case to control ratio.
#' @param N      [\code{int}]\cr
#'               The total sample size.
#'
#' @return A \code{data.table} object, containing the counts of cases and controls
#' in each class of H, L or O, their respective ratios in and not in the class and
#' the chi-square test statistic.
Step2 <- function(counts, mu, N) {

  model_counts <- counts[, .(controls = sum(controls), cases = sum(cases)), by = label]

  model_counts[, n1 := controls + cases]
  model_counts[, n0 := N - n1]

  model_counts[, y1 := cases/n1]
  model_counts[, y0 := (sum(cases) - cases)/n0]

  model_counts[, S := 1/(mu*(1-mu))*(n1*(y1-mu)^2+n0*(y0-mu)^2)]

  return(model_counts)

}

#' Get the next successive model
#'
#' This function returns the next feature combination, i.e. model, given the
#' current model, the order and the total number of features.
#'
#' @param model [\code{integer}]\cr
#'              The current model.
#' @param k     [\code{int}]\cr
#'              Total number of features.
#' @param order [\code{int}]\cr
#'              Interaction order.
#' @param j     [\code{int}]\cr
#'              The model position to increase. If called, this is usually the
#'              interaction order.
#'
#' @return A \code{integer} vector of features. If the last model is reached,
#' \code{FALSE} is returned.
NextModel <- function(model, k, order, j = order) {

  # Input checks
  assertions <- checkmate::makeAssertCollection()
  checkmate::assertInt(k, add = assertions, lower = 1L)
  checkmate::assertInteger(model, add = assertions, lower = 0L, upper = k)
  checkmate::assertInt(order, add = assertions, lower = 1L)
  checkmate::assertInt(j, add = assertions, lower = 1L)
  checkmate::reportAssertions(assertions)

  # Check if last model is reached
  if (all(model == rev(k - (0L:(order - 1L))))) {

    return(FALSE)

  }

  if (model[j] < k - (order - j)) {
    # Current index to increase is less than maximum at this index
    model[j] <- model[j] + 1L

  } else {
    # Current index is at its maximum
    model <- NextModel(model, j = j - 1L, k = k, order = order)

    # Restart iteration in current index
    model[j] <- model[j - 1L] + 1L

  }

  return(model)

}