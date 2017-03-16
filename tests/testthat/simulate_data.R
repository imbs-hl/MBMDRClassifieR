
## ================================================================= ##
## Simulate SNP data ----
## ================================================================= ##
simulate_snps <- function(n.obs = 100, n.snp = 100, maf = 0.25) {

  assertions <- checkmate::makeAssertCollection()

  checkmate::assertCount(n.obs, na.ok = FALSE, null.ok = FALSE, positive = TRUE,
                         add = assertions)
  checkmate::assertCount(n.snp, na.ok = FALSE, null.ok = FALSE, positive = TRUE,
                         add = assertions)

  checkmate::assertNumeric(maf, lower = 0, upper = 0.5, min.len = 1, max.len = n.snp,
                           finite = TRUE,
                           any.missing = FALSE,
                           null.ok = FALSE,
                           add = assertions)

  checkmate::reportAssertions(assertions)

  if (length(maf) == 1) {
    maf <- rep(maf, n.snp)
  }

  ## Create SNP matrix
  X <- data.frame(sapply(maf, function(x) {
    sample(c(0, 1, 2), n.obs, replace = TRUE,
           prob = c((1-x)^2, 2*(1-x)*x, x^2))
  }))
  colnames(X) <- sprintf(paste0("SNP%0", nchar(n.snp), "d"), 1:n.snp)
  rownames(X) <- sprintf(paste0("ID%0", nchar(n.obs), "d"), 1:n.obs)

  return(X)

}

## ================================================================= ##
## Simulate response data ----
## ================================================================= ##
simulate_response <- function(X, beta0, beta = 0, type = "binary") {

  checkmate::assert(checkmate::checkClass(X, "matrix"),
                    checkmate::checkClass(X, "data.frame"))

  assertions <- checkmate::makeAssertCollection()

  checkmate::assertChoice(type, choices = c("binary", "quant", "surv"),
                          null.ok = FALSE)

  checkmate::assertNumeric(beta, min.len = 1, max.len = ncol(X) + 1,
                           finite = TRUE,
                           any.missing = FALSE,
                           null.ok = FALSE,
                           add = assertions)

  if (!missing(beta0)) {
    checkmate::assertNumeric(beta0, lower = 0, upper = 0.5, len = 1,
                             finite = TRUE,
                             any.missing = FALSE,
                             null.ok = FALSE,
                             add = assertions)
  }

  checkmate::reportAssertions(assertions)

  n.obs <- nrow(X)

  if (missing(beta0)) {
    beta0 <- NULL
  }
  if (length(beta) == 1) {
    beta <- rep(beta, n.snp)
  }

  ## Create response and return
  if (type == "binary") {
    z <- cbind(1, as.matrix(X)) %*% matrix(c(beta0, beta), ncol = 1)
    prob_pheno <- 1/(1+exp(-z))
    y <- factor(rbinom(n.obs, 1, prob_pheno))
  } else if (type == "quant") {
    y <- cbind(1, as.matrix(X)) %*% matrix(c(beta0, beta), ncol = 1) + rnorm(n.obs)
  } else if (type == "surv") {
    time <- cbind(1, as.matrix(X)) %*% matrix(c(beta0, beta), ncol = 1) + rnorm(n.obs)
    time <- time - min(time)
    status <- rbinom(n.obs, 1, 0.8)
    y <- data.frame(time = time, status = status)
  }

  return(y)

}

## ================================================================= ##
## Simulate data with two interacting SNPs ----
## ================================================================= ##
simulate_interaction_data <- function(n,
                                      type = "binary", prop_cases = 0.5,
                                      num_effect_snps, num_noise_snps,
                                      maf_interaction, maf_effect, maf_noise,
                                      model, beta, beta_effect) {

  assertions <- checkmate::makeAssertCollection()

  checkmate::assertCount(n, na.ok = FALSE, null.ok = FALSE, positive = TRUE,
                         add = assertions)

  checkmate::assertCount(num_effect_snps, na.ok = FALSE, null.ok = FALSE, positive = FALSE,
                         add = assertions)

  checkmate::assertCount(num_noise_snps, na.ok = FALSE, null.ok = FALSE, positive = FALSE,
                         add = assertions)

  checkmate::assertChoice(type, choices = c("binary", "quant", "surv"),
                          null.ok = FALSE)

  checkmate::assertChoice(model, choices = c("no_interaction", "synergistic",
                                             "interaction_only", "modifier_snp1",
                                             "modifier_snp2", "redundant"),
                          null.ok = FALSE, add = assertions)

  checkmate::assertNumeric(prop_cases, lower = 0, upper = 1, len = 1,
                           finite = type != "binary",
                           any.missing = type != "binary",
                           null.ok = type != "binary",
                           add = assertions)

  checkmate::assertNumeric(maf_interaction, lower = 0, upper = 0.5, len = 1,
                           finite = TRUE,
                           any.missing = FALSE,
                           null.ok = FALSE,
                           add = assertions)

  checkmate::assertNumeric(maf_effect, lower = 0, upper = 0.5, len = 1,
                           finite = TRUE,
                           any.missing = FALSE,
                           null.ok = FALSE,
                           add = assertions)

  if(!missing(maf_noise)) {
    checkmate::assertNumeric(maf_noise, lower = 0, upper = 0.5, len = 1,
                             finite = TRUE,
                             any.missing = FALSE,
                             null.ok = FALSE,
                             add = assertions)
  }

  checkmate::assertNumeric(beta, len = 1,
                           finite = TRUE,
                           any.missing = FALSE,
                           null.ok = FALSE,
                           add = assertions)

  checkmate::assertNumeric(beta_effect, len = 1,
                           finite = TRUE,
                           any.missing = FALSE,
                           null.ok = FALSE,
                           add = assertions)

  checkmate::reportAssertions(assertions)

  num_snps <- 2 + num_effect_snps + num_noise_snps

  ## Create model
  ## Intercept, Interaction SNP 1, Interaction SNP 2, Interaction, effect SNPs
  beta_model <- switch (model,
                        "no_interaction" = c(beta, beta, 0),
                        "synergistic" = c(beta, beta, beta),
                        "interaction_only" = c(0, 0, beta),
                        "modifier_snp1" = c(0, beta, beta),
                        "modifier_snp2" = c(beta, 0, beta),
                        "redundant" = c(beta, beta, -beta))

  ## Adjust beta_0 to get prop_cases as proportion of cases
  ## Additive model -> 2*maf^2
  if (type == "binary") {
    bx_interaction <- beta_model[1]*(2*maf_interaction^2 + 2*maf_interaction*(1-maf_interaction)) + # main effect of interaction SNP 1
      beta_model[2]*(2*maf_interaction^2 + 2*maf_interaction*(1-maf_interaction)) + # main effect of interaction SNP 2
      beta_model[3]*maf_interaction^2*(1+maf_interaction)^2 # interaction effect of interaction SNPs 1 and 2
    bx_effect <- num_effect_snps*beta_effect*(2*maf_effect^2 + 2*maf_effect*(1-maf_effect))
    beta_0 <- -log((1-prop_cases)/prop_cases) - bx_interaction - bx_effect
  } else {
    beta_0 <- 0
  }

  beta_model <- c(beta_0, beta_model, rep(beta_effect, num_effect_snps))

  ## Simulate interaction SNPs
  interaction_snps <- as.matrix(simulate_snps(n.obs = n, n.snp = 2,
                                              maf = maf_interaction))
  colnames(interaction_snps) <- sprintf("interactionSNP%d", 1:2)

  ## Simulate effect SNPs
  if (!missing(num_effect_snps) & num_effect_snps > 0) {
    effect_snps <- as.matrix(simulate_snps(n.obs = n, n.snp = num_effect_snps,
                                           maf = maf_effect))
    colnames(effect_snps) <- sprintf(paste0("effectSNP%0", nchar(num_effect_snps), "d"), 1:num_effect_snps)
  } else {
    effect_snps <- NULL
  }

  ## Simulate noise SNPs
  if (missing(maf_noise)) {
    maf_noise <- runif(num_noise_snps, max = 0.5)
  }
  if (!missing(num_noise_snps) & num_noise_snps > 0) {
    noise_snps <- as.matrix(simulate_snps(n.obs = n, n.snp = num_noise_snps, maf = maf_noise))
    colnames(noise_snps) <- sprintf(paste0("noiseSNP%0", nchar(num_noise_snps), "d"), 1:num_noise_snps)
  } else {
    noise_snps <- NULL
  }

  snps <- cbind(interaction_snps, "interaction" = apply(interaction_snps, 1, prod), effect_snps, noise_snps)

  ## Simulate phenotype, 0,1,2 coding -> additive model
  pheno <- simulate_response(snps[, 1:(num_effect_snps + 3)],
                             beta = beta_model, type = type)

  ## Create and return dataset
  data.frame(pheno, cbind(interaction_snps, effect_snps, noise_snps))
}
