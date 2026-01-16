#' Multidimensional Mode-Based Estimation for Mendelian Randomization
#'
#' Performs mode-based estimation (MBE) for multiple correlated outcomes simultaneously
#' using multidimensional kernel density estimation. This approach identifies the modal
#' causal effect across multiple traits while accounting for pleiotropy.
#'
#' @param BetaXG Numeric vector of SNP-exposure associations
#' @param BetaYG_matrix Matrix of SNP-outcome associations (rows = SNPs, columns = outcomes/traits)
#' @param seBetaXG Numeric vector of standard errors for BetaXG
#' @param seBetaYG_matrix Matrix of standard errors for BetaYG_matrix
#' @param phi Numeric vector of bandwidth multipliers for sensitivity analysis (default: c(1, 0.5, 0.25))
#' @param n_boot Integer, number of bootstrap iterations for standard error estimation (default: 10000)
#' @param alpha Numeric, significance level for confidence intervals (default: 0.05)
#'
#' @return A data frame containing:
#' \item{Method}{Method name}
#' \item{phi}{Bandwidth multiplier used}
#' \item{Estimate.1, Estimate.2, ...}{Causal effect estimates for each outcome}
#' \item{SE.1, SE.2, ...}{Standard errors for each outcome}
#' \item{CI_low.1, CI_low.2, ...}{Lower confidence interval bounds}
#' \item{CI_upp.1, CI_upp.2, ...}{Upper confidence interval bounds}
#' \item{P.1, P.2, ...}{P-values for each outcome}
#'
#' @details
#' The function uses weighted multidimensional kernel density estimation to find the
#' mode of the joint distribution of causal effects across multiple outcomes. Bootstrap
#' resampling is used to estimate standard errors. Multiple bandwidth parameters (phi)
#' can be tested for sensitivity analysis.
#'
#' @importFrom ks kde
#' @importFrom stats rnorm qnorm pt mad sd
#' @export
#'
IBMODE <- function(BetaXG,
                           BetaYG_matrix,
                           seBetaXG,
                           seBetaYG_matrix,
                           phi = c(1, 0.5, 0.25),
                           n_boot = 1e4,
                           alpha = 0.05) {

  # Check if ks package is available
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Package 'ks' is required. Please install it with: install.packages('ks')")
  }

  # Input validation
  if (!is.matrix(BetaYG_matrix)) {
    BetaYG_matrix <- as.matrix(BetaYG_matrix)
  }
  if (!is.matrix(seBetaYG_matrix)) {
    seBetaYG_matrix <- as.matrix(seBetaYG_matrix)
  }

  if (length(BetaXG) != nrow(BetaYG_matrix)) {
    stop("BetaXG length must equal number of rows in BetaYG_matrix")
  }
  if (length(seBetaXG) != nrow(seBetaYG_matrix)) {
    stop("seBetaXG length must equal number of rows in seBetaYG_matrix")
  }
  if (!all(dim(BetaYG_matrix) == dim(seBetaYG_matrix))) {
    stop("BetaYG_matrix and seBetaYG_matrix must have the same dimensions")
  }

  n_outcomes <- ncol(BetaYG_matrix)
  n_snps <- nrow(BetaYG_matrix)

  #--------------------------------------#
  # Function to compute the point estimate
  #--------------------------------------#
  compute_mode_estimate <- function(BetaIV_matrix, seBetaIV_matrix) {

    # Compute bandwidth matrix (diagonal)
    S <- diag(n_outcomes)
    for (i in 1:n_outcomes) {
      S[i, i] <- sd(BetaIV_matrix[, i]) / n_snps^(1/6)
    }

    # Compute weights based on inverse product of standard errors
    weights <- 1 / apply(seBetaIV_matrix, 1, prod)
    weights <- weights / sum(weights)

    # Store estimates for each phi
    estimates <- matrix(NA, nrow = length(phi), ncol = n_outcomes)

    for (idx in seq_along(phi)) {
      cur_phi <- phi[idx]
      H <- cur_phi * S %*% t(S)  # Bandwidth matrix

      # Compute kernel density estimate
      kde_result <- ks::kde(x = BetaIV_matrix, H = H, compute.cont = TRUE, w = weights)

      # Find mode (maximum density)
      mode_index <- which(kde_result$estimate == max(kde_result$estimate), arr.ind = TRUE)

      # Extract mode coordinates
      for (j in 1:n_outcomes) {
        estimates[idx, j] <- kde_result$eval.points[[j]][mode_index[j]]
      }
    }

    return(estimates)
  }

  #------------------------------------------#
  # Function to estimate SEs through bootstrap
  #------------------------------------------#
  bootstrap_estimates <- function(BetaIV_matrix, seBetaIV_matrix) {

    beta_boot <- array(NA, dim = c(n_boot, length(phi), n_outcomes))

    for (i in 1:n_boot) {
      # Bootstrap: sample from normal distribution for each SNP-outcome pair
      BetaIV_boot <- matrix(NA, nrow = n_snps, ncol = n_outcomes)
      for (j in 1:n_outcomes) {
        BetaIV_boot[, j] <- rnorm(n_snps,
                                  mean = BetaIV_matrix[, j],
                                  sd = seBetaIV_matrix[, j])
      }

      # Compute mode estimate for bootstrap sample
      beta_boot[i, , ] <- compute_mode_estimate(BetaIV_boot, seBetaIV_matrix)
    }

    return(beta_boot)
  }

  #--------------------------------------#
  # Main computation
  #--------------------------------------#

  # Compute ratio estimates (Wald ratios)
  BetaIV_matrix <- sweep(BetaYG_matrix, 1, BetaXG, FUN = "/")

  # Compute standard errors using delta method
  seBetaIV_matrix <- sqrt(
    (seBetaYG_matrix^2) / (BetaXG^2) +
      (BetaYG_matrix^2 * seBetaXG^2) / (BetaXG^4)
  )

  # Point estimates
  beta_MBE <- compute_mode_estimate(BetaIV_matrix, seBetaIV_matrix)

  # Bootstrap for standard errors
  message(paste0("Running ", n_boot, " bootstrap iterations..."))
  beta_MBE_boot <- bootstrap_estimates(BetaIV_matrix, seBetaIV_matrix)

  # Compute statistics for each outcome
  se_MBE <- matrix(NA, nrow = length(phi), ncol = n_outcomes)
  CIlow_MBE <- matrix(NA, nrow = length(phi), ncol = n_outcomes)
  CIupp_MBE <- matrix(NA, nrow = length(phi), ncol = n_outcomes)
  P_MBE <- matrix(NA, nrow = length(phi), ncol = n_outcomes)

  for (i in 1:length(phi)) {
    for (j in 1:n_outcomes) {
      # Standard error (using MAD for robustness)
      se_MBE[i, j] <- mad(beta_MBE_boot[, i, j])

      # Confidence intervals
      CIlow_MBE[i, j] <- beta_MBE[i, j] - qnorm(1 - alpha / 2) * se_MBE[i, j]
      CIupp_MBE[i, j] <- beta_MBE[i, j] + qnorm(1 - alpha / 2) * se_MBE[i, j]

      # P-value (two-tailed t-test)
      t_stat <- beta_MBE[i, j] / se_MBE[i, j]
      P_MBE[i, j] <- 2 * pt(abs(t_stat), df = n_snps - 1, lower.tail = FALSE)
    }
  }

  #--------------------------------------#
  # Format results
  #--------------------------------------#

  # Create column names for outcomes
  outcome_names <- if (!is.null(colnames(BetaYG_matrix))) {
    colnames(BetaYG_matrix)
  } else {
    paste0("Outcome", 1:n_outcomes)
  }

  # Build results data frame
  Results <- data.frame(
    Method = rep("Multidimensional MBE", length(phi)),
    phi = phi,
    stringsAsFactors = FALSE
  )

  # Add estimates for each outcome
  for (j in 1:n_outcomes) {
    Results[[paste0("Estimate_", outcome_names[j])]] <- beta_MBE[, j]
    Results[[paste0("SE_", outcome_names[j])]] <- se_MBE[, j]
    Results[[paste0("CI_low_", outcome_names[j])]] <- CIlow_MBE[, j]
    Results[[paste0("CI_upp_", outcome_names[j])]] <- CIupp_MBE[, j]
    Results[[paste0("P_", outcome_names[j])]] <- P_MBE[, j]
  }

  return(Results)
}
