#' Calculate Coheterogeneity Q Statistics
#'
#' Computes coheterogeneity Q and Q-correlation matrices for multiple traits
#' in Mendelian Randomization analyses. This quantifies the similarity of
#' causal effect estimates across different outcome traits.
#'
#' @param BetaXG Numeric vector of SNP-exposure associations (optional if BetaIV_matrix provided)
#' @param BetaYG_matrix Matrix of SNP-outcome associations, rows = SNPs, columns = traits (optional if BetaIV_matrix provided)
#' @param seBetaXG Numeric vector of standard errors for BetaXG (optional)
#' @param seBetaYG_matrix Matrix of standard errors for BetaYG_matrix (optional)
#' @param BetaIV_matrix Matrix of instrumental variable estimates (BetaYG/BetaXG), rows = SNPs, columns = traits
#' @param seBetaIV_matrix Matrix of standard errors for BetaIV_matrix (optional)
#'
#' @return A list containing:
#' \item{Q_matrix}{Matrix of coheterogeneity Q statistics between trait pairs}
#' \item{Q_corr_matrix}{Matrix of Q-correlation values between trait pairs (normalized Q)}
#'
#' @export
coheterogeneity_Q <- function(BetaXG = NULL,
                              BetaYG_matrix = NULL,
                              seBetaXG = NULL,
                              seBetaYG_matrix = NULL,
                              BetaIV_matrix = NULL,
                              seBetaIV_matrix = NULL) {

  # Input validation: Check if BetaIV_matrix is provided; otherwise, compute it
  if (is.null(BetaIV_matrix)) {
    if (is.null(BetaXG) || is.null(BetaYG_matrix)) {
      stop("Either provide BetaIV_matrix or both BetaXG and BetaYG_matrix.")
    }

    # Convert to matrix if needed
    if (!is.matrix(BetaYG_matrix)) {
      BetaYG_matrix <- as.matrix(BetaYG_matrix)
    }

    # Dimension checks
    if (nrow(BetaYG_matrix) != length(BetaXG)) {
      stop("BetaYG_matrix must have the same number of rows as the length of BetaXG.")
    }

    # Compute IV estimates (Wald ratio)
    BetaIV_matrix <- sweep(BetaYG_matrix, 1, BetaXG, FUN = "/")
  } else {
    # Ensure BetaIV_matrix is a matrix
    if (!is.matrix(BetaIV_matrix)) {
      BetaIV_matrix <- as.matrix(BetaIV_matrix)
    }
  }

  # Compute seBetaIV_matrix if not provided but component standard errors exist
  if (is.null(seBetaIV_matrix)) {
    if (!is.null(seBetaXG) && !is.null(seBetaYG_matrix)) {
      if (!is.matrix(seBetaYG_matrix)) {
        seBetaYG_matrix <- as.matrix(seBetaYG_matrix)
      }

      if (nrow(seBetaYG_matrix) != length(seBetaXG)) {
        stop("seBetaYG_matrix must have the same number of rows as the length of seBetaXG.")
      }

      # Delta method for standard error of ratio
      seBetaIV_matrix <- sqrt(
        (seBetaYG_matrix^2) / (BetaXG^2) +
          (BetaYG_matrix^2 * seBetaXG^2) / (BetaXG^4)
      )
    } else {
      message("Standard errors not provided. Using uniform weights for Q statistics.")
    }
  } else {
    if (!is.matrix(seBetaIV_matrix)) {
      seBetaIV_matrix <- as.matrix(seBetaIV_matrix)
    }
  }

  # Get number of traits
  num_traits <- ncol(BetaIV_matrix)
  num_snps <- nrow(BetaIV_matrix)

  # Initialize output matrices
  Q_matrix <- matrix(NA, nrow = num_traits, ncol = num_traits)
  Q_corr_matrix <- matrix(NA, nrow = num_traits, ncol = num_traits)

  # Add trait names if available
  if (!is.null(colnames(BetaIV_matrix))) {
    rownames(Q_matrix) <- colnames(Q_matrix) <- colnames(BetaIV_matrix)
    rownames(Q_corr_matrix) <- colnames(Q_corr_matrix) <- colnames(BetaIV_matrix)
  }

  # Compute coheterogeneity Q and Q-correlation for all trait pairs
  for (i in 1:num_traits) {
    for (j in 1:num_traits) {
      if (i == j) {
        # Diagonal: self-comparison
        Q_matrix[i, j] <- NA
        Q_corr_matrix[i, j] <- 1
      } else {
        # Extract beta estimates for trait pair
        beta1 <- BetaIV_matrix[, i]
        beta2 <- BetaIV_matrix[, j]

        # Remove missing values
        valid_idx <- complete.cases(beta1, beta2)

        if (sum(valid_idx) < 2) {
          warning(paste0("Insufficient valid data for trait pair (", i, ", ", j, "). Skipping."))
          next
        }

        beta1 <- beta1[valid_idx]
        beta2 <- beta2[valid_idx]

        # Calculate means
        beta1_mean <- mean(beta1)
        beta2_mean <- mean(beta2)

        # Determine weights
        if (!is.null(seBetaIV_matrix)) {
          se1 <- seBetaIV_matrix[valid_idx, i]
          se2 <- seBetaIV_matrix[valid_idx, j]

          # Avoid division by zero
          valid_se_idx <- (se1 > 0) & (se2 > 0) & is.finite(se1) & is.finite(se2)

          if (sum(valid_se_idx) < 2) {
            warning(paste0("Insufficient valid SE data for trait pair (", i, ", ", j, "). Using uniform weights."))
            weights <- rep(1, length(beta1))
          } else {
            beta1 <- beta1[valid_se_idx]
            beta2 <- beta2[valid_se_idx]
            se1 <- se1[valid_se_idx]
            se2 <- se2[valid_se_idx]

            weights <- 1 / (se1 * se2)
          }
        } else {
          weights <- rep(1, length(beta1))
        }

        # Compute coheterogeneity Q (weighted covariance)
        Q_matrix[i, j] <- sum(weights * (beta1 - beta1_mean) * (beta2 - beta2_mean))

        # Compute Q-correlation (normalized Q)
        sd_beta1 <- sqrt(sum(weights * (beta1 - beta1_mean)^2))
        sd_beta2 <- sqrt(sum(weights * (beta2 - beta2_mean)^2))

        if (sd_beta1 > 0 && sd_beta2 > 0) {
          Q_corr_matrix[i, j] <- Q_matrix[i, j] / (sd_beta1 * sd_beta2)
        } else {
          Q_corr_matrix[i, j] <- NA
        }
      }
    }
  }

  return(list(
    Q_matrix = Q_matrix,
    Q_corr_matrix = Q_corr_matrix,
    num_traits = num_traits,
    num_snps = num_snps
  ))
}
