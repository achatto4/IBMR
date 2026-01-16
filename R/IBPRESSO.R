#' MR-PRESSO with Instrument Borrowing
#'
#' Performs Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO)
#' test extended for instrument borrowing scenarios. Detects and corrects for
#' horizontal pleiotropy using an auxiliary trait to identify outlier instruments.
#'
#' @param BetaOutcome Character, column name for primary outcome SNP effects
#' @param BetaExposure Character or vector, column name(s) for exposure SNP effects
#' @param BetaAux Character, column name for auxiliary trait SNP effects (borrowed instrument)
#' @param SdOutcome Character, column name for standard errors of outcome effects
#' @param SdExposure Character or vector, column name(s) for standard errors of exposure effects
#' @param SdAux Character, column name for standard errors of auxiliary trait effects
#' @param data Data frame containing all the above columns
#' @param OUTLIERtest Logical, whether to perform outlier detection (default: FALSE)
#' @param DISTORTIONtest Logical, whether to test if outliers cause distortion (default: FALSE)
#' @param SignifThreshold Numeric, significance threshold for tests (default: 0.05)
#' @param NbDistribution Integer, number of null distributions to generate (default: 1000)
#' @param seed Integer, random seed for reproducibility (default: NULL)
#' @param n_cores Integer, number of cores for parallel processing (default: 1)
#'
#' @return A list containing:
#' \item{raw_beta}{Raw causal effect estimate before outlier correction}
#' \item{raw_se}{Standard error of raw estimate}
#' \item{corrected_beta}{Outlier-corrected causal effect estimate (if outliers detected)}
#' \item{corrected_se}{Standard error of corrected estimate}
#' \item{p_value}{Global test p-value for horizontal pleiotropy}
#' \item{p_distort}{P-value for distortion test (if DISTORTIONtest = TRUE)}
#' \item{outlier_idx}{Indices of detected outlier instruments}
#' \item{n_instruments}{Number of instruments used}
#' \item{n_outliers}{Number of outliers detected}
#'
#' @details
#' This function extends MR-PRESSO to handle instrument borrowing by incorporating
#' an auxiliary trait. It uses leave-one-out cross-validation and Mahalanobis
#' distance in a bivariate residual space (outcome + auxiliary trait) to detect
#' pleiotropic outliers. The global test assesses overall horizontal pleiotropy,
#' while the distortion test evaluates whether removing outliers significantly
#' changes the causal estimate.
#'
#' @importFrom MASS cov.rob
#' @importFrom stats lm coef var predict quantile mahalanobis rnorm complete.cases
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n_snps <- 50
#' dat <- data.frame(
#'   beta_exposure = rnorm(n_snps, 0.1, 0.02),
#'   beta_outcome = rnorm(n_snps, 0.05, 0.03),
#'   beta_aux = rnorm(n_snps, 0.04, 0.025),
#'   se_exposure = abs(rnorm(n_snps, 0.01, 0.002)),
#'   se_outcome = abs(rnorm(n_snps, 0.015, 0.003)),
#'   se_aux = abs(rnorm(n_snps, 0.012, 0.002))
#' )
#'
#' # Run MR-PRESSO with instrument borrowing
#' results <- mrpresso_ib(
#'   BetaOutcome = "beta_outcome",
#'   BetaExposure = "beta_exposure",
#'   BetaAux = "beta_aux",
#'   SdOutcome = "se_outcome",
#'   SdExposure = "se_exposure",
#'   SdAux = "se_aux",
#'   data = dat,
#'   OUTLIERtest = TRUE,
#'   DISTORTIONtest = TRUE,
#'   NbDistribution = 1000
#' )
#' }
IBPRESSO <- function(
    BetaOutcome,
    BetaExposure,
    BetaAux,
    SdOutcome,
    SdExposure,
    SdAux,
    data,
    OUTLIERtest = FALSE,
    DISTORTIONtest = FALSE,
    SignifThreshold = 0.05,
    NbDistribution = 1000,
    seed = NULL,
    n_cores = 1
) {

  # Check required packages
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Please install it with: install.packages('MASS')")
  }

  if (n_cores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for multi-core processing. Please install it with: install.packages('parallel')")
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Input validation
  required_cols <- c(BetaExposure, BetaOutcome, BetaAux, SdExposure, SdOutcome, SdAux)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }

  # Prepare core dataset
  core <- data[, required_cols, drop = FALSE]
  core <- core[complete.cases(core), ]

  if (nrow(core) == 0) {
    stop("No complete cases found in the data")
  }

  # Compute weights
  core$Weights <- 1 / (core[[SdOutcome]] * core[[SdAux]])
  core$Weights1 <- 1 / (core[[SdOutcome]]^2)

  n_instruments <- nrow(core)
  n_params <- length(BetaExposure)

  if (n_instruments <= n_params + 2) {
    stop(paste0("Insufficient instruments. Need at least ", n_params + 3,
                " instruments, but only ", n_instruments, " available."))
  }

  message(paste0("Running MR-PRESSO-IB with ", n_instruments, " instruments..."))

  #--------------------------------------#
  # Helper: Robust covariance estimation
  #--------------------------------------#
  compute_robust_cov <- function(mat) {
    covmat <- tryCatch({
      MASS::cov.rob(mat, method = "mcd")$cov
    }, error = function(e) NULL)

    # Fallback to diagonal if robust estimation fails
    if (is.null(covmat) || any(is.na(covmat)) || det(covmat) <= 0) {
      v1 <- var(mat[, 1], na.rm = TRUE)
      v2 <- var(mat[, 2], na.rm = TRUE)
      covmat <- matrix(c(v1, 0, 0, v2), 2, 2)
    }
    return(covmat)
  }

  #--------------------------------------#
  # Helper: Leave-one-out predictions
  #--------------------------------------#
  predict_loo <- function(B_loo, X) {
    n_params <- ncol(X)
    Bp <- B_loo[, 1:n_params, drop = FALSE]
    Ba <- B_loo[, (n_params + 1):(2 * n_params), drop = FALSE]

    pred_primary <- rowSums(Bp * X)
    pred_aux <- rowSums(Ba * X)

    list(pred_primary = pred_primary, pred_aux = pred_aux)
  }

  #--------------------------------------#
  # Helper: Compute RSS using leave-one-out
  #--------------------------------------#
  compute_rss_loo <- function(d, return_details = FALSE) {
    # Prepare weighted data
    X <- as.matrix(d[, BetaExposure, drop = FALSE]) * sqrt(d$Weights)
    Y_primary <- d[[BetaOutcome]] * sqrt(d$Weights)
    Y_aux <- d[[BetaAux]] * sqrt(d$Weights)

    n <- nrow(X)
    n_params <- ncol(X)

    # Leave-one-out coefficient estimation
    B_loo <- matrix(NA, nrow = n, ncol = 2 * n_params)

    for (i in seq_len(n)) {
      X_loo <- X[-i, , drop = FALSE]
      Y_primary_loo <- Y_primary[-i]
      Y_aux_loo <- Y_aux[-i]

      # Solve for coefficients
      XtX_inv <- tryCatch(
        solve(t(X_loo) %*% X_loo),
        error = function(e) NULL
      )

      if (!is.null(XtX_inv)) {
        beta_primary <- XtX_inv %*% t(X_loo) %*% Y_primary_loo
        beta_aux <- XtX_inv %*% t(X_loo) %*% Y_aux_loo
        B_loo[i, ] <- c(beta_primary, beta_aux)
      }
    }

    # Compute predictions and residuals
    predictions <- predict_loo(B_loo, X)
    residuals <- cbind(
      Y_primary - predictions$pred_primary,
      Y_aux - predictions$pred_aux
    )

    # Robust covariance of residuals
    cov_residuals <- compute_robust_cov(residuals)

    # Mahalanobis distance (multivariate outlier detection)
    mahal_dist <- mahalanobis(
      residuals,
      center = c(0, 0),
      cov = cov_residuals,
      inverted = FALSE
    )

    RSS <- sum(mahal_dist, na.rm = TRUE)

    if (return_details) {
      list(
        RSS = RSS,
        B_loo = B_loo,
        cov_residuals = cov_residuals,
        residuals = residuals,
        mahal_dist = mahal_dist
      )
    } else {
      RSS
    }
  }

  #--------------------------------------#
  # Helper: Generate random dataset under null
  #--------------------------------------#
  generate_null_data <- function(d) {
    n <- nrow(d)

    # Simulate exposure values
    X_sim <- rnorm(n, d[[BetaExposure]], d[[SdExposure]])

    # Fit leave-one-out models
    models_primary <- lapply(seq_len(n), function(i) {
      lm(
        as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)),
        weights = Weights,
        data = d[-i, ]
      )
    })

    models_aux <- lapply(seq_len(n), function(i) {
      lm(
        as.formula(paste0(BetaAux, " ~ -1 + ", BetaExposure)),
        weights = Weights,
        data = d[-i, ]
      )
    })

    # Simulate outcomes based on LOO predictions
    Y_primary_sim <- mapply(
      function(model, x_val, se) {
        pred <- predict(model, newdata = setNames(data.frame(x_val), BetaExposure))
        rnorm(1, pred, se)
      },
      models_primary, X_sim, d[[SdOutcome]]
    )

    Y_aux_sim <- mapply(
      function(model, x_val, se) {
        pred <- predict(model, newdata = setNames(data.frame(x_val), BetaExposure))
        rnorm(1, pred, se)
      },
      models_aux, X_sim, d[[SdAux]]
    )

    # Create simulated dataset
    sim_data <- d
    sim_data[[BetaExposure]] <- X_sim
    sim_data[[BetaOutcome]] <- Y_primary_sim
    sim_data[[BetaAux]] <- Y_aux_sim
    sim_data$Weights <- 1 / (d[[SdOutcome]] * d[[SdAux]])

    return(sim_data)
  }

  #--------------------------------------#
  # Main computation: Generate null distribution
  #--------------------------------------#
  message(paste0("Generating ", NbDistribution, " null distributions..."))

  if (n_cores > 1) {
    simulated_data <- parallel::mclapply(
      1:NbDistribution,
      function(i) generate_null_data(core),
      mc.cores = n_cores
    )
    RSS_null <- unlist(parallel::mclapply(
      simulated_data,
      compute_rss_loo,
      mc.cores = n_cores
    ))
  } else {
    simulated_data <- replicate(NbDistribution, generate_null_data(core), simplify = FALSE)
    RSS_null <- sapply(simulated_data, compute_rss_loo)
  }

  #--------------------------------------#
  # Observed statistics
  #--------------------------------------#
  obs_results <- compute_rss_loo(core, return_details = OUTLIERtest)
  RSS_observed <- if (OUTLIERtest) obs_results$RSS else obs_results

  # Global test p-value
  p_global <- mean(RSS_null >= RSS_observed)

  # Fit full model (without outlier removal)
  fit_full <- lm(
    as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)),
    weights = Weights1,
    data = core
  )

  beta_raw <- coef(fit_full)[1]
  se_raw <- summary(fit_full)$coefficients[1, "Std. Error"]

  # Initialize corrected estimates
  beta_corrected <- NA
  se_corrected <- NA
  outlier_indices <- NULL
  n_outliers <- 0

  #--------------------------------------#
  # Outlier detection and correction
  #--------------------------------------#
  if (OUTLIERtest && p_global < SignifThreshold) {
    message("Significant pleiotropy detected. Identifying outliers...")

    mahal_dist <- obs_results$mahal_dist
    outlier_threshold <- quantile(mahal_dist, 1 - SignifThreshold, na.rm = TRUE)
    outlier_indices <- which(mahal_dist > outlier_threshold)
    n_outliers <- length(outlier_indices)

    if (n_outliers > 0 && n_outliers < nrow(core)) {
      message(paste0("Found ", n_outliers, " outlier(s). Re-estimating without outliers..."))

      fit_corrected <- lm(
        as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)),
        weights = Weights1,
        data = core[-outlier_indices, ]
      )

      beta_corrected <- coef(fit_corrected)[1]
      se_corrected <- summary(fit_corrected)$coefficients[1, "Std. Error"]
    } else {
      message("No valid outliers to remove or too many outliers detected.")
    }
  }

  #--------------------------------------#
  # Distortion test
  #--------------------------------------#
  p_distortion <- NA

  if (DISTORTIONtest && !is.na(beta_corrected)) {
    message("Running distortion test...")

    bias_observed <- (beta_raw - beta_corrected) / abs(beta_corrected)

    # Generate null distribution of bias
    bias_null <- replicate(NbDistribution, {
      random_indices <- sample(seq_len(nrow(core)), n_outliers)
      fit_random <- update(fit_full, data = core[-random_indices, ])
      beta_random <- coef(fit_random)[1]
      (beta_raw - beta_random) / abs(beta_random)
    })

    p_distortion <- mean(abs(bias_null) >= abs(bias_observed), na.rm = TRUE)
  }

  #--------------------------------------#
  # Return results
  #--------------------------------------#
  message("Analysis complete.")

  results <- list(
    raw_beta = beta_raw,
    raw_se = se_raw,
    corrected_beta = beta_corrected,
    corrected_se = se_corrected,
    p_value = p_global,
    p_distort = p_distortion,
    outlier_idx = outlier_indices,
    n_instruments = n_instruments,
    n_outliers = n_outliers
  )

  class(results) <- c("mrpresso_ib", "list")
  return(results)
}

#' Print method for mrpresso_ib
#' @param x An object of class mrpresso_ib
#' @param ... Additional arguments (not used)
#' @export
print.mrpresso_ib <- function(x, ...) {
  cat("MR-PRESSO with Instrument Borrowing Results\n")
  cat("============================================\n\n")
  cat(sprintf("Number of instruments: %d\n", x$n_instruments))
  cat(sprintf("Global test p-value: %.4f\n\n", x$p_value))

  cat("Raw estimate (before outlier correction):\n")
  cat(sprintf("  Beta: %.4f (SE: %.4f)\n\n", x$raw_beta, x$raw_se))

  if (!is.na(x$corrected_beta)) {
    cat(sprintf("Number of outliers detected: %d\n", x$n_outliers))
    cat("Corrected estimate (after outlier removal):\n")
    cat(sprintf("  Beta: %.4f (SE: %.4f)\n\n", x$corrected_beta, x$corrected_se))

    if (!is.na(x$p_distort)) {
      cat(sprintf("Distortion test p-value: %.4f\n", x$p_distort))
    }
  } else {
    cat("No outliers detected or correction not performed.\n")
  }
}
