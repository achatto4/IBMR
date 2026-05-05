#' Calculate Coheterogeneity Using a Guarded Closed-Form Method
#'
#' Computes pairwise coheterogeneity correlations across multiple outcome traits
#' using the bias-corrected moment estimator and closed-form plug-in standard
#' error described for instrument-borrowing Mendelian randomization.
#'
#' @param BetaXG Numeric vector of SNP-exposure associations.
#' @param BetaYG_matrix Matrix of SNP-outcome associations, rows = SNPs,
#'   columns = traits. Optional if `BetaIV_matrix` is provided.
#' @param seBetaXG Numeric vector of standard errors for `BetaXG`.
#' @param seBetaYG_matrix Matrix of standard errors for `BetaYG_matrix`.
#'   Optional if `seBetaIV_matrix` is provided.
#' @param BetaIV_matrix Matrix of ratio estimates (`BetaYG / BetaXG`), rows =
#'   SNPs, columns = traits. Used to reconstruct `BetaYG_matrix` when needed.
#' @param seBetaIV_matrix Matrix of standard errors for `BetaIV_matrix`. Used
#'   to reconstruct `seBetaYG_matrix` when needed.
#' @param ldsc_intercepts Optional matrix of LDSC intercepts used for pairwise
#'   outcome covariance adjustment.
#' @param SNP_keep Optional logical or integer index vector specifying the SNPs
#'   to retain.
#' @param use_ldsc Logical; if `TRUE`, use `ldsc_intercepts` when supplied.
#' @param alpha Significance level retained in the returned guard settings.
#' @param eps Small positive constant used for numerical stabilization.
#' @param bx_min Minimum absolute exposure effect size allowed.
#' @param F_min Minimum first-stage F statistic threshold, where
#'   `F = (BetaXG / seBetaXG)^2`.
#' @param min_K_pair Minimum number of valid SNPs required for a trait pair.
#' @param return_diagnostics Logical; if `TRUE`, include a diagnostics list.
#'
#' @return A list containing:
#' \item{rho}{Pairwise coheterogeneity correlation matrix.}
#' \item{se}{Pairwise closed-form standard error matrix for `rho`.}
#' \item{z_statistic}{Pairwise z-statistic matrix.}
#' \item{wald_statistic}{Pairwise Wald statistic matrix.}
#' \item{p_value}{Pairwise p-value matrix.}
#' \item{K}{Matrix of valid SNP counts per trait pair.}
#' \item{flag}{Matrix of per-pair diagnostic flags.}
#' \item{method}{Method label for the estimator.}
#' \item{guards}{List of guardrail settings used in the analysis.}
#' \item{diagnostics}{Optional detailed diagnostics by trait pair.}
#'
#' @export
coheterogeneity_Q <- function(
    BetaXG,
    BetaYG_matrix = NULL,
    seBetaXG,
    seBetaYG_matrix = NULL,
    BetaIV_matrix = NULL,
    seBetaIV_matrix = NULL,
    ldsc_intercepts = NULL,
    SNP_keep = NULL,
    use_ldsc = TRUE,
    alpha = 0.05,
    eps = 1e-12,
    bx_min = 1e-6,
    F_min = 0,
    min_K_pair = 50,
    return_diagnostics = FALSE
) {
  stopifnot(!is.null(BetaXG), !is.null(seBetaXG))
  K0 <- length(BetaXG)

  if (length(seBetaXG) != K0) {
    stop("seBetaXG must have the same length as BetaXG.")
  }

  if (is.null(SNP_keep)) {
    SNP_keep <- rep(TRUE, K0)
  } else if (is.logical(SNP_keep)) {
    stopifnot(length(SNP_keep) == K0)
  } else {
    keep <- rep(FALSE, K0)
    keep[SNP_keep] <- TRUE
    SNP_keep <- keep
  }

  if (is.null(BetaYG_matrix)) {
    if (is.null(BetaIV_matrix)) {
      stop("Provide BetaYG_matrix or BetaIV_matrix.")
    }
    BetaYG_matrix <- sweep(as.matrix(BetaIV_matrix), 1, BetaXG, FUN = "*")
  } else {
    BetaYG_matrix <- as.matrix(BetaYG_matrix)
  }

  if (is.null(seBetaYG_matrix)) {
    if (is.null(seBetaIV_matrix)) {
      stop("Provide seBetaYG_matrix or seBetaIV_matrix.")
    }
    seBetaYG_matrix <- sweep(as.matrix(seBetaIV_matrix), 1, abs(BetaXG), FUN = "*")
  } else {
    seBetaYG_matrix <- as.matrix(seBetaYG_matrix)
  }

  if (nrow(BetaYG_matrix) != K0 || nrow(seBetaYG_matrix) != K0) {
    stop("Outcome matrices must have one row per SNP in BetaXG.")
  }
  if (!all(dim(BetaYG_matrix) == dim(seBetaYG_matrix))) {
    stop("BetaYG_matrix and seBetaYG_matrix must have the same dimensions.")
  }

  BetaXG <- BetaXG[SNP_keep]
  seBetaXG <- seBetaXG[SNP_keep]
  BetaYG_matrix <- BetaYG_matrix[SNP_keep, , drop = FALSE]
  seBetaYG_matrix <- seBetaYG_matrix[SNP_keep, , drop = FALSE]

  J <- ncol(BetaYG_matrix)
  if (J < 2) {
    stop("At least two outcome traits are required.")
  }

  Fstat <- (BetaXG / seBetaXG)^2
  ok_x <- is.finite(BetaXG) & is.finite(seBetaXG) & seBetaXG > 0 &
    abs(BetaXG) > bx_min & is.finite(Fstat) & (Fstat > F_min)

  rho_matrix <- matrix(NA_real_, J, J)
  se_matrix <- matrix(NA_real_, J, J)
  z_matrix <- matrix(NA_real_, J, J)
  wald_matrix <- matrix(NA_real_, J, J)
  p_matrix <- matrix(NA_real_, J, J)
  K_matrix <- matrix(0L, J, J)
  flag_matrix <- matrix("", J, J)

  if (!is.null(colnames(BetaYG_matrix))) {
    dimnames(rho_matrix) <- list(colnames(BetaYG_matrix), colnames(BetaYG_matrix))
    dimnames(se_matrix) <- dimnames(rho_matrix)
    dimnames(z_matrix) <- dimnames(rho_matrix)
    dimnames(wald_matrix) <- dimnames(rho_matrix)
    dimnames(p_matrix) <- dimnames(rho_matrix)
    dimnames(K_matrix) <- dimnames(rho_matrix)
    dimnames(flag_matrix) <- dimnames(rho_matrix)
  }

  diagnostics_list <- if (return_diagnostics) {
    list()
  } else {
    NULL
  }

  rho_from_betas_guarded <- function(bx, by1, by2, sebx, sey1, sey2, I12) {
    ok <- ok_x &
      complete.cases(by1, by2, sey1, sey2) &
      is.finite(by1) & is.finite(by2) &
      is.finite(sey1) & is.finite(sey2) &
      sey1 > 0 & sey2 > 0

    if (sum(ok) < min_K_pair) {
      return(list(rho = NA_real_, K = sum(ok), flag = "K<min", cache = NULL))
    }

    bx <- bx[ok]
    by1 <- by1[ok]
    by2 <- by2[ok]
    sebx <- sebx[ok]
    sey1 <- sey1[ok]
    sey2 <- sey2[ok]

    vbx <- sebx^2
    bx2 <- bx^2 + eps
    bx4 <- bx2^2

    theta1 <- by1 / bx
    theta2 <- by2 / bx

    sigma1 <- (sey1^2) / bx2 + (by1^2) * vbx / bx4
    sigma2 <- (sey2^2) / bx2 + (by2^2) * vbx / bx4
    cov_by12 <- I12 * sey1 * sey2
    sigma12 <- cov_by12 / bx2 + (by1 * by2) * vbx / bx4

    ok_sigma <- is.finite(theta1) & is.finite(theta2) &
      is.finite(sigma1) & is.finite(sigma2) & is.finite(sigma12) &
      sigma1 > 0 & sigma2 > 0

    if (sum(ok_sigma) < min_K_pair) {
      return(list(rho = NA_real_, K = sum(ok_sigma), flag = "bad sigma", cache = NULL))
    }

    bx <- bx[ok_sigma]
    by1 <- by1[ok_sigma]
    by2 <- by2[ok_sigma]
    sebx <- sebx[ok_sigma]
    sey1 <- sey1[ok_sigma]
    sey2 <- sey2[ok_sigma]
    theta1 <- theta1[ok_sigma]
    theta2 <- theta2[ok_sigma]
    sigma1 <- sigma1[ok_sigma]
    sigma2 <- sigma2[ok_sigma]
    sigma12 <- sigma12[ok_sigma]
    cov_by12 <- cov_by12[ok_sigma]

    Kp <- length(theta1)
    w <- 1 / sqrt((sigma1 + eps) * (sigma2 + eps))
    if (any(!is.finite(w)) || sum(w) <= 0) {
      return(list(rho = NA_real_, K = Kp, flag = "bad w", cache = NULL))
    }
    w <- w / sum(w)

    delta1 <- theta1 - sum(w * theta1)
    delta2 <- theta2 - sum(w * theta2)

    C12 <- sum(w * (delta1 * delta2 - sigma12))
    tau1_sq <- max(sum(w * (delta1^2 - sigma1)), 0)
    tau2_sq <- max(sum(w * (delta2^2 - sigma2)), 0)

    if (!is.finite(C12) || !is.finite(tau1_sq) || !is.finite(tau2_sq)) {
      return(list(rho = NA_real_, K = Kp, flag = "moment nonfinite", cache = NULL))
    }

    if (tau1_sq <= eps || tau2_sq <= eps) {
      return(list(rho = NA_real_, K = Kp, flag = "tau0", cache = NULL))
    }

    tau1 <- sqrt(tau1_sq)
    tau2 <- sqrt(tau2_sq)
    bound <- tau1 * tau2
    C12_clamped <- sign(C12) * min(abs(C12), bound)
    rho <- C12_clamped / bound

    if (!is.finite(rho)) {
      return(list(rho = NA_real_, K = Kp, flag = "rho nonfinite", cache = NULL))
    }
    rho <- max(-1, min(1, rho))

    cache <- list(
      bx = bx,
      by1 = by1,
      by2 = by2,
      sebx = sebx,
      sey1 = sey1,
      sey2 = sey2,
      cov_by12 = cov_by12,
      theta1 = theta1,
      theta2 = theta2,
      delta1 = delta1,
      delta2 = delta2,
      w = w,
      tau1 = tau1,
      tau2 = tau2,
      C12 = C12,
      C12_clamped = C12_clamped,
      I12 = I12
    )

    list(rho = rho, K = Kp, flag = "OK", cache = cache)
  }

  se_closed_form <- function(cache, rho) {
    bx <- cache$bx
    by1 <- cache$by1
    by2 <- cache$by2
    sebx <- cache$sebx
    sey1 <- cache$sey1
    sey2 <- cache$sey2
    cov_by12 <- cache$cov_by12
    theta1 <- cache$theta1
    theta2 <- cache$theta2
    delta1 <- cache$delta1
    delta2 <- cache$delta2
    w <- cache$w
    tau1 <- cache$tau1
    tau2 <- cache$tau2

    D1 <- delta1 - rho * (tau1 / tau2) * delta2
    D2 <- delta2 - rho * (tau2 / tau1) * delta1

    var_terms <- (theta1 * D2 + theta2 * D1)^2 * sebx^2 +
      D2^2 * sey1^2 +
      D1^2 * sey2^2 +
      2 * D1 * D2 * cov_by12

    var_contrib <- ((w / bx)^2) * var_terms
    var_rho <- sum(var_contrib) / (tau1^2 * tau2^2)

    if (!is.finite(var_rho) || var_rho <= 0) {
      return(NA_real_)
    }

    sqrt(var_rho)
  }

  for (j in seq_len(J - 1)) {
    for (l in (j + 1):J) {
      I12 <- 0
      if (use_ldsc && !is.null(ldsc_intercepts)) {
        I12 <- ldsc_intercepts[j, l]
        if (!is.finite(I12)) {
          I12 <- 0
        }
      }

      out <- rho_from_betas_guarded(
        bx = BetaXG,
        by1 = BetaYG_matrix[, j],
        by2 = BetaYG_matrix[, l],
        sebx = seBetaXG,
        sey1 = seBetaYG_matrix[, j],
        sey2 = seBetaYG_matrix[, l],
        I12 = I12
      )

      rho_matrix[j, l] <- rho_matrix[l, j] <- out$rho
      K_matrix[j, l] <- K_matrix[l, j] <- out$K
      flag_matrix[j, l] <- flag_matrix[l, j] <- out$flag

      if (out$flag != "OK") {
        next
      }

      s <- se_closed_form(out$cache, out$rho)
      se_matrix[j, l] <- se_matrix[l, j] <- s

      if (!is.finite(s) || s <= 0) {
        flag_matrix[j, l] <- flag_matrix[l, j] <- "se NA"
        next
      }

      z <- out$rho / s
      wald <- z^2
      pval <- 2 * pnorm(-abs(z))

      z_matrix[j, l] <- z_matrix[l, j] <- z
      wald_matrix[j, l] <- wald_matrix[l, j] <- wald
      p_matrix[j, l] <- p_matrix[l, j] <- pval

      if (return_diagnostics) {
        diagnostics_list[[paste0(j, "_", l)]] <- list(
          K = out$K,
          I12 = I12,
          rho = out$rho,
          se = s,
          z = z,
          wald = wald,
          p = pval,
          flag = out$flag,
          C12 = out$cache$C12,
          C12_clamped = out$cache$C12_clamped,
          tau1 = out$cache$tau1,
          tau2 = out$cache$tau2
        )
      }
    }
  }

  diag(rho_matrix) <- 1
  diag(se_matrix) <- 0
  diag(z_matrix) <- Inf
  diag(wald_matrix) <- Inf
  diag(p_matrix) <- 0
  diag(flag_matrix) <- "diag"

  res <- list(
    rho = rho_matrix,
    se = se_matrix,
    z_statistic = z_matrix,
    wald_statistic = wald_matrix,
    p_value = p_matrix,
    K = K_matrix,
    flag = flag_matrix,
    method = "closed_form_coheterogeneity_guarded",
    guards = list(
      bx_min = bx_min,
      F_min = F_min,
      min_K_pair = min_K_pair,
      alpha = alpha
    )
  )

  if (return_diagnostics) {
    res$diagnostics <- diagnostics_list
  }

  class(res) <- c("coheterogeneity", "list")
  res
}
