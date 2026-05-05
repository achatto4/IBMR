#' Calculate Coheterogeneity Using a Guarded Theoretical Method
#'
#' Computes pairwise coheterogeneity correlations across multiple outcome traits
#' using a guarded delta-exact approach. The implementation supports optional
#' SNP filtering, weak-instrument filtering, winsorization of ratio estimates,
#' LDSC intercept adjustment, and per-pair diagnostics.
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
#'   covariance adjustment.
#' @param SNP_keep Optional logical or integer index vector specifying the SNPs
#'   to retain.
#' @param use_ldsc Logical; if `TRUE`, use `ldsc_intercepts` when supplied.
#' @param alpha Significance level retained for interface compatibility.
#' @param step_mult Scalar used to scale finite-difference steps for the delta
#'   method standard error.
#' @param eps Small positive constant used for numerical stabilization.
#' @param bx_min Minimum absolute exposure effect size allowed.
#' @param F_min Minimum first-stage F statistic threshold, where
#'   `F = (BetaXG / seBetaXG)^2`.
#' @param winsor_theta Optional cap applied to ratio estimates to reduce
#'   instability from very small exposure effects. Set to `NULL` to disable.
#' @param min_K_pair Minimum number of valid SNPs required for a trait pair.
#' @param return_diagnostics Logical; if `TRUE`, include a diagnostics list.
#'
#' @return A list containing:
#' \item{rho}{Pairwise coheterogeneity correlation matrix.}
#' \item{se}{Pairwise standard error matrix for `rho`.}
#' \item{z_statistic}{Pairwise z-statistic matrix.}
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
    Kp <- length(bx)

    vbx <- sebx^2
    bx2 <- bx^2 + eps
    bx4 <- bx2^2

    th1 <- by1 / bx
    th2 <- by2 / bx

    if (!is.null(winsor_theta)) {
      th1 <- pmax(pmin(th1, winsor_theta), -winsor_theta)
      th2 <- pmax(pmin(th2, winsor_theta), -winsor_theta)
      by1 <- th1 * bx
      by2 <- th2 * bx
    }

    s1 <- (sey1^2) / bx2 + (by1^2) * vbx / bx4
    s2 <- (sey2^2) / bx2 + (by2^2) * vbx / bx4

    cov_by12 <- I12 * sey1 * sey2
    s12 <- cov_by12 / bx2 + (by1 * by2) * vbx / bx4

    ok2 <- is.finite(s1) & is.finite(s2) & is.finite(s12) & s1 > 0 & s2 > 0

    if (sum(ok2) < min_K_pair) {
      return(list(rho = NA_real_, K = sum(ok2), flag = "bad sigma", cache = NULL))
    }

    bx_final <- bx[ok2]
    by1_final <- by1[ok2]
    by2_final <- by2[ok2]
    sebx_final <- sebx[ok2]
    sey1_final <- sey1[ok2]
    sey2_final <- sey2[ok2]
    th1 <- th1[ok2]
    th2 <- th2[ok2]
    s1 <- s1[ok2]
    s2 <- s2[ok2]
    s12 <- s12[ok2]

    Kp <- length(th1)
    w <- 1 / sqrt((s1 + eps) * (s2 + eps))
    if (any(!is.finite(w)) || sum(w) <= 0) {
      return(list(rho = NA_real_, K = Kp, flag = "bad w", cache = NULL))
    }
    w <- w / sum(w)

    d1 <- th1 - sum(w * th1)
    d2 <- th2 - sum(w * th2)

    M11 <- sum(w * d1^2)
    M22 <- sum(w * d2^2)
    M12 <- sum(w * d1 * d2)

    B11 <- sum(w * (1 - w) * s1)
    B22 <- sum(w * (1 - w) * s2)
    B12 <- sum(w * (1 - w) * s12)

    S11 <- M11 - B11
    S22 <- M22 - B22
    S12 <- M12 - B12

    if (!is.finite(S11) || !is.finite(S22)) {
      return(list(rho = NA_real_, K = Kp, flag = "Sii nonfinite", cache = NULL))
    }

    if (S11 <= eps || S22 <= eps) {
      return(list(rho = NA_real_, K = Kp, flag = "tau0", cache = NULL))
    }

    bound <- sqrt(S11 * S22)
    S12c <- sign(S12) * min(abs(S12), bound)

    rho <- S12c / bound
    if (!is.finite(rho)) {
      return(list(rho = NA_real_, K = Kp, flag = "rho nonfinite", cache = NULL))
    }
    rho <- max(-1, min(1, rho))

    list(
      rho = rho,
      K = Kp,
      flag = "OK",
      cache = list(
        bx = bx_final,
        by1 = by1_final,
        by2 = by2_final,
        sebx = sebx_final,
        sey1 = sey1_final,
        sey2 = sey2_final,
        I12 = I12
      )
    )
  }

  numeric_grad <- function(f, x, step) {
    f0 <- f(x)
    if (!is.finite(f0)) {
      return(rep(NA_real_, length(x)))
    }

    g <- numeric(length(x))
    for (i in seq_along(x)) {
      h <- step[i]
      xp <- x
      xm <- x
      xp[i] <- xp[i] + h
      xm[i] <- xm[i] - h
      fp <- f(xp)
      fm <- f(xm)

      if (!is.finite(fp) || !is.finite(fm)) {
        return(rep(NA_real_, length(x)))
      }
      g[i] <- (fp - fm) / (2 * h)
    }
    g
  }

  se_delta_exact <- function(bx, by1, by2, sebx, sey1, sey2, I12) {
    Kp <- length(bx)
    if (Kp < 10) {
      warning(sprintf("Only %d SNPs for SE estimation; unstable", Kp))
      return(NA_real_)
    }

    x <- c(bx, by1, by2)

    f <- function(xx) {
      bbx <- xx[1:Kp]
      bby1 <- xx[(Kp + 1):(2 * Kp)]
      bby2 <- xx[(2 * Kp + 1):(3 * Kp)]

      vvbx <- sebx^2
      bbx2 <- bbx^2 + eps
      bbx4 <- bbx2^2

      tth1 <- bby1 / bbx
      tth2 <- bby2 / bbx

      ss1 <- (sey1^2) / bbx2 + (bby1^2) * vvbx / bbx4
      ss2 <- (sey2^2) / bbx2 + (bby2^2) * vvbx / bbx4
      ss12 <- (I12 * sey1 * sey2) / bbx2 + (bby1 * bby2) * vvbx / bbx4

      if (any(!is.finite(ss1)) || any(!is.finite(ss2)) ||
          any(ss1 <= 0) || any(ss2 <= 0)) {
        return(NA_real_)
      }

      ww <- 1 / sqrt((ss1 + eps) * (ss2 + eps))
      if (any(!is.finite(ww)) || sum(ww) <= 0) {
        return(NA_real_)
      }
      ww <- ww / sum(ww)

      dd1 <- tth1 - sum(ww * tth1)
      dd2 <- tth2 - sum(ww * tth2)

      MM11 <- sum(ww * dd1^2)
      MM22 <- sum(ww * dd2^2)
      MM12 <- sum(ww * dd1 * dd2)

      BB11 <- sum(ww * (1 - ww) * ss1)
      BB22 <- sum(ww * (1 - ww) * ss2)
      BB12 <- sum(ww * (1 - ww) * ss12)

      SS11 <- MM11 - BB11
      SS22 <- MM22 - BB22
      SS12 <- MM12 - BB12

      if (!is.finite(SS11) || !is.finite(SS22) || SS11 <= eps || SS22 <= eps) {
        return(NA_real_)
      }

      bound <- sqrt(SS11 * SS22)
      SS12c <- sign(SS12) * min(abs(SS12), bound)
      rho_val <- SS12c / bound

      if (!is.finite(rho_val)) {
        return(NA_real_)
      }
      max(-1, min(1, rho_val))
    }

    step <- c(sebx, sey1, sey2) * step_mult
    step[!is.finite(step) | step == 0] <- 1e-8

    g <- numeric_grad(f, x, step)
    if (any(!is.finite(g))) {
      warning("Non-finite gradient in SE computation")
      return(NA_real_)
    }

    cov_by12 <- I12 * sey1 * sey2
    var <- 0

    for (k in seq_len(Kp)) {
      gk <- c(g[k], g[Kp + k], g[2 * Kp + k])
      Ok <- matrix(c(
        sebx[k]^2, 0, 0,
        0, sey1[k]^2, cov_by12[k],
        0, cov_by12[k], sey2[k]^2
      ), 3, 3, byrow = TRUE)

      var_k <- as.numeric(t(gk) %*% Ok %*% gk)
      if (!is.finite(var_k)) {
        warning(sprintf("Non-finite variance contribution at SNP %d", k))
        return(NA_real_)
      }
      var <- var + var_k
    }

    se_val <- sqrt(max(var, 0))
    if (!is.finite(se_val) || se_val <= 0) {
      warning("Invalid final SE")
      return(NA_real_)
    }

    se_val
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

      cache <- out$cache
      if (is.null(cache) || is.null(cache$bx) || length(cache$bx) < min_K_pair) {
        flag_matrix[j, l] <- flag_matrix[l, j] <- "cache invalid"
        next
      }

      s <- tryCatch(
        se_delta_exact(
          bx = cache$bx,
          by1 = cache$by1,
          by2 = cache$by2,
          sebx = cache$sebx,
          sey1 = cache$sey1,
          sey2 = cache$sey2,
          I12 = cache$I12
        ),
        error = function(e) {
          warning(sprintf("SE error for pair (%d,%d): %s", j, l, e$message))
          NA_real_
        }
      )

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
          p = pval,
          flag = out$flag
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
    p_value = p_matrix,
    K = K_matrix,
    flag = flag_matrix,
    method = "theoretical_delta_exact_guarded",
    guards = list(
      bx_min = bx_min,
      F_min = F_min,
      winsor_theta = winsor_theta,
      min_K_pair = min_K_pair,
      step_mult = step_mult,
      alpha = alpha
    )
  )

  if (return_diagnostics) {
    res$diagnostics <- diagnostics_list
  }

  class(res) <- c("coheterogeneity", "list")
  res
}
