#' IB-PRESSO: instrument borrowing for outlier-robust Mendelian randomization
#'
#' Extends MR-PRESSO (Verbanck et al., 2018) by scoring each instrument on a
#' primary outcome AND an auxiliary outcome simultaneously. Every other element
#' of the procedure -- leave-one-out IVW fitting through the origin, the global
#' heterogeneity bootstrap, the Bonferroni correction across instruments, and
#' the outlier-corrected IVW estimate -- is MR-PRESSO's, unchanged, so that any
#' difference in performance is attributable to the auxiliary trait and to
#' nothing else.
#'
#' @section Per-SNP statistic:
#' Let \eqn{r_{1j}} and \eqn{r_{2j}} be the leave-one-out residuals of the
#' primary and auxiliary outcomes at instrument \eqn{j}, on the raw scale of the
#' reported effect estimates. Conditional on the leave-one-out slopes
#' \eqn{b_{1,-j}} and \eqn{b_{2,-j}} their sampling covariance is known:
#' \deqn{V_j = \begin{pmatrix} s_{Y_1j}^2 + b_{1,-j}^2 s_{Xj}^2 &
#'   \rho_j s_{Y_1j} s_{Y_2j} + b_{1,-j} b_{2,-j} s_{Xj}^2 \\
#'   \rho_j s_{Y_1j} s_{Y_2j} + b_{1,-j} b_{2,-j} s_{Xj}^2 &
#'   s_{Y_2j}^2 + b_{2,-j}^2 s_{Xj}^2 \end{pmatrix},}
#' and the statistic is the generalized squared residual
#' \deqn{T_j = (r_{1j}, r_{2j}) V_j^{-1} (r_{1j}, r_{2j})^\top
#'  = \frac{v_{2j} r_{1j}^2 - 2 c_{12,j} r_{1j} r_{2j} + v_{1j} r_{2j}^2}
#'         {v_{1j} v_{2j} - c_{12,j}^2}.}
#' MR-PRESSO's statistic is the primary term alone. Two dependence sources enter
#' \eqn{V_j}, each exactly once: outcome-GWAS sample overlap through
#' \eqn{\rho_j}, and shared exposure noise through \eqn{b_{1,-j} b_{2,-j}
#' s_{Xj}^2}, which is present even when \eqn{\rho_j = 0}.
#'
#' @section Why the covariance is not estimated:
#' \eqn{V_j} is built from the reported GWAS standard errors and a cross-trait
#' LDSC intercept, never from the observed residuals. Estimating a
#' \eqn{2 \times 2} residual covariance (robustly or otherwise) fails in three
#' ways: the invalid instruments inflate the estimate along precisely the
#' direction that identifies them; the resulting Mahalanobis distance is bounded
#' by \eqn{1/\pi} where \eqn{\pi} is the fraction of invalid instruments, so at
#' \eqn{\pi = 0.5} no instrument can exceed 2 against a cutoff of 5.99; and the
#' minimum covariance determinant has a breakdown point of 50 percent.
#'
#' @section Reference distributions:
#' \eqn{T_j} is pivotal -- a quadratic form in a bivariate normal with known
#' covariance -- so conditional on the leave-one-out slopes its null is exactly
#' \eqn{\chi^2_2} and the per-SNP p-value is analytic by default. This removes
#' the bootstrap resolution floor of \eqn{K/\code{NbDistribution}}, which at
#' \eqn{K = 1270} and 10000 draws is 0.127 and makes flagging impossible at any
#' sane threshold. The GLOBAL test keeps the bootstrap, because the
#' leave-one-out slopes are refitted on every draw and \eqn{\sum_j T_j} has no
#' closed form. Setting `OutlierTestMethod = "bootstrap"` restores the previous
#' per-SNP behaviour; both sets of p-values are always returned.
#'
#' Neither route propagates the uncertainty in the leave-one-out slopes, so on
#' self-fitted null data \eqn{T_j} is inflated by roughly \eqn{1 + 2/K}
#' (measured mean 2.03 at \eqn{K = 100}, 2.02 at \eqn{K = 200}, against 2).
#'
#' @section Weights:
#' Two weights, with a deliberate division of labour.
#' \eqn{W_j = 1/(s_{Y_1j} s_{Y_2j})} fits the leave-one-out slopes, a single
#' weight applied to both outcomes so that neither trait is privileged; it does
#' not enter the statistic, since \eqn{V_j} carries all the scaling. The causal
#' estimate is an IVW fit weighted by \eqn{1/s_{Y_1j}^2}, which is MR-PRESSO's
#' own weight; the final estimator is therefore identical to MR-PRESSO's and the
#' two methods differ only in which instruments they retain.
#'
#' @param BetaOutcome Character, column name for primary outcome SNP effects.
#' @param BetaExposure Character, column name for exposure SNP effects (one only).
#' @param BetaAux Character, column name for auxiliary trait SNP effects.
#' @param SdOutcome Character, column name for standard errors of outcome effects.
#' @param SdExposure Character, column name for standard errors of exposure effects.
#' @param SdAux Character, column name for standard errors of auxiliary effects.
#' @param data Data frame containing all the above columns, one row per instrument.
#' @param OUTLIERtest Logical; perform the per-SNP outlier test and report the
#'   corrected estimate. Default `TRUE`.
#' @param DISTORTIONtest Logical; test whether the removed instruments distort
#'   the estimate more than a random set of the same size. Default `FALSE`.
#' @param SignifThreshold Numeric significance threshold, Bonferroni-corrected
#'   across instruments for the per-SNP test. Default 0.05.
#' @param NbDistribution Integer, number of bootstrap null datasets for the
#'   GLOBAL test. Default 10000.
#' @param seed Integer random seed, or `NULL`.
#' @param n_cores Retained for backward compatibility and ignored.
#' @param CorOutcomeAux Cross-outcome sampling correlation \eqn{\rho_j}: the
#'   cross-trait LDSC intercept on the CORRELATION scale,
#'   \eqn{I_{12}/\sqrt{I_{11} I_{22}}}. Scalar, or a vector of length
#'   `nrow(data)` or of the number of complete instruments. Must lie strictly
#'   inside \eqn{(-1, 1)}; a raw LDSC intercept generally does not, and is not a
#'   correlation. Default 0, which reproduces independent outcome GWAS.
#' @param OutlierTestMethod `"chisq"` (default) for the analytic per-SNP null,
#'   or `"bootstrap"` for the previous Monte Carlo route.
#'
#' @return An object of class `mrpresso_ib`, a list with components:
#' \item{raw_beta, raw_se}{IVW estimate and standard error using all instruments.}
#' \item{corrected_beta, corrected_se}{Outlier-corrected estimate, or `NA` when
#'   no outlier was removed. Callers should fall back to the raw estimate.}
#' \item{p_value}{Global heterogeneity test p-value (bootstrap).}
#' \item{p_distort}{Distortion test p-value, or `NA`.}
#' \item{outlier_idx}{Row indices of the removed instruments.}
#' \item{outlyingness}{Observed per-SNP statistic \eqn{T_j}.}
#' \item{outlier_pvals}{Bonferroni-corrected per-SNP p-values that were used for
#'   flagging; a copy of whichever of the next two `outlier_test_method` names.}
#' \item{outlier_test_method}{`"chisq"` or `"bootstrap"`.}
#' \item{outlier_pvals_chisq}{Analytic per-SNP p-values. Always computed.}
#' \item{outlier_pvals_boot}{Bootstrap per-SNP p-values, or all `NA` when the
#'   analytic route was used.}
#' \item{outcome_correlation}{The \eqn{\rho_j} actually used, after alignment.}
#' \item{covariance_model}{Free-text description of the covariance assumption.}
#' \item{n_instruments, n_outliers}{Counts.}
#'
#' @references
#' Verbanck M, Chen CY, Neale B, Do R (2018). Detection of widespread horizontal
#' pleiotropy in causal relationships inferred from Mendelian randomization
#' between complex traits and diseases. \emph{Nature Genetics} 50, 693-698.
#'
#' @importFrom stats lm coef complete.cases rnorm as.formula pchisq
#' @export
#'
#' @examples
#' set.seed(1)
#' K <- 40
#' bx <- runif(K, 0.02, 0.05)
#' g  <- c(rnorm(12, 0.012, 0.004), rep(0, K - 12))
#' d  <- data.frame(
#'   BetaExposure = rnorm(K, bx, 1/sqrt(1e5)),  SdExposure = 1/sqrt(1e5),
#'   BetaOutcome  = rnorm(K, 0.1 * bx + g, 1/sqrt(5e4)), SdOutcome = 1/sqrt(5e4),
#'   BetaAux      = rnorm(K, 0.3 * bx + g, 1/sqrt(1e5)), SdAux     = 1/sqrt(1e5)
#' )
#' IBPRESSO("BetaOutcome", "BetaExposure", "BetaAux",
#'          "SdOutcome", "SdExposure", "SdAux",
#'          data = d, NbDistribution = 2000, seed = 1, CorOutcomeAux = 0.3)
IBPRESSO <- function(
    BetaOutcome,
    BetaExposure,
    BetaAux,
    SdOutcome,
    SdExposure,
    SdAux,
    data,
    OUTLIERtest     = TRUE,
    DISTORTIONtest  = FALSE,
    SignifThreshold = 0.05,
    NbDistribution  = 10000,
    seed            = NULL,
    n_cores         = 1,       # retained for compatibility; ignored
    ## Cross-outcome sampling correlation: the cross-trait LDSC intercept on
    ## the correlation scale, I12 / sqrt(I11 * I22).  Scalar, or a vector of
    ## length nrow(data).  The default 0 assumes independent outcome GWAS.
    CorOutcomeAux   = 0,
    ## Reference distribution for the PER-SNP test only; the global test is
    ## always the bootstrap.  "chisq" uses the exact pivotal null and is the
    ## default.  "bootstrap" reproduces the previous behaviour and is retained
    ## as a secondary option for comparison.
    OutlierTestMethod = c("chisq", "bootstrap")
) {

  if (!is.null(seed)) set.seed(seed)

  OutlierTestMethod <- match.arg(OutlierTestMethod)
  use_boot_snp      <- identical(OutlierTestMethod, "bootstrap")

  if (length(BetaExposure) != 1L || length(SdExposure) != 1L)
    stop("mr_presso_ib handles a single exposure; BetaExposure and SdExposure ",
         "must each name exactly one column.")

  need <- c(BetaOutcome, BetaExposure, BetaAux, SdOutcome, SdExposure, SdAux)
  miss <- setdiff(need, colnames(data))
  if (length(miss))
    stop("missing columns: ", paste(miss, collapse = ", "))

  ## Keep the complete-case INDEX, not just the filtered frame, so a per-SNP
  ## CorOutcomeAux supplied on the original rows can be aligned to it.
  d0 <- data[, need, drop = FALSE]
  cc <- stats::complete.cases(d0)
  d  <- d0[cc, , drop = FALSE]

  K <- nrow(d)
  if (K <= 3)
    stop("not enough instruments: need at least 4 complete rows, got ", K, ".")

  #--------------------------------------------------------------------------#
  # Cross-outcome sampling correlation
  #--------------------------------------------------------------------------#
  rho <- CorOutcomeAux
  ## A vector aligned to the ORIGINAL rows is subset first, so that a caller
  ## who supplies per-SNP values does not have to know about the complete-case
  ## filter.  A scalar (the normal RDA case) is recycled afterwards.
  if (length(rho) == nrow(data)) rho <- rho[cc]
  if (length(rho) == 1L)         rho <- rep(rho, K)

  if (!is.numeric(rho) || length(rho) != K || any(!is.finite(rho)))
    stop("CorOutcomeAux must be finite and have length 1, nrow(data), ",
         "or the number of complete instruments.")

  ## Not clamped: an out-of-range value indicates an unnormalised intercept or
  ## an unusable LDSC fit, either of which the caller should resolve.
  if (any(abs(rho) >= 1))
    stop("CorOutcomeAux must be strictly between -1 and 1; the supplied value ",
         "does not define a positive-definite outcome sampling covariance. ",
         "A RAW cross-trait LDSC intercept is not a correlation unless the ",
         "univariate intercepts equal 1 -- use I12 / sqrt(I11 * I22), i.e. ",
         "ldsc_cor12() rather than ldsc_I12().")

  ## Same warning MR-PRESSO issues, and for the same reason: the bootstrap tail
  ## is a multiple of 1/NbDistribution, so once K / NbDistribution exceeds
  ## SignifThreshold the Bonferroni-corrected p-value can clear the threshold
  ## only by being exactly zero.  Bootstrap route only; the analytic
  ## chi-square p-value has no resolution floor, which is why it is the default.
  if (OUTLIERtest && use_boot_snp && K / NbDistribution > SignifThreshold)
    warning("NbDistribution is small relative to the number of instruments: ",
            "K / NbDistribution = ", signif(K / NbDistribution, 3), " exceeds ",
            "SignifThreshold = ", SignifThreshold, ", so the per-SNP outlier ",
            "test is below its nominal level. Use NbDistribution >= ",
            ceiling(K / SignifThreshold), ".", call. = FALSE)

  #--------------------------------------------------------------------------#
  # Raw-scale summary statistics.  The statistic is built on the raw scale
  # because V_j is expressed in the units of the reported effect estimates.
  #--------------------------------------------------------------------------#
  bx <- d[[BetaExposure]]
  by <- d[[BetaOutcome]]
  ba <- d[[BetaAux]]

  sx <- d[[SdExposure]]
  sy <- d[[SdOutcome]]
  sa <- d[[SdAux]]

  if (any(!is.finite(c(bx, by, ba, sx, sy, sa))) ||
      any(sx <= 0) || any(sy <= 0) || any(sa <= 0))
    stop("Effect estimates must be finite and all standard errors must be ",
         "positive.")

  ## Symmetric weight, used ONLY to fit the leave-one-out slopes.
  W <- 1 / (sy * sa)

  ## MR-PRESSO's primary-outcome weight, used ONLY for the final causal fit.
  d$.ibmr_w1 <- 1 / sy^2

  #--------------------------------------------------------------------------#
  # Leave-one-out IVW slopes through the origin.
  # Algebraically identical to regressing sqrt(W)-scaled variables, in closed
  # form rather than by refitting K models.
  #--------------------------------------------------------------------------#
  loo <- function(x, y) {
    total_xx <- sum(W * x^2)
    total_xy <- sum(W * x * y)

    den <- total_xx - W * x^2
    num <- total_xy - W * x * y

    tol <- 100 * .Machine$double.eps * max(1, total_xx)

    if (any(!is.finite(den)) || any(den <= tol))
      stop("Leave-one-out slope is undefined because there is insufficient ",
           "weighted exposure variation.")

    num / den
  }

  sl_y <- loo(bx, by)
  sl_a <- loo(bx, ba)

  #--------------------------------------------------------------------------#
  # Analytic residual covariance V_j and the 2x2 quadratic form.
  #
  # V_parts() builds the three distinct entries plus the determinant; quad()
  # evaluates the form.  They are separated so the PER-SNP branch, whose V is
  # constant across bootstrap draws, can build it once instead of 10,000 times.
  # stat_cov() is the combined entry point, used where the slopes vary.
  #--------------------------------------------------------------------------#
  V_parts <- function(slope_y, slope_a) {
    v_y  <- sy^2 + slope_y^2 * sx^2
    v_a  <- sa^2 + slope_a^2 * sx^2
    c_ya <- rho * sy * sa + slope_y * slope_a * sx^2

    det_v <- v_y * v_a - c_ya^2

    tol <- 100 * .Machine$double.eps *
      pmax(v_y * v_a, .Machine$double.xmin)

    if (any(!is.finite(v_y)) || any(!is.finite(v_a)) ||
        any(!is.finite(c_ya)) || any(!is.finite(det_v)) ||
        any(det_v <= tol))
      stop("The residual sampling covariance is not positive definite. ",
           "Check CorOutcomeAux and the supplied standard errors.")

    list(v_y = v_y, v_a = v_a, c_ya = c_ya, det_v = det_v)
  }

  quad <- function(r_y, r_a, P)
    (P$v_a * r_y^2 - 2 * P$c_ya * r_y * r_a + P$v_y * r_a^2) / P$det_v

  stat_cov <- function(r_y, r_a, slope_y, slope_a)
    quad(r_y, r_a, V_parts(slope_y, slope_a))

  V_obs <- V_parts(sl_y, sl_a)
  T_obs <- quad(by - sl_y * bx, ba - sl_a * bx, V_obs)

  #--------------------------------------------------------------------------#
  # Parametric bootstrap null, following MR-PRESSO's construction:
  # outcomes are drawn about the OBSERVED leave-one-out fit, and the residual is
  # then formed against the RESAMPLED exposure.
  #
  # Two null quantities, because mr_presso builds its two tests differently:
  #   per-SNP  : residual formed with the OBSERVED leave-one-out slope
  #              (mr_presso.R, `Exp <- randomSNP[,Y] - randomSNP[,X]*RSSobs[[2]][SNV]`)
  #   global   : leave-one-out slopes REFITTED on each bootstrap dataset
  #              (mr_presso.R, `RSSexp <- sapply(randomData, getRSS_LOO, ...)`)
  # Slopes are refitted to match mr_presso's construction.
  #
  # THREE independent standard normals per draw, one per GWAS.  The exposure is
  # assumed independent of both outcomes; the two outcomes are coupled through
  # rho by the standard Cholesky construction (rho z_y + sqrt(1-rho^2) z_a),
  # which gives corr(e_y, e_a) = rho exactly.  The shared exposure error enters
  # separately, through the single bx_b used in BOTH residuals -- so neither
  # dependence source is double counted.
  #--------------------------------------------------------------------------#
  ## The NbDistribution x K matrix is allocated ONLY on the bootstrap route.
  ## At K = 1270 and NbDistribution = 10000 it is 127 MB, which is the single
  ## largest allocation in the function.
  T_null   <- if (use_boot_snp) matrix(NA_real_, NbDistribution, K) else NULL
  RSS_null <- numeric(NbDistribution)               # global : refitted slopes

  rho_scale <- sqrt(1 - rho^2)

  for (b in seq_len(NbDistribution)) {

    z_x <- stats::rnorm(K)
    z_y <- stats::rnorm(K)
    z_a <- stats::rnorm(K)

    bx_b <- bx + sx * z_x
    by_b <- sl_y * bx + sy * z_y
    ba_b <- sl_a * bx + sa * (rho * z_y + rho_scale * z_a)

    ## Per-SNP null: observed slopes, hence V_obs, fixed across draws.
    if (use_boot_snp)
      T_null[b, ] <- quad(by_b - sl_y * bx_b, ba_b - sl_a * bx_b, V_obs)

    ## Global null: refit both slope vectors, and rebuild V from them so the
    ## null dataset is treated exactly as the observed one was.
    sl_y_b <- loo(bx_b, by_b)
    sl_a_b <- loo(bx_b, ba_b)

    RSS_null[b] <- sum(stat_cov(by_b - sl_y_b * bx_b,
                                ba_b - sl_a_b * bx_b,
                                sl_y_b, sl_a_b))
  }

  ## GLOBAL test: bootstrap, strict '>' as in mr_presso.  Not analytic, because
  ## the leave-one-out slopes are refitted on every draw.
  p_global <- mean(RSS_null > sum(T_obs), na.rm = TRUE)

  #--------------------------------------------------------------------------#
  # PER-SNP p-values.  Both are always computed when available, and both are
  # returned; OutlierTestMethod selects which one FLAGS.
  #
  # Bonferroni across instruments in both cases, exactly as MR-PRESSO does
  # (mr_presso.R, line 83), so the two are on the same scale and directly
  # comparable to MR-PRESSO's.
  #--------------------------------------------------------------------------#
  p_snp_chisq <- pmin(stats::pchisq(T_obs, df = 2, lower.tail = FALSE) * K, 1)

  p_snp_boot <- if (use_boot_snp) {
    pmin(colMeans(sweep(T_null, 2, T_obs, ">"), na.rm = TRUE) * K, 1)
  } else {
    rep(NA_real_, K)
  }

  p_snp <- if (use_boot_snp) p_snp_boot else p_snp_chisq

  #--------------------------------------------------------------------------#
  # Causal estimate.  lm() with MR-PRESSO's weights, so the standard error
  # follows the same convention (residual scale included) and the two methods
  # remain comparable on type-I error.
  #--------------------------------------------------------------------------#
  form     <- stats::as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure))
  fit_full <- stats::lm(form, weights = d$.ibmr_w1, data = d)
  beta_raw <- unname(stats::coef(fit_full)[1])
  se_raw   <- unname(summary(fit_full)$coefficients[1, "Std. Error"])

  outlier_indices <- integer(0)
  beta_corrected  <- NA_real_
  se_corrected    <- NA_real_

  if (OUTLIERtest && p_global < SignifThreshold) {
    outlier_indices <- which(p_snp <= SignifThreshold)
    keep <- setdiff(seq_len(K), outlier_indices)
    if (length(outlier_indices) > 0 && length(keep) >= 3) {
      fit_trim       <- stats::lm(form, weights = d$.ibmr_w1[keep], data = d[keep, ])
      beta_corrected <- unname(stats::coef(fit_trim)[1])
      se_corrected   <- unname(summary(fit_trim)$coefficients[1, "Std. Error"])
    }
  }
  n_outliers <- length(outlier_indices)

  #--------------------------------------------------------------------------#
  # Distortion test
  #--------------------------------------------------------------------------#
  p_distortion <- NA_real_
  if (DISTORTIONtest && !is.na(beta_corrected) && n_outliers > 0) {
    bias_observed <- (beta_raw - beta_corrected) / abs(beta_corrected)
    bias_null <- replicate(NbDistribution, {
      drop_i <- sample(seq_len(K), n_outliers)
      b_rand <- unname(stats::coef(stats::lm(form, weights = d$.ibmr_w1[-drop_i],
                                             data = d[-drop_i, ]))[1])
      (beta_raw - b_rand) / abs(b_rand)
    })
    p_distortion <- mean(abs(bias_null) >= abs(bias_observed), na.rm = TRUE)
  }

  results <- list(
    raw_beta       = beta_raw,
    raw_se         = se_raw,
    corrected_beta = beta_corrected,   # NA when nothing was removed
    corrected_se   = se_corrected,     # NA when nothing was removed
    p_value        = p_global,
    p_distort      = p_distortion,
    outlier_idx    = outlier_indices,
    outlyingness   = T_obs,
    outlier_pvals  = p_snp,
    n_instruments  = K,
    n_outliers     = n_outliers,

    ## ---- audit fields -----------------------------------------------------
    outcome_correlation = rho,
    ## Which route FLAGGED.  outlier_pvals above is a copy of whichever of the
    ## next two this names.
    outlier_test_method = OutlierTestMethod,
    ## Analytic per-SNP p-values.  Always computed.  T_obs is NOT claimed to be
    ## chi-square in general -- the claim is that its null, conditional on the
    ## observed leave-one-out slopes, is exactly chi^2_2.
    outlier_pvals_chisq = p_snp_chisq,
    ## Bootstrap per-SNP p-values.  All NA unless OutlierTestMethod =
    ## "bootstrap", since the NbDistribution x K matrix is not built otherwise.
    outlier_pvals_boot  = p_snp_boot,
    covariance_model    = paste0("analytic residual covariance: outcome LDSC ",
                                 "intercept + shared exposure error")
  )

  class(results) <- c("mrpresso_ib", "list")
  results
}
