#' Recover MR-PRESSO's outlier-corrected estimate
#'
#' `MRPRESSO::mr_presso` fits its outlier-corrected model only inside the branch
#'
#' ```
#' if (DISTORTIONtest & OUTLIERtest) { ... mod_noOutliers <- lm(...) ... }
#' ```
#'
#' and builds the reported `"Outlier-corrected"` row from
#' `exists("mod_noOutliers")`. With `DISTORTIONtest = FALSE` that object is never
#' created, so the corrected row is `NA` **even when instruments were flagged**,
#' and the function emits the misleading warning `"No outlier were identified"`.
#' Callers that fall back to the `"Raw"` row therefore report an inverse-variance
#' weighted fit that still contains every instrument the outlier test rejected.
#'
#' This function takes MR-PRESSO's own outlier list and refits, using the same
#' weighted IVW through the origin that [IBPRESSO()] uses for its corrected
#' estimate. The two methods then share a single refit code path, so the claim
#' that they compute an identical final estimator given the same retained set
#' holds by construction rather than by coincidence. It also avoids setting
#' `DISTORTIONtest = TRUE`, which would add a second full bootstrap.
#'
#' @section Threshold:
#' `Outlier Test$Pvalue` is **already** Bonferroni-corrected inside `mr_presso`
#' (`OutlierTest$Pvalue * nrow(data)`), so instruments are flagged at
#' `p <= SignifThreshold`, *not* `K * p <= SignifThreshold`. Applying the factor
#' again would double-correct. Exact zeros are returned by `mr_presso` as the
#' character string `"<K/NbDistribution"` and are parsed accordingly.
#'
#' @param presso The object returned by [MRPRESSO::mr_presso()].
#' @param data The same data frame that was passed to `mr_presso`.
#' @param BetaOutcome,BetaExposure,SdOutcome Column names, as passed to `mr_presso`.
#' @param SignifThreshold Significance threshold for the per-SNP outlier test.
#'
#' @return A list with `b`, `se`, `n_outliers`, and `corrected` (logical). The
#'   raw fit is returned when the global test was not significant (no
#'   `Outlier Test` element), when nothing was flagged, or when fewer than three
#'   instruments would survive.
#'
#' @importFrom stats lm coef as.formula
#' @export
mrpresso_refit <- function(presso, data, BetaOutcome, BetaExposure, SdOutcome,
                           SignifThreshold = 0.05) {
  ot   <- presso[["MR-PRESSO results"]][["Outlier Test"]]
  drop <- integer(0)
  if (!is.null(ot)) {
    ## CRITICAL.  mr_presso replaces an EXACTLY ZERO p-value with the STRING
    ##     paste0("<", nrow(data)/NbDistribution)
    ## so "<0.0455" does not mean "somewhere below 0.0455" -- it means the
    ## empirical p-value was 0, i.e. the observed statistic beat every null
    ## draw, the most extreme outcome possible.  Parsing that string as the
    ## numeric bound is wrong in principle, and once K/NbDistribution exceeds
    ## SignifThreshold it silently discards precisely the MOST significant
    ## instruments -- at K = 1270 and NbDistribution = 10000 the bound is
    ## 0.127, so every exact zero would be treated as non-significant.
    p_chr   <- as.character(ot$Pvalue)
    is_zero <- startsWith(p_chr, "<")
    p_num   <- suppressWarnings(as.numeric(p_chr))
    flagged <- which(is_zero | (!is.na(p_num) & p_num <= SignifThreshold))
    if (length(flagged)) {
      ## mr_presso sets row.names(OutlierTest) <- row.names(data) after its
      ## complete-case filter, so map by name; if any name fails to match
      ## (non-standard row names on the input), fall back to position.
      idx <- match(rownames(ot)[flagged], rownames(data))
      drop <- if (anyNA(idx)) flagged else idx
    }
  }
  keep <- setdiff(seq_len(nrow(data)), drop)
  used_correction <- length(drop) > 0 && length(keep) >= 3
  if (!used_correction) keep <- seq_len(nrow(data))

  d <- data[keep, , drop = FALSE]
  d$.mrp_w <- 1 / d[[SdOutcome]]^2
  fit <- lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)),
            weights = d$.mrp_w, data = d)
  list(b          = unname(coef(fit)[1]),
       se         = unname(summary(fit)$coefficients[1, "Std. Error"]),
       n_outliers = if (used_correction) length(drop) else 0L,
       corrected  = used_correction)
}
