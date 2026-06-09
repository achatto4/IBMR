# Calculate Coheterogeneity with Fixed- or Random-Weight Standard Errors

Computes pairwise coheterogeneity correlations across multiple outcome
traits using the bias-corrected moment estimator for
instrument-borrowing Mendelian randomization. The standard error of the
estimate can be obtained two ways:

## Usage

``` r
coheterogeneity_Q(
  BetaXG,
  BetaYG_matrix = NULL,
  seBetaXG,
  seBetaYG_matrix = NULL,
  BetaIV_matrix = NULL,
  seBetaIV_matrix = NULL,
  ldsc_intercepts = NULL,
  SNP_keep = NULL,
  use_ldsc = TRUE,
  se_weights = c("fixed", "random", "both"),
  grad_rel_step = 0.001,
  alpha = 0.05,
  eps = 1e-12,
  bx_min = 1e-06,
  F_min = 0,
  min_K_pair = 50,
  return_diagnostics = FALSE
)
```

## Arguments

- BetaXG:

  Numeric vector of SNP-exposure associations.

- BetaYG_matrix:

  Matrix of SNP-outcome associations, rows = SNPs, columns = traits.
  Optional if `BetaIV_matrix` is provided.

- seBetaXG:

  Numeric vector of standard errors for `BetaXG`.

- seBetaYG_matrix:

  Matrix of standard errors for `BetaYG_matrix`. Optional if
  `seBetaIV_matrix` is provided.

- BetaIV_matrix:

  Matrix of ratio estimates (`BetaYG / BetaXG`), rows = SNPs, columns =
  traits. Used to reconstruct `BetaYG_matrix` when needed.

- seBetaIV_matrix:

  Matrix of standard errors for `BetaIV_matrix`. Used to reconstruct
  `seBetaYG_matrix` when needed.

- ldsc_intercepts:

  Optional matrix of LDSC intercepts used for pairwise outcome
  covariance adjustment.

- SNP_keep:

  Optional logical or integer index vector specifying the SNPs to
  retain.

- use_ldsc:

  Logical; if `TRUE`, use `ldsc_intercepts` when supplied.

- se_weights:

  Character; how the inverse-variance weights are treated when computing
  the standard error. One of `"fixed"` (closed-form plug-in SE, weights
  held fixed; the default), `"random"` (numerical delta-method SE that
  propagates the data-dependence of the weights), or `"both"` (return
  both).

- grad_rel_step:

  Numeric; relative finite-difference step used by the random-weight
  numerical SE, expressed as a multiple of each association's standard
  error. Only used when `se_weights` is `"random"` or `"both"`.

- alpha:

  Significance level retained in the returned guard settings.

- eps:

  Small positive constant used for numerical stabilization.

- bx_min:

  Minimum absolute exposure effect size allowed.

- F_min:

  Minimum first-stage F statistic threshold, where
  `F = (BetaXG / seBetaXG)^2`.

- min_K_pair:

  Minimum number of valid SNPs required for a trait pair.

- return_diagnostics:

  Logical; if `TRUE`, include a diagnostics list.

## Value

A list containing:

- rho:

  Pairwise coheterogeneity correlation matrix.

- se:

  Pairwise standard error matrix for `rho`, computed with the primary
  weighting (`"fixed"` for `se_weights` `"fixed"` or `"both"`,
  `"random"` for `se_weights` `"random"`).

- z_statistic:

  Pairwise z-statistic matrix for the primary SE.

- wald_statistic:

  Pairwise Wald statistic matrix for the primary SE.

- p_value:

  Pairwise p-value matrix for the primary SE.

- se_fixed, z_fixed, wald_fixed, p_value_fixed:

  Fixed-weight (closed-form) results; returned only when
  `se_weights = "both"`.

- se_random, z_random, wald_random, p_value_random:

  Random-weight (numerical) results; returned only when
  `se_weights = "both"`.

- K:

  Matrix of valid SNP counts per trait pair.

- flag:

  Matrix of per-pair diagnostic flags.

- method:

  Method label for the estimator.

- se_weights:

  The weighting option used for the standard error.

- guards:

  List of guardrail settings used in the analysis.

- diagnostics:

  Optional detailed diagnostics by trait pair.

## Details

- `"fixed"`:

  The closed-form plug-in SE (Theorem 1). The inverse-variance weights
  are treated as fixed (non-random) when the estimator is linearized.
  This is the original, computationally cheap SE.

- `"random"`:

  A numerical delta-method SE that differentiates the estimator through
  the data-dependent inverse-variance weights, i.e. the weights are
  treated as random (estimated). This is the exact first-order SE; the
  closed form omits the weight-estimation term.

The two coincide when the per-SNP ratios are homogeneous and diverge as
the coheterogeneity dispersion grows, the fixed-weight SE being the
conservative (larger) of the two. Use `se_weights = "both"` to return
them together.
