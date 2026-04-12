# Calculate Coheterogeneity Using a Guarded Theoretical Method

Computes pairwise coheterogeneity correlations across multiple outcome
traits using a guarded delta-exact approach. The implementation supports
optional SNP filtering, weak-instrument filtering, winsorization of
ratio estimates, LDSC intercept adjustment, and per-pair diagnostics.

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
  alpha = 0.05,
  step_mult = 0.001,
  eps = 1e-12,
  bx_min = 1e-06,
  F_min = 0,
  winsor_theta = 20,
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

  Optional matrix of LDSC intercepts used for pairwise covariance
  adjustment.

- SNP_keep:

  Optional logical or integer index vector specifying the SNPs to
  retain.

- use_ldsc:

  Logical; if `TRUE`, use `ldsc_intercepts` when supplied.

- alpha:

  Significance level retained for interface compatibility.

- step_mult:

  Scalar used to scale finite-difference steps for the delta method
  standard error.

- eps:

  Small positive constant used for numerical stabilization.

- bx_min:

  Minimum absolute exposure effect size allowed.

- F_min:

  Minimum first-stage F statistic threshold, where
  `F = (BetaXG / seBetaXG)^2`.

- winsor_theta:

  Optional cap applied to ratio estimates to reduce instability from
  very small exposure effects. Set to `NULL` to disable.

- min_K_pair:

  Minimum number of valid SNPs required for a trait pair.

- return_diagnostics:

  Logical; if `TRUE`, include a diagnostics list.

## Value

A list containing:

- rho:

  Pairwise coheterogeneity correlation matrix.

- se:

  Pairwise standard error matrix for `rho`.

- z_statistic:

  Pairwise z-statistic matrix.

- p_value:

  Pairwise p-value matrix.

- K:

  Matrix of valid SNP counts per trait pair.

- flag:

  Matrix of per-pair diagnostic flags.

- method:

  Method label for the estimator.

- guards:

  List of guardrail settings used in the analysis.

- diagnostics:

  Optional detailed diagnostics by trait pair.
