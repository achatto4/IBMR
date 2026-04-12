# Multidimensional Mode-Based Estimation for Mendelian Randomization

Performs mode-based estimation (MBE) for multiple correlated outcomes
simultaneously using multidimensional kernel density estimation. This
approach identifies the modal causal effect across multiple traits while
accounting for pleiotropy.

## Usage

``` r
IBMODE(
  BetaXG,
  BetaYG_matrix,
  seBetaXG,
  seBetaYG_matrix,
  phi = c(1, 0.5, 0.25),
  n_boot = 10000,
  alpha = 0.05
)
```

## Arguments

- BetaXG:

  Numeric vector of SNP-exposure associations

- BetaYG_matrix:

  Matrix of SNP-outcome associations (rows = SNPs, columns =
  outcomes/traits)

- seBetaXG:

  Numeric vector of standard errors for BetaXG

- seBetaYG_matrix:

  Matrix of standard errors for BetaYG_matrix

- phi:

  Numeric vector of bandwidth multipliers for sensitivity analysis
  (default: c(1, 0.5, 0.25))

- n_boot:

  Integer, number of bootstrap iterations for standard error estimation
  (default: 10000)

- alpha:

  Numeric, significance level for confidence intervals (default: 0.05)

## Value

A data frame containing:

- Method:

  Method name

- phi:

  Bandwidth multiplier used

- Estimate.1, Estimate.2, ...:

  Causal effect estimates for each outcome

- SE.1, SE.2, ...:

  Standard errors for each outcome

- CI_low.1, CI_low.2, ...:

  Lower confidence interval bounds

- CI_upp.1, CI_upp.2, ...:

  Upper confidence interval bounds

- P.1, P.2, ...:

  P-values for each outcome

## Details

The function uses weighted multidimensional kernel density estimation to
find the mode of the joint distribution of causal effects across
multiple outcomes. Bootstrap resampling is used to estimate standard
errors. Multiple bandwidth parameters (phi) can be tested for
sensitivity analysis.
