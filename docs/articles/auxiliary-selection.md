# Auxiliary Trait Selection with Coheterogeneity

## Overview

This vignette describes the first step in the `IBMR` workflow: selecting
an auxiliary outcome trait for a given exposure and primary outcome.

The central motivation is that closely related outcomes may share
overlapping valid instruments, or may reflect similar mechanisms of
instrument invalidity, with respect to a common exposure. If this shared
structure can be identified, an auxiliary trait may be used to improve
downstream MR estimation and outlier detection.

`IBMR` operationalizes this idea through
[`coheterogeneity_Q()`](https://achatto4.github.io/IBMR/reference/coheterogeneity_Q.md),
which computes pairwise coheterogeneity between outcome traits using
SNP-level summary statistics.

## The Data Structure

At minimum, the coheterogeneity analysis requires:

- `BetaXG`: SNP-exposure associations
- `seBetaXG`: standard errors for `BetaXG`
- `BetaYG_matrix`: a matrix of SNP-outcome associations
- `seBetaYG_matrix`: a matrix of standard errors corresponding to
  `BetaYG_matrix`

Each row corresponds to a SNP, and the rows must be aligned across all
vectors and matrices.

In the common screening setting, the columns of `BetaYG_matrix` are:

- one primary outcome
- one or more candidate auxiliary outcomes

## Toy Example Data

The package includes a toy dataset called `toy_ibmr_example`.

``` r
data("toy_ibmr_example")

names(toy_ibmr_example)
#> [1] "BetaXG"                "seBetaXG"              "BetaYG_matrix"        
#> [4] "seBetaYG_matrix"       "primary_trait"         "candidate_auxiliaries"
#> [7] "recommended_auxiliary" "dat_ibpresso_aux1"     "simulation_setup"
toy_ibmr_example$primary_trait
#> [1] "primary_trait"
toy_ibmr_example$candidate_auxiliaries
#> [1] "aux_trait_1" "aux_trait_2"
```

The object contains:

- one exposure vector
- one primary outcome
- two candidate auxiliary traits
- a small formatted data frame for
  [`IBPRESSO()`](https://achatto4.github.io/IBMR/reference/IBPRESSO.md)
- metadata describing the simulation setup

The simulation is deliberately simple. It was constructed so that
`aux_trait_1` shares more latent structure with the primary outcome than
`aux_trait_2`, making it the more natural auxiliary trait in this
example.

``` r
str(toy_ibmr_example$simulation_setup)
#> List of 11
#>  $ seed                        : num 123
#>  $ n_snps                      : num 80
#>  $ exposure_mean               : num 0.08
#>  $ exposure_sd                 : num 0.02
#>  $ exposure_se                 : num 0.01
#>  $ outcome_se                  : num 0.02
#>  $ primary_effect_on_exposure  : num 0.2
#>  $ aux1_effect_on_exposure     : num 0.25
#>  $ aux1_shared_component_weight: num 0.8
#>  $ aux2_effect_on_exposure     : num -0.05
#>  $ shared_component_sd         : num 0.015
```

## Running the Coheterogeneity Analysis

We begin by extracting the summary-statistics objects.

``` r
BetaXG <- toy_ibmr_example$BetaXG
seBetaXG <- toy_ibmr_example$seBetaXG
BetaYG_matrix <- toy_ibmr_example$BetaYG_matrix
seBetaYG_matrix <- toy_ibmr_example$seBetaYG_matrix

colnames(BetaYG_matrix)
#> [1] "primary_trait" "aux_trait_1"   "aux_trait_2"
```

We now compute coheterogeneity across the primary outcome and the
candidate auxiliary traits.

``` r
cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix,
  F_min = 5,
  min_K_pair = 20
)
```

The principal outputs are:

- `rho`: pairwise coheterogeneity estimates
- `p_value`: pairwise significance levels
- `flag`: diagnostic status for each pair
- `K`: number of SNPs contributing to each pairwise calculation

``` r
round(cohet_res$rho, 3)
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait             1          NA          NA
#> aux_trait_1              NA           1          NA
#> aux_trait_2              NA          NA           1
round(cohet_res$p_value, 3)
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait             0          NA          NA
#> aux_trait_1              NA           0          NA
#> aux_trait_2              NA          NA           0
cohet_res$flag
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait "diag"        "tau0"      "tau0"     
#> aux_trait_1   "tau0"        "diag"      "tau0"     
#> aux_trait_2   "tau0"        "tau0"      "diag"
cohet_res$K
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait             0          80          80
#> aux_trait_1              80           0          80
#> aux_trait_2              80          80           0
```

## Interpreting the Output

The diagonal entries of `rho` are equal to 1 by definition. The
quantities of primary interest are the off-diagonal entries involving
the primary outcome.

``` r
primary_name <- toy_ibmr_example$primary_trait
cohet_res$rho[primary_name, ]
#> primary_trait   aux_trait_1   aux_trait_2 
#>             1            NA            NA
cohet_res$p_value[primary_name, ]
#> primary_trait   aux_trait_1   aux_trait_2 
#>             0            NA            NA
cohet_res$flag[primary_name, ]
#> primary_trait   aux_trait_1   aux_trait_2 
#>        "diag"        "tau0"        "tau0"
```

In this toy example, `aux_trait_1` should usually exhibit stronger
coheterogeneity with the primary outcome than `aux_trait_2`.

A practical interpretation framework is:

- larger `abs(rho)` indicates a stronger shared heterogeneity pattern
- smaller `p_value` provides stronger evidence against no
  coheterogeneity
- `flag` should be checked to ensure the pairwise estimate is usable

The most important diagnostic flags are:

- `"OK"`: the pairwise comparison was successfully estimated
- `"tau0"`: the corrected heterogeneity estimate was not identifiable
- `"se NA"`: the coheterogeneity estimate was computed, but its standard
  error could not be reliably estimated

## Ranking Candidate Auxiliary Traits

The screening step is usually carried out by ranking candidate
auxiliaries for the primary outcome.

``` r
ranking <- data.frame(
  aux_trait = colnames(cohet_res$rho),
  rho = cohet_res$rho[primary_name, ],
  p_value = cohet_res$p_value[primary_name, ],
  flag = cohet_res$flag[primary_name, ],
  stringsAsFactors = FALSE
)

ranking <- subset(ranking, aux_trait != primary_name)
ranking$abs_rho <- abs(ranking$rho)
ranking <- ranking[order(-ranking$abs_rho), ]
ranking
#>               aux_trait rho p_value flag abs_rho
#> aux_trait_1 aux_trait_1  NA      NA tau0      NA
#> aux_trait_2 aux_trait_2  NA      NA tau0      NA
```

A practical selection rule is:

1.  discard pairs with unusable flags
2.  prioritize pairs with finite and informative `p_value`
3.  among the remaining candidates, select the trait with the largest
    `abs(rho)`

In this toy dataset, the intended selected auxiliary trait is:

``` r
toy_ibmr_example$recommended_auxiliary
#> [1] "aux_trait_1"
```

## Notes on Guardrails

[`coheterogeneity_Q()`](https://achatto4.github.io/IBMR/reference/coheterogeneity_Q.md)
includes several guardrails to stabilize estimation:

- optional SNP restriction through `SNP_keep`
- weak-instrument filtering through `F_min`
- small-exposure filtering through `bx_min`
- optional winsorization through `winsor_theta`
- pair-level diagnostic flags

These are especially useful in large-scale screening settings where some
traits may have sparse or unstable SNP-level signal.

## Next Step

Once an auxiliary trait has been selected, it can be carried into
downstream instrument-borrowing analysis using:

- [`IBMODE()`](https://achatto4.github.io/IBMR/reference/IBMODE.md) for
  multidimensional mode-based estimation
- [`IBPRESSO()`](https://achatto4.github.io/IBMR/reference/IBPRESSO.md)
  for outlier detection and correction with an auxiliary trait

The end-to-end downstream workflow is illustrated in the companion
vignette:

- `instrument-borrowing-workflow`
