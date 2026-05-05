# Instrument Borrowing Workflow

## Overview

This vignette illustrates the full `IBMR` workflow:

1.  load SNP-level summary statistics
2.  screen candidate auxiliary traits using coheterogeneity
3.  select the most informative auxiliary trait
4.  use the selected auxiliary trait in
    [`IBMODE()`](https://achatto4.github.io/IBMR/reference/IBMODE.md)
5.  use the same auxiliary trait in
    [`IBPRESSO()`](https://achatto4.github.io/IBMR/reference/IBPRESSO.md)

The purpose is to demonstrate how auxiliary-trait selection and
downstream instrument borrowing fit together in a single analysis
pipeline.

## Load the Toy Example

``` r

data("toy_ibmr_example")

BetaXG <- toy_ibmr_example$BetaXG
seBetaXG <- toy_ibmr_example$seBetaXG
BetaYG_matrix <- toy_ibmr_example$BetaYG_matrix
seBetaYG_matrix <- toy_ibmr_example$seBetaYG_matrix
primary_name <- toy_ibmr_example$primary_trait
candidate_aux <- toy_ibmr_example$candidate_auxiliaries
```

The toy example contains one primary outcome and two candidate auxiliary
outcomes. The data were simulated so that `aux_trait_1` is usually the
more informative auxiliary trait for the primary outcome.

## Step 1: Screen Auxiliary Traits with Coheterogeneity

``` r

cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix,
  F_min = 5,
  min_K_pair = 20
)

round(cohet_res$rho, 3)
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait         1.000       0.622      -0.036
#> aux_trait_1           0.622       1.000      -0.141
#> aux_trait_2          -0.036      -0.141       1.000
round(cohet_res$p_value, 3)
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait         0.000           0       0.271
#> aux_trait_1           0.000           0       0.000
#> aux_trait_2           0.271           0       0.000
cohet_res$flag
#>               primary_trait aux_trait_1 aux_trait_2
#> primary_trait "diag"        "OK"        "OK"       
#> aux_trait_1   "OK"          "diag"      "OK"       
#> aux_trait_2   "OK"          "OK"        "diag"
```

We focus on the row corresponding to the primary outcome.

``` r

rho_primary <- cohet_res$rho[primary_name, ]
p_primary <- cohet_res$p_value[primary_name, ]
flag_primary <- cohet_res$flag[primary_name, ]

rho_primary
#> primary_trait   aux_trait_1   aux_trait_2 
#>    1.00000000    0.62170790   -0.03630623
p_primary
#> primary_trait   aux_trait_1   aux_trait_2 
#>  0.000000e+00 1.024886e-116  2.713419e-01
flag_primary
#> primary_trait   aux_trait_1   aux_trait_2 
#>        "diag"          "OK"          "OK"
```

## Step 2: Rank Candidate Auxiliary Traits

``` r

ranking <- data.frame(
  aux_trait = candidate_aux,
  rho = rho_primary[candidate_aux],
  p_value = p_primary[candidate_aux],
  flag = flag_primary[candidate_aux],
  stringsAsFactors = FALSE
)

ranking$abs_rho <- abs(ranking$rho)
ranking <- ranking[order(-ranking$abs_rho), ]
ranking
#>               aux_trait         rho       p_value flag    abs_rho
#> aux_trait_1 aux_trait_1  0.62170790 1.024886e-116   OK 0.62170790
#> aux_trait_2 aux_trait_2 -0.03630623  2.713419e-01   OK 0.03630623
```

The selected auxiliary trait can be taken as the candidate with the
largest absolute coheterogeneity among those with acceptable diagnostic
flags.

``` r

chosen_aux <- ranking$aux_trait[1]
chosen_aux
#> [1] "aux_trait_1"
```

For the packaged toy example, the intended selected auxiliary trait is:

``` r

toy_ibmr_example$recommended_auxiliary
#> [1] "aux_trait_1"
```

## Step 3: Run `IBMODE()`

Once the auxiliary trait has been selected, we subset the outcome
matrices to the primary outcome and the chosen auxiliary trait.

``` r

BetaYG_mode <- BetaYG_matrix[, c(primary_name, chosen_aux), drop = FALSE]
seBetaYG_mode <- seBetaYG_matrix[, c(primary_name, chosen_aux), drop = FALSE]
```

We then run
[`IBMODE()`](https://achatto4.github.io/IBMR/reference/IBMODE.md).

``` r

ibmode_res <- IBMODE(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_mode,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_mode,
  phi = c(1, 0.5),
  n_boot = 200
)

ibmode_res
```

This returns joint mode-based estimates for the primary and auxiliary
outcomes. In routine analyses, the number of bootstrap replicates should
be chosen to balance computational cost and stability.

## Step 4: Run `IBPRESSO()`

[`IBPRESSO()`](https://achatto4.github.io/IBMR/reference/IBPRESSO.md)
expects a data frame rather than matrices. The toy dataset includes a
ready-made example for the recommended auxiliary trait.

``` r

dat_ibpresso <- toy_ibmr_example$dat_ibpresso_aux1
head(dat_ibpresso)
#>   beta_exposure se_exposure beta_primary se_primary   beta_aux se_aux
#> 1    0.08879049       0.005  0.017732229      0.008 0.01965337  0.008
#> 2    0.09539645       0.005 -0.021832887      0.008 0.01928450  0.008
#> 3    0.13117417       0.005  0.015588900      0.008 0.02884544  0.008
#> 4    0.10141017       0.005  0.003214524      0.008 0.04533534  0.008
#> 5    0.10258575       0.005  0.093197322      0.008 0.08665289  0.008
#> 6    0.13430130       0.005  0.012196435      0.008 0.01398760  0.008
```

The analysis can be run as follows.

``` r

ibpresso_res <- IBPRESSO(
  BetaOutcome = "beta_primary",
  BetaExposure = "beta_exposure",
  BetaAux = "beta_aux",
  SdOutcome = "se_primary",
  SdExposure = "se_exposure",
  SdAux = "se_aux",
  data = dat_ibpresso,
  OUTLIERtest = TRUE,
  DISTORTIONtest = FALSE,
  NbDistribution = 200,
  seed = 123
)

ibpresso_res
```

This analysis uses the selected auxiliary trait to improve outlier
detection in the joint primary-auxiliary residual space.

## Interpreting the Workflow

The key logic of the instrument-borrowing approach is:

1.  identify a related auxiliary trait with meaningful coheterogeneity
2.  preserve the primary outcome as the principal target of inference
3.  use the auxiliary trait only to stabilize downstream robust MR
    analysis

The auxiliary trait is therefore not an alternative endpoint. Instead,
it is an informative companion trait used to strengthen inference for
the primary outcome.

## Practical Recommendations

- Ensure that SNP order is aligned across all inputs.
- Harmonize alleles before constructing summary-statistics matrices.
- Examine `flag` and `K` from
  [`coheterogeneity_Q()`](https://achatto4.github.io/IBMR/reference/coheterogeneity_Q.md)
  before selecting an auxiliary trait.
- Use more than one candidate auxiliary trait whenever scientifically
  reasonable, since screening is most useful when several plausible
  auxiliaries are available.
- Increase `n_boot` in
  [`IBMODE()`](https://achatto4.github.io/IBMR/reference/IBMODE.md) and
  `NbDistribution` in
  [`IBPRESSO()`](https://achatto4.github.io/IBMR/reference/IBPRESSO.md)
  for final analyses.

## Summary

`IBMR` supports a two-stage strategy:

- first, identify a suitable auxiliary trait using coheterogeneity
- second, use that auxiliary trait to improve robust MR estimation and
  outlier detection

This framework is particularly useful when related outcome traits are
expected to share informative instrument-validity structure relative to
a common exposure.
