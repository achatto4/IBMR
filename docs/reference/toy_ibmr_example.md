# Toy Summary-Statistics Example for IBMR

A small simulated summary-statistics example illustrating the intended
`IBMR` workflow with one exposure, one primary outcome, and two
candidate auxiliary outcomes.

## Usage

``` r
data(toy_ibmr_example)
```

## Format

A named list with the following entries:

- BetaXG:

  Numeric vector of SNP-exposure associations.

- seBetaXG:

  Numeric vector of standard errors for `BetaXG`.

- BetaYG_matrix:

  Matrix of SNP-outcome associations for the primary and auxiliary
  traits.

- seBetaYG_matrix:

  Matrix of standard errors matching `BetaYG_matrix`.

- primary_trait:

  Name of the primary outcome column.

- candidate_auxiliaries:

  Names of candidate auxiliary outcome columns.

- recommended_auxiliary:

  Auxiliary trait intended to rank highest in the toy example.

- dat_ibpresso_aux1:

  Data frame formatted for
  [`IBPRESSO()`](https://achatto4.github.io/IBMR/reference/IBPRESSO.md)
  using the recommended auxiliary trait.

- simulation_setup:

  List describing the simulation settings used to generate the toy data.

## Details

The data were simulated so that `aux_trait_1` shares a moderate latent
heterogeneity component with the primary outcome, whereas `aux_trait_2`
is much less aligned. This makes `aux_trait_1` the more natural
auxiliary trait for instrument borrowing in the example workflow.
