#' Toy Summary-Statistics Example for IBMR
#'
#' A small simulated summary-statistics example illustrating the intended
#' `IBMR` workflow with one exposure, one primary outcome, and two candidate
#' auxiliary outcomes.
#'
#' The data were simulated so that `aux_trait_1` shares a moderate latent
#' heterogeneity component with the primary outcome, whereas `aux_trait_2` is
#' much less aligned. This makes `aux_trait_1` the more natural auxiliary trait
#' for instrument borrowing in the example workflow.
#'
#' @format A named list with the following entries:
#' \describe{
#'   \item{BetaXG}{Numeric vector of SNP-exposure associations.}
#'   \item{seBetaXG}{Numeric vector of standard errors for `BetaXG`.}
#'   \item{BetaYG_matrix}{Matrix of SNP-outcome associations for the primary and
#'   auxiliary traits.}
#'   \item{seBetaYG_matrix}{Matrix of standard errors matching
#'   `BetaYG_matrix`.}
#'   \item{primary_trait}{Name of the primary outcome column.}
#'   \item{candidate_auxiliaries}{Names of candidate auxiliary outcome columns.}
#'   \item{recommended_auxiliary}{Auxiliary trait intended to rank highest in
#'   the toy example.}
#'   \item{dat_ibpresso_aux1}{Data frame formatted for `IBPRESSO()` using the
#'   recommended auxiliary trait.}
#'   \item{simulation_setup}{List describing the simulation settings used to
#'   generate the toy data.}
#' }
#'
#' @usage data(toy_ibmr_example)
#'
"toy_ibmr_example"


#' Illustrative Instrument-Borrowing Example
#'
#' A simulated summary-statistics example for the `IBMR` workflow, framed as the
#' effect of body mass index (BMI) on coronary artery disease (CAD) with type 2
#' diabetes (T2D) as the auxiliary trait. The data are simulated (not real GWAS
#' data): the exposure has a genuine positive effect on both outcomes, and a
#' subset of invalid instruments exhibit pleiotropy shared across CAD and T2D,
#' so T2D is an informative auxiliary trait for instrument borrowing. Trait names
#' are used only to make the workflow concrete.
#'
#' @format A named list with the following entries:
#' \describe{
#'   \item{BetaXG}{Numeric vector of SNP-exposure associations.}
#'   \item{seBetaXG}{Standard errors for `BetaXG`.}
#'   \item{BetaYG_matrix}{Matrix of SNP-outcome associations; columns `CAD`
#'   (primary) and `T2D` (auxiliary).}
#'   \item{seBetaYG_matrix}{Standard errors matching `BetaYG_matrix`.}
#'   \item{primary_trait}{`"CAD"`.}
#'   \item{candidate_auxiliaries}{`"T2D"`.}
#'   \item{recommended_auxiliary}{`"T2D"`.}
#'   \item{dat_ibpresso_aux1}{Data frame formatted for `IBPRESSO()` with CAD as
#'   the primary outcome and T2D as the borrowed auxiliary trait.}
#'   \item{source_note}{Description of how the example was generated.}
#' }
#'
#' @usage data(ibmr_example)
#'
"ibmr_example"
