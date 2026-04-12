#' Toy Summary-Statistics Example for IBMR
#'
#' A small simulated summary-statistics example illustrating the intended
#' `IBMR` workflow with one exposure, one primary outcome, and two candidate
#' auxiliary outcomes.
#'
#' The data were simulated so that `aux_trait_1` shares more structure with the
#' primary outcome than `aux_trait_2`, making it the more natural auxiliary
#' trait for instrument borrowing in the example workflow.
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
