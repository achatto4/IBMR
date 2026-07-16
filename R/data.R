#' Illustrative Instrument-Borrowing Example
#'
#' A simulated summary-statistics example for the `IBMR` workflow with one
#' exposure, one primary outcome, and two candidate auxiliary outcomes. The data
#' are simulated (not real GWAS data): the exposure has a genuine positive effect
#' on the primary outcome, `Auxiliary_1` shares pleiotropic (invalid-instrument)
#' structure with the primary outcome and is therefore informative for instrument
#' borrowing, while `Auxiliary_2` shares little and is not.
#'
#' @format A named list with the following entries:
#' \describe{
#'   \item{BetaXG}{Numeric vector of SNP-exposure associations.}
#'   \item{seBetaXG}{Standard errors for `BetaXG`.}
#'   \item{BetaYG_matrix}{Matrix of SNP-outcome associations; columns `Primary`,
#'   `Auxiliary_1`, and `Auxiliary_2`.}
#'   \item{seBetaYG_matrix}{Standard errors matching `BetaYG_matrix`.}
#'   \item{primary_trait}{`"Primary"`.}
#'   \item{candidate_auxiliaries}{`c("Auxiliary_1", "Auxiliary_2")`.}
#'   \item{recommended_auxiliary}{`"Auxiliary_1"`.}
#'   \item{dat_ibpresso_aux1}{Data frame formatted for `IBPRESSO()` with the
#'   primary outcome and the recommended auxiliary outcome (`Auxiliary_1`).}
#'   \item{source_note}{Description of how the example was generated.}
#' }
#'
#' @usage data(ibmr_example)
#'
"ibmr_example"
