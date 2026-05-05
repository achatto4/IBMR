% Generated manually to match R/data.R
\name{toy_ibmr_example}
\docType{data}
\title{Toy Summary-Statistics Example for IBMR}
\format{
A named list with the following entries:
\describe{
  \item{BetaXG}{Numeric vector of SNP-exposure associations.}
  \item{seBetaXG}{Numeric vector of standard errors for \code{BetaXG}.}
  \item{BetaYG_matrix}{Matrix of SNP-outcome associations for the primary and auxiliary traits.}
  \item{seBetaYG_matrix}{Matrix of standard errors matching \code{BetaYG_matrix}.}
  \item{primary_trait}{Name of the primary outcome column.}
  \item{candidate_auxiliaries}{Names of candidate auxiliary outcome columns.}
  \item{recommended_auxiliary}{Auxiliary trait intended to rank highest in the toy example.}
  \item{dat_ibpresso_aux1}{Data frame formatted for \code{IBPRESSO()} using the recommended auxiliary trait.}
  \item{simulation_setup}{List describing the simulation settings used to generate the toy data.}
}
}
\usage{data(toy_ibmr_example)}
\description{
A small simulated summary-statistics example illustrating the intended
\code{IBMR} workflow with one exposure, one primary outcome, and two candidate
auxiliary outcomes.
}
\details{
The data were simulated so that \code{aux_trait_1} shares a moderate latent
heterogeneity component with the primary outcome, whereas \code{aux_trait_2} is
much less aligned. This makes \code{aux_trait_1} the more natural auxiliary
trait for instrument borrowing in the example workflow.
}
\keyword{datasets}
