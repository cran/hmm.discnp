\name{scovmat}
\alias{scovmat}
\title{
   Simulation based covariance matrix.
}
\description{
   Produces an estimate of the covariance matrix of the parameter
   estimates in a model fitted by \code{hmm.discnp}.  Uses a method
   based on simulation (or \dQuote{parametric bootstrapping}).
}
\usage{
scovmat(object, seed = NULL, nsim=100, verbose = TRUE)
}
\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by \code{\link{hmm}()}.
}
  \item{seed}{
  Integer scalar serving as a seed for the random number generator.
  If left \code{NULL} the seed itself is chosen randomly from the
  set of integers between 1 and \eqn{10^5}{1e5}.
}
  \item{nsim}{
  A positive integer.  The number of simulations upon which
  the covariance matrix estimate will be based.
}
  \item{verbose}{
  Logical scalar; if \code{TRUE}, iteration counts will be
  printed out during each of the simulation and model-fitting
  stages.
}
}
\details{
  This function is currently applicable only to models fitted to
  univariate data.  The covariance matrix produced is for the
  \dQuote{raw} parameters (entries of \code{tpm} with the
  last column dropped --- since the rows sum to 1, and the
  entries of \code{Rho} with the last row dropped --- since
  the columns sum to 1.
}
\value{
   A (positive definite) matrix which is an estimate of the
   covariance of the parameter estimates from the fitted model
   specified by \code{object}.  It has row and column labels
   which indicate the parameters to which its entries pertain,
   in a reasonably perspicuous manner.

   This matrix has an attribute \code{seed} (the random number
   generation seed that was used) so that the calculations can
   be reproduced.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}
\seealso{
\code{\link{squantCI}()} \code{link{rhmm}()} \code{link{hmm)}()}
}
\examples{
\dontrun{
y   <- list(lindLandFlows$deciles,ftLiardFlows$deciles)
fit <- hmm(y,K=3)
ccc <- scovmat(fit,nsim=100)
}
}
\keyword{ utility }
