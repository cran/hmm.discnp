\name{squantCI}
\alias{squantCI}
\title{
   Simulation-quantile based confidence intervals.
}
\description{
   Calculates estimates of confidence intervals for the parameters of a
   model fitted by \code{hmm.discnp}.  Uses a method based quantiles
   of estimates produced by simulation (or \dQuote{parametric
   bootstrapping}).
}
\usage{
squantCI(object, seed = NULL, alpha = 0.05, nsim=100, verbose = TRUE)
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
  \item{alpha}{
  Positive real number strictly between 0 and 1.  A set of
  \code{100*(1-alpha)}\% confidence intervals will be produced.
}
  \item{nsim}{
  A positive integer.  The number of simulations upon which
  the confidence interval estimates will be based.
}
  \item{verbose}{
  Logical scalar; if \code{TRUE}, iteration counts will be
  printed out during each of the simulation and model-fitting
  stages.
}
}
\details{
  This function is currently applicable only to models fitted to
  univariate data.  The confidence intervals calculated are for the
  \dQuote{raw} parameters (entries of \code{tpm} with the
  last column dropped --- since the rows sum to 1, and the
  entries of \code{Rho} with the last row dropped --- since
  the columns sum to 1.
}
\value{
   A 2-by-\code{npar} matrix (where \code{npar} is the number of
   \dQuote{independent} parameters in the model) whose columns
   form the estimated confidence intervals.  The column labels
   indicate the parameters to which each column pertains, in a
   reasonably perspicuous manner.  The row labels indicate the
   relevant quantiles in percentages.

   This matrix has an attribute \code{seed} (the random number
   generation seed that was used) so that the calculations can
   be reproduced.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}
\seealso{
\code{\link{scovmat}()} \code{link{rhmm}()} \code{link{hmm)}()}
}
\examples{
\dontrun{
y   <- list(lindLandFlows$deciles,ftLiardFlows$deciles)
fit <- hmm(y,K=3)
CIs <- squantCI(fit,nsim=100)
}
}
\keyword{ utility }
