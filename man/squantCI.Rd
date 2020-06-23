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
squantCI(object, expForm = TRUE, seed = NULL, alpha = 0.05,
         nsim=100, verbose = TRUE)
}
\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by \code{\link{hmm}()}.
}
  \item{expForm}{
  Logical scalar.  Should the confidence intervals produced
  be for the parameters expressed in \dQuote{exponential}
  (or \dQuote{smooth} or \dQuote{logistic}) form?
  If \code{expForm=FALSE} then the parameters considered are
  \dQuote{raw} probabilities, with redundancies (last column of
  \code{tpm}; last row of \code{Rho}) removed.
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
  univariate data.  If there are \emph{predictors} in the model,
  then only the exponential form of the parameters may be used,
  i.e. \code{expForm} \emph{must} be \code{TRUE}.
}
\value{
   A \code{npar}-by-2 matrix (where \code{npar} is the number
   of \dQuote{independent} parameters in the model) whose rows
   form the estimated confidence intervals.  (The first entry of
   each row is the lower bound of a confidence interval for the
   corresponding parameter, and the second entry is the upper bound.
   The row labels indicate the parameters to which each row pertains,
   in a reasonably perspicuous manner.  The column labels indicate
   the relevant quantiles in percentages.

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
