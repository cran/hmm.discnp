\name{fitted.hmm.discnp}
\Rdversion{1.1}
\alias{fitted.hmm.discnp}
\title{
  Fitted values of a discrete non-parametric hidden Markov
  model.
}
\description{
  Calculates the fitted values of a discrete non-parametric hidden
  Markov model as the conditional expectations of the observations,
  given the entire observation sequence (and the estimated parameters
  of the model).
}
\usage{
\method{fitted}{hmm.discnp}(object, ...)
}
\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{\dots}{
  Not used.
}
}
\details{
  This function is essentially a wrapper for \code{\link{sp}()}.
  This latter function gets called with \code{mean=TRUE}.
  The fitted values are the conditional expected values of the
  observations, given the entire observation sequence.  Note that the
  \dQuote{conditional expected value} concept makes sense \emph{only}
  if the observations are \emph{numeric}. If this is not the case,
  an error is thrown.  The observation sequence must be present in
  \code{object} (which will be the case if \code{object} was returned
  by \code{hmm()}).  If it is not present an error is thrown.
}
\value{
  If the observations (stored in \code{object}) consist of a
  single sequence, then the object returned is the corresponding
  sequence of fitted values.  If the observations consist of a
  list of observation sequences, then the object returned is a
  list of corresponding sequences of fitted values.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
  \url{http://www.stat.auckland.ac.nz/~rolf}
}
\seealso{ \code{\link{sp}()} }
\examples{
# See the help for sim.hmm() for how to generate y.num.
\dontrun{
fit.num <- hmm(y.num,K=2,verb=TRUE)
fv      <- fitted(fit.num)
}
}

\keyword{models}
