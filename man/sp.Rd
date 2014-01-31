\name{sp}
\Rdversion{1.1}
\alias{sp}
\title{
  Calculate the conditional state probabilities.
}
\description{
  Returns the probabilities that the underlying hidden
  state is equal to each of the possible state values,
  at each time point, given the observation sequence.
  Also can return the fitted conditional means, if requested,
  given that the observations are numeric.
}
\usage{
sp(y, object = NULL, tpm, Rho, ispd=NULL, means = FALSE)
}
\arguments{
  \item{y}{The observations on the basis of which the probabilities
  of the underlying hidden states are to be calculated.  May be
  a sequence of observations, or a list each component of which
  constitutes a (replicate) sequence of observations.  If \code{y} is
  missing it is set equal to the \code{y} component of \code{object},
  given that that object and that component exist.  Otherwise an
  error is given.
}
  \item{object}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{tpm}{
  The transition probability matrix for the underlying hidden
  Markov chain.  Ignored if \code{object} is not \code{NULL}.
  Ignored if \code{object} is not \code{NULL} (in which
  case \code{tpm} is extracted from \code{object}).
}
  \item{Rho}{
  The matrix of probabilities specifying the distribution of
  the observations, given the underlying state.  The rows of this
  matrix correspond to the possible values of the observations, the
  columns to the states.  Ignored if \code{object} is not \code{NULL}
  (in which case \code{Rho} is extracted from \code{object}).
}
  \item{ispd}{
  Vector specifying the initial state probability distribution
  of the underlying hidden Markov chain.
  Ignored if \code{object} is not \code{NULL} (in which
  case \code{ispd} is extracted from \code{object}).
  If both \code{object} and \code{ispd} are NULL then
  \code{ispd} is calculated to be the stationary distribution
  of the chain as determined by \code{tpm}.
}
  \item{means}{
  A logical scalar; if \code{means} is \code{TRUE} then the
  conditional expected value of the observations (given the
  observation sequence) is calculated at each time point.
  If \code{means} is \code{TRUE} and the observation values
  are not numeric, then an error is given.
}
}
\details{
  Then conditional mean value at time \eqn{t} is calculated
  as \deqn{\sum_k \gamma_t(k) \mu_k}{SUM_k gamma_t(k)*mu_k}
  where \eqn{\gamma_t(k)}{gamma_t(k)} is the conditional
  probability (given the observations) that the hidden Markov
  chain is in state \eqn{k} at time \eqn{t}, and \eqn{\mu_k}{mu_k}
  is the expected value of an observation given that the chain
  is in state \eqn{k}.
}
\value{
    If \code{means} is \code{TRUE} then the returned value is
    a list with components
    \item{probs}{The conditional probabilities of the states
    at each time point.}
    \item{means}{The conditional expectations of the observations
                 at each time point.}

    Otherwise the returned value consists of \code{probs} as
    described above.

    If there is a single vector of observations \code{y} then
    \code{probs} is a matrix whose rows correspond to the states
    of the hidden Markov chain, and whose columns correspond to
    the observation times.  If the observations consist of a
    list of observation vectors, then \code{probs} is a list
    of such matrices, one for each vector of observations.

    Likewise for the \code{means} component of the list returned
    when the argument \code{means} is \code{TRUE}.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
  \url{http://www.stat.auckland.ac.nz/~rolf}
}

\seealso{\code{\link{hmm}()}, \code{\link{mps}()},
         \code{\link{viterbi}()}, \code{\link{pr}()},
         \code{\link{fitted.hmm.discnp}()}}
\examples{
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
y.num   <- sim.hmm(rep(300,20),P,R)
fit.num <- hmm(y.num,K=2,verb=TRUE)
cpe1    <- sp(object=fit.num,means=TRUE)    # Using the estimated parameters.
cpe2    <- sp(y.num,tpm=P,Rho=R,means=TRUE) # Using the ``true'' parameters.
}

\keyword{models}
