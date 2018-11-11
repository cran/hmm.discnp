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
sp(y, model = NULL, tpm=NULL, Rho=NULL, ispd=NULL, X=NULL,
   addIntercept=NULL, means = FALSE, warn=TRUE)
}
\arguments{
  \item{y}{The observations on the basis of which the probabilities
  of the underlying hidden states are to be calculated.  May be
  a (one or two column) matrix of observations, or a list each
  component of which is such a matrix.  If \code{y} is missing it is
  set equal to the \code{y} component of \code{model}, given that
  that argument is non-\code{NULL} and that that component exist.
  Otherwise an error is given.
}
  \item{model}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{tpm}{
  The transition probability matrix for the underlying hidden
  Markov chain.  Ignored if \code{model} is not \code{NULL}.
  Ignored if \code{model} is not \code{NULL} (in which
  case \code{tpm} is extracted from \code{model}).
}
  \item{Rho}{
  An object specifying the distribution of the observations, given
  the underlying state.  I.e. the \dQuote{emission} probabilities.
  See \code{\link{hmm}()}.  Ignored if \code{model} is not
  \code{NULL} (in which case \code{Rho} is extracted from
  \code{model}).
}
  \item{ispd}{
  Vector specifying the initial state probability distribution
  of the underlying hidden Markov chain.
  Ignored if \code{model} is not \code{NULL} (in which
  case \code{ispd} is extracted from \code{model}).
  If both \code{model} and \code{ispd} are NULL then
  \code{ispd} is calculated to be the stationary distribution
  of the chain as determined by \code{tpm}.
}
  \item{X}{An optional \emph{numeric} matrix, or a list of such
  matrices, of \emph{predictors}.  The use of such predictors
  is (currently, at least) applicable only in the univariate
  emissions setting, and then only if \code{Rho} is supplied in
  the \dQuote{newstyle} format.  If \code{X} is a list it must
  be of the same length as \code{y} and all entries of this list
  must have the same number of columns.  The number of rows of each
  entry must be equal to the length of the corresponding entry of
  \code{y}. If \code{X} is a matrix then \code{y} should be a vector
  or one-column matrix (or a list with a single entry equal to such).
}
  \item{addIntercept}{
  Logical scalar.  See the documentation of \code{\link{hmm}()}.
  }
  \item{means}{
  A logical scalar; if \code{means} is \code{TRUE} then the
  conditional expected value of (each column of) the observations
  (given the observation sequence) is calculated at each time point.
  If \code{means} is \code{TRUE} and the observation values are
  not numeric, then an error is given.
}
  \item{warn}{Logical scalar; should a warning be issued if \code{Rho}
  hasn't got relevant dimension names?  (Note that if this is so,
  then the corresponding dimension names are formed from the sorted
  unique values of \code{y} or of the appropriate column(s) of \code{y}.
  And if \emph{this} is so, then the user should be sure that the
  ordering of the entries of \code{Rho} corresponds properly to the
  the sorted unique values of \code{y}.)  This argument is passed
  to the utility function \code{check.yval()} which actually issues
  the warning if \code{warn=TRUE}.
  }
}
\details{
  The conditional mean value at time \eqn{t} is calculated
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

    Clearly this makes sense only if the observations are numeric.

    Otherwise the returned value consists of \code{probs} as
    described above.

    If there is a single matrix of observations \code{y} then
    \code{probs} is a matrix whose rows correspond to the states
    of the hidden Markov chain, and whose columns correspond to
    the observation times.  If the observations consist of a
    list of observation vectors, then \code{probs} is a list
    of such matrices, one for each matrix of observations.

    Likewise for the \code{means} component of the list returned
    when the argument \code{means} is \code{TRUE}.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{\code{\link{hmm}()}, \code{\link{mps}()},
         \code{\link{viterbi}()}, \code{\link{pr}()},
         \code{\link{fitted.hmm.discnp}()}}
\examples{
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
y.num   <- rhmm(ylengths=rep(300,20),nsim=1,tpm=P,Rho=R,drop=TRUE)
fit.num <- hmm(y.num,K=2,verb=TRUE,keep.y=TRUE,itmax=10)
cpe1    <- sp(model=fit.num,means=TRUE) # Using the estimated parameters.
cpe2    <- sp(y.num,tpm=P,Rho=R,means=TRUE,
              warn=FALSE) # Using the ``true'' parameters.
# The foregoing would issue a warning that Rho had no row names
# were it not for the fact that "warn" has been set to FALSE.
}

\keyword{models}
