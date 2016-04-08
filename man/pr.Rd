\name{pr}
\Rdversion{1.1}
\alias{pr}
\title{
  Probability of state sequences.
}
\description{
  Calculates the conditional probability of one or more state sequences,
  given the corresponding observations sequences (and the model
  parameters.
}
\usage{
pr(s,y,object=NULL,tpm,Rho,ispd=NULL)
}
\arguments{
  \item{s}{
  A sequence of states of the underlying Markov chain, or
  a list of such sequences.
}
  \item{y}{
  A sequence of observations from a hidden Markov model, corresponding
  to the state sequence \code{s}, or a list of such sequences
  corresponding to the state sequences in the list \code{s}.  If
  \code{y} is or consists of a single sequence, it is assumed to
  correspond to each of the state sequences in \code{s} in turn.
  Otherwise the length of the list \code{y} must be the same as the
  length of the list \code{s} (otherwise an error is given).  If
  \code{y} is missing, it is extracted from \code{object} provided
  that \code{object} and its \code{y} component are not \code{NULL}.
  Otherwise an error is given.
}
  \item{object}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{tpm}{
  The transition probability matrix of the chain.  Ignored (and
  extracted from \code{object} instead) if \code{object} is not
  \code{NULL}.
}
  \item{Rho}{
  The matrix of probabilities specifying the distribution of the
  observations, given the underlying state.  The rows of this matrix
  correspond to the possible values of the observations, the columns
  to the states.  Ignored (and extracted from \code{object} instead)
  if \code{object} is not \code{NULL}.
}
  \item{ispd}{
  The vector specifying the initial state probability distribution
  of the Markov chain.  Ignored (and extracted from \code{object}
  instead) if \code{object} is not \code{NULL}.  If both \code{ispd}
  and \code{object} are NULL then \code{ispd} is taken to be the
  stationary distribution of the chain, calculated from \code{tpm}.
}
}
\value{
  The probability of \code{s} given \code{y}, or a vector
  of such probabilities if \code{s} and \code{y} are lists.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\section{Warning}{
  The conditional probabilities will be tiny if the sequences
  involved are of any substantial length.  Underflow may be a
  problem.  The implementation of the calculations is not
  sophisticated.
}

\seealso{\code{\link{hmm}()}, \code{\link{mps}()},
         \code{\link{viterbi}()}, \code{\link{sp}()},
         \code{\link{fitted.hmm.discnp}()}}

\examples{
# See the help for sim.hmm() for how to generate y.num.
\dontrun{
fit.num <- hmm(y.num,K=2,verb=TRUE)
# Using fitted parmeters.
s.vit.1   <- viterbi(y.num,fit.num)
pr.vit.1  <- pr(s.vit.1,object=fit.num)
# Using true parameters from which y.num was generated.
s.vit.2   <- viterbi(y.num,tpm=P,Rho=R)
pr.vit.2  <- pr(s.vit.2,y.num,tpm=P,Rho=R)
}}
\keyword{models}
