\name{pr}
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
pr(s, y, model=NULL, tpm, Rho, ispd=NULL, warn=TRUE)
}
\arguments{
  \item{s}{
  A sequence of states of the underlying Markov chain, or
  a list of such sequences.
}
  \item{y}{
  A one or two column matrix of observations from a hidden Markov
  model, corresponding to the state sequence \code{s}, or a list
  of such matrices corresponding to the state sequences in the
  list \code{s}.  If \code{y} consists of a single matrix,
  it is assumed to correspond to each of the state sequences in
  \code{s} in turn.  Otherwise the length of the list \code{y}
  must be the same as the length of the list \code{s}. (If not,
  then an error is given).  If \code{y} is missing, it is extracted
  from \code{model} provided that \code{model} and its \code{y}
  component are not \code{NULL}.  Otherwise an error is given.
}
  \item{model}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{tpm}{
  The transition probability matrix of the chain.  Ignored (and
  extracted from \code{model} instead) if \code{model} is not
  \code{NULL}.
}
  \item{Rho}{
  An object specifying the \dQuote{emission} probabilities of
  observations, given the underlying state.  See \code{\link{hmm}()}.
  Ignored (and extracted from \code{model} instead) if \code{model}
  is not \code{NULL}.
}
  \item{ispd}{
  The vector specifying the initial state probability distribution
  of the Markov chain.  Ignored (and extracted from \code{model}
  instead) if \code{model} is not \code{NULL}.  If both \code{ispd}
  and \code{model} are NULL then \code{ispd} is taken to be the
  stationary distribution of the chain, calculated from \code{tpm}.
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
# See the help for rhmm() for how to generate y.num.
\dontrun{
fit.num <- hmm(y.num,K=2,verb=TRUE)
# Using fitted parmeters.
s.vit.1   <- viterbi(y.num,fit.num)
pr.vit.1  <- pr(s.vit.1,model=fit.num)
# Using true parameters from which y.num was generated.
s.vit.2   <- viterbi(y.num,tpm=P,Rho=R)
pr.vit.2  <- pr(s.vit.2,y.num,tpm=P,Rho=R)
}
}
\keyword{models}
