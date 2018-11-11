\name{mps}
\alias{mps}
\title{ Most probable states. }
\description{
  Calculates the most probable hidden state underlying each
  observation.
}
\usage{
mps(y, model = NULL, tpm, Rho, ispd=NULL, warn=TRUE)
}
\arguments{
  \item{y}{The observations for which the underlying most
  probable hidden states are required.  May be a sequence
  of observations, or a list each component of which constitutes
  a (replicate) sequence of observations.  If \code{y} is
  missing it is set equal to the \code{y} component of
  \code{model}, given that that object and that component
  exist.  Otherwise an error is given.}
  \item{model}{ An object describing a fitted hidden Markov
   model, as returned by \code{\link{hmm}()}.  In order to
   make any kind of sense, \code{model} should bear some
   reasonable relationship to \code{y}.}
  \item{tpm}{ The transition probability matrix for a hidden
  Markov model; ignored if \code{model} is non-null. Should
  bear some reasonable relationship to \code{y}.}
  \item{Rho}{An object specifying the probability distributions of
  the observations (\dQuote{emission} probabilities) for a hidden
  Markov model.  See \code{\link{hmm}()}.  Ignored if \code{model}
  is non-null. Should bear some reasonable relationship to \code{y}.}
  \item{ispd}{A vector specifying the initial state probability
  distribution for a hidden Markov model, or a matrix each of whose
  columns are trivial (\dQuote{delta function}) vectors specifying
  the \dQuote{most probable} initial state for each observation
  sequence.

  This argument is ignored if \code{model} is non-null. It should
  bear some reasonable relationship to \code{y}.  If both \code{ispd}
  and \code{model} are \code{NULL} then \code{ispd} is taken to
  be the stationary distribution of the chain, calculated from
  \code{tpm}.
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
  For each \eqn{t} the maximum value of \eqn{\gamma_t(i)}{gamma_t(i)},
  i.e. of the (estimated) probability that the state at time \eqn{t}
  is equal to \eqn{i}, is calculated, and the value of the state
  with the  corresponding index is returned.
}
\value{
  If \code{y} is a single observation sequence, then the
  value is a vector of corresponding most probable states.

  If \code{y} is a list of replicate sequences, then the value is
  a list, the \eqn{j}-th entry of which constitutes the vector of
  most probable states underlying the \eqn{j}-th replicate sequence.
}

\references{
Rabiner, L. R., "A tutorial on hidden Markov models and
selected applications in speech recognition," Proc. IEEE vol. 77,
pp. 257 -- 286, 1989.
}

\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\section{Warning}{

The \emph{sequence of most probable states} as calculated by this
function will not in general be the \emph{most probable sequence of
states}.  It may not even be a \emph{possible} sequence of states.
This function looks at the state probabilities separately for each
time \eqn{t}, and not at the states in their sequential context.

To obtain the most probable sequence of states use
\code{\link{viterbi}()}.
}

\seealso{\code{\link{hmm}()}, \code{\link{rhmm}()},
         \code{\link{viterbi}()}}

\examples{
# See the help for rhmm() for how to generate y.num.
\dontrun{
fit.num <- hmm(y.num,K=2,verb=TRUE)
s.1 <- mps(y.num,fit.num)
s.2 <- mps(y.num,tpm=P,ispd=c(0.25,0.75),Rho=R) # P and R as in the help
                                                # for rhmm().
# The order of the states has gotten swapped; (3 - s.1[,1]) is much
# more similar to s.2[,1] than is s.1[,1].
}
}
\keyword{models}
