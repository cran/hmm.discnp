\name{mps}
\alias{mps}
\title{ Most probable states. }
\description{
  Calculates the most probable hidden state underlying each
  observation.
}
\usage{
mps(y, object = NULL, tpm, Rho, ispd, yval = NULL)
}
\arguments{
  \item{y}{The observations for which the underlying most
  probable hidden states are required.  May be a sequence
  of observations, or a matrix each column of which constitutes
  a (replicate) sequence of observations.}
  \item{object}{ An object describing a fitted hidden Markov
   model, as returned by \code{\link{hmm}()}.  In order to
   make any kind of sense, \code{object} should bear some
   reasonable relationship to \code{y}.}
  \item{tpm}{ The transition probability matrix for a hidden
  Markov model; ignored if \code{object} is non-null. Should
  bear some reasonable relationship to \code{y}.}
  \item{Rho}{ A matrix specifying the probability distributions
  of the observations for a hidden
  Markov model; ignored if \code{object} is non-null. Should
  bear some reasonable relationship to \code{y}.}
  \item{ispd}{ The initial state probability distribution for a hidden
  Markov model; ignored if \code{object} is non-null. Should
  bear some reasonable relationship to \code{y}.}
  \item{yval}{ The set of unique values of the observations;
  calculated from the observations \code{y} if left \code{NULL}.}
}
\details{
  For each \eqn{t} the maximum value of
  \eqn{\gamma_t(i)}{gamma_t(i)}, i.e. of the (estimated) probability
  that the state at time \eqn{t} is equal to \eqn{i}, is calculated,
  and the corresponding index returned.  These indices are interpreted
  as the values of the (most probable) states.  I.e. the states
  are assumed to be 1, 2, \dots, \eqn{K}, for some \eqn{K}.
}
\value{
  If \code{y} is a single observation sequence, then the
  value is a vector of corresponding most probable states.

  If \code{y} is a matrix of replicate sequences, then the value is
  a matrix, the \eqn{j}-th column of which constitutes the vector of
  most probable states underlying the \eqn{j}-th replicate sequence.
}

\references{
Rabiner, L. R., "A tutorial on hidden Markov models and
selected applications in speech recognition," Proc. IEEE vol. 77,
pp. 257 -- 286, 1989.
}

\author{Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}

\section{Warning}{

The \emph{sequence} of most probable states as calculated by this
function will \emph{not} in general be the most probable sequence of
states.  It may not even be a \emph{possible} sequence of states.
This function looks at the state probabilities separately for each
time \eqn{t}, and not at the states in their sequential context.

To obtain the most probable sequence of states use
\code{\link{viterbi}()}.
}

\seealso{\code{\link{hmm}()}, \code{\link{sim.hmm}()},
         \code{\link{viterbi}()}}

\examples{
# See the help for sim.hmm() for how to generate y.sim.
\dontrun{
try <- hmm(y.sim,K=2,verb=TRUE)
sss.1 <- mps(y.sim,try)
sss.2 <- mps(y.sim,tpm=P,ispd=c(0.25,0.75),Rho=R) # P and R as in the help
                                                  # for sim.hmm().
# The order of the states has gotten swapped; 3-sss.1[,1] is much
# more similar to sss.2[,1] than is sss.1[,1].
}
}
\keyword{models}
