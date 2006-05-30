\name{viterbi}
\alias{viterbi}
\title{ Most probable state sequence. }
\description{
  Calculates ``the'' most probable state sequence underlying
  each of one or more replicate observation sequences.
}
\usage{
viterbi(y, object = NULL, tpm, Rho, ispd, yval = NULL)
}
\arguments{
  \item{y}{The observations for which the underlying most
  probable hidden states are required.  May be a sequence
  of observations, or a matrix each column of which constitutes
  a (replicate) sequence of observations.}
  \item{object}{ An object describing a hidden Markov model, as
  fitted to the data set \code{y} by \code{\link{hmm}()}.}
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
  Applies the Viterbi algorithm to calculate ``the'' most probable
  robable state sequence underlying each observation sequences.
}
\value{
  If \code{y} consists of a single observation sequence, the
  value is the underlying most probable observation sequence,
  or a matrix whose columns consist of such sequences if there
  is more than one (equally) most probable sequence.

  If \code{y} consists of a matrix of observation sequences, the
  value is a list each entry of which is of the form described
  above.
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
There \emph{may} be more than one equally most probable
state sequence underlying a given observation sequence.  This
phenomenon appears to be unlikely to occur in practice.
}

\seealso{\code{\link{hmm}()}, \code{\link{sim.hmm}()},
         \code{\link{mps}()}}
\examples{
# See the help for sim.hmm() for how to generate y.sim.
\dontrun{
try <- hmm(y.sim,K=2,verb=TRUE)
vvv.1 <- viterbi(y.sim,try)
vvv.2 <- viterbi(y.sim,tpm=P,ispd=c(0.25,0.75),Rho=R) # P and R as in the
                                                      # help for sim.hmm().
# The order of the states has gotten swapped; 3-vvv.1[[1]] is much
# more similar to vvv.2[[1]] than is vvv.1[[1]].

}
}
\keyword{models}
