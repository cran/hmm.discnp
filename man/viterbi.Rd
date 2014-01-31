\name{viterbi}
\alias{viterbi}
\title{ Most probable state sequence. }
\description{
  Calculates ``the'' most probable state sequence underlying
  each of one or more replicate observation sequences.
}
\usage{
viterbi(y, object = NULL, tpm, Rho, ispd=NULL,log=FALSE)
}
\arguments{
  \item{y}{The observations for which the underlying most probable
  hidden states are required.  May be a sequence of observations,
  or a list each entry of which constitutes an independent sequence
  of observations.  If \code{y} is missing (and if \code{object}
  is not \code{NULL}) then \code{y} is extracted from \code{object},
  provided that the \code{y} component of \code{object} is present.
  Otherwise an error is given.}
  \item{object}{ An object describing a hidden Markov model, as
  fitted to the data set \code{y} by \code{\link{hmm}()}.}
  \item{tpm}{ The transition probability matrix for a hidden
  Markov model; ignored if \code{object} is non-null.}
  \item{Rho}{ A matrix specifying the probability distributions
  of the observations for a hidden Markov model; ignored if
  \code{object} is non-null. Should bear some reasonable relationship
  to \code{y}.  If \code{Rho} has row names then these must
  include all values of the observations.  If \code{Rho} does
  not have row names then the values of the observations must
  be integers between \code{1} and \code{nrow(Rho)}.}
  \item{ispd}{ The initial state probability distribution for a hidden
  Markov model; ignored if \code{object} is non-null. Should
  bear some reasonable relationship to \code{y}. If \code{object}
  and \code{ispd} are both \code{NULL} then \code{ispd} is set
  equal to the stationary distribution calculated from \code{tpm}.}
  \item{log}{Logical scalar.  Should logarithms be used in the
  recursive calculations of the probabilities involved in the
  Viterbi algorithm, so as to avoid underflow?  If \code{log} is
  \code{FALSE} then underflow is avoided instead by a normalization
  procedure.  The quantity \code{delta} (see Rabiner 1989, page 264)
  is replaced by \code{delta/sum(delta)} at each step.  It should
  actually make no difference whether \code{log} is set to \code{TRUE}.
  I just included the option because I could.  Also the \code{HMM}
  package uses the logarithm approach so setting \code{log=TRUE}
  might be of interest if comparisons are to be made between the
  results of the two packages.}
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

  If \code{y} consists of a list of observation sequences, the
  value is a list each entry of which is of the form described
  above.
}

\references{
Rabiner, L. R., "A tutorial on hidden Markov models and
selected applications in speech recognition," Proc. IEEE vol. 77,
pp. 257 -- 286, 1989.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
  \url{http://www.stat.auckland.ac.nz/~rolf}
}

\section{Warning}{
There \emph{may} be more than one equally most probable
state sequence underlying a given observation sequence.  This
phenomenon appears to be unlikely to occur in practice.
}

\section{Thanks}{
The correction made to the code so as to avoid underflow
problems was made due to an inquiry and suggestion from
Owen Marshall.
}

\seealso{\code{\link{hmm}()}, \code{\link{sim.hmm}()},
         \code{\link{mps}()}, \code{\link{pr}()},
         \code{\link[HMM]{viterbi}()}}
\examples{
# See the help for sim.hmm() for how to generate y.num and y.let.
\dontrun{
fit.num <- hmm(y.num,K=2,verb=TRUE)
v.1 <- viterbi(object=fit.num)
v.2 <- viterbi(y.num,tpm=P,Rho=R) # P and R as in the
                                  # help for sim.hmm().
# The order of the states has gotten swapped; 3-v.1[[1]] is much
# more similar to v.2[[1]] than is v.1[[1]].

fit.let <- hmm(y.let,K=2,verb=TRUE)
v.3 <- viterbi(object=fit.let)    # Works.
v.4 <- viterbi(y.let,tpm=P,Rho=R) # Throws an error (R has no row names.)
}
}
\keyword{models}
