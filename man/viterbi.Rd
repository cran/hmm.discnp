\name{viterbi}
\alias{viterbi}
\title{ Most probable state sequence. }
\description{
  Calculates ``the'' most probable state sequence underlying
  each of one or more replicate observation sequences.
}
\usage{
viterbi(y, model = NULL, tpm, Rho, ispd=NULL,log=FALSE, warn=TRUE)
}
\arguments{
  \item{y}{The observations for which the most probable sequence(s)
  of underlying hidden states are required.  May be a sequence of
  observations in the form of a vector or a one or two column matrix,
  or a list each component of which constitutes a (replicate)
  sequence of observations.  It may also be an object of class
  \code{"multipleHmmDataSets"} as returned by \code{\link{rhmm}()}
  with \code{nsim>1}.

  If \code{y} is missing, it is extracted from
  \code{model} (whence it will \emph{not} be of class
  \code{"multipleHmmDataSets"}!) provided that \code{model} and
  its \code{y} component are not \code{NULL}.  Otherwise an error
  is given.
}
  \item{model}{ An object describing a hidden Markov model, as
  fitted to the data set \code{y} by \code{\link{hmm}()}.}
  \item{tpm}{ The transition probability matrix for a hidden
  Markov model; ignored if \code{model} is non-null.}
  \item{Rho}{An object specifying the probability distributions
  of the observations for a hidden Markov model.  See
  \code{\link{hmm}()}. Ignored if \code{model} is non-null. Should
  bear some reasonable relationship to \code{y}.

  If \code{Rho} has dimension names (or if its entries have
  dimension names in the case where \code{Rho} is a list) then the
  appropriate dimension names must include all corresponding values
  of the observations.  If a relevant vector of dimension names is
  \code{NULL} then it is formed as the sort unique values of the
  approprate columns of the observation matrices.   In this case the
  corresponding dimensions must match the number of unique values.}

  \item{ispd}{ The initial state probability distribution for a hidden
  Markov model; ignored if \code{model} is non-null. Should
  bear some reasonable relationship to \code{y}. If \code{model}
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
  might be of interest if comparisons are to be made between
  results from the two packages.
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

  If \code{y} is of class \code{"multipleHmmDataSets"} then the
  value returned is a list of lists of the sort described above.
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
There \emph{may} be more than one equally most probable state
sequence underlying a given observation sequence.  This phenomenon
can occur but appears to be unlikely to do so in practice.
}

\section{Thanks}{
The correction made to the code so as to avoid underflow
problems was made due to an inquiry and suggestion from
Owen Marshall.
}

\seealso{\code{\link{hmm}()}, \code{\link{rhmm}()},
         \code{\link{mps}()}, \code{\link{pr}()}}
\examples{
# See the help for logLikHmm() for how to generate y.num and y.let.
\dontrun{
fit.num     <- hmm(y.num,K=2,verb=TRUE,keep.y=TRUE)
v.1         <- viterbi(model=fit.num)
rownames(R) <- 1:5 # Avoids a (harmless) warning.
v.2         <- viterbi(y.num,tpm=P,Rho=R)
# P and R as in the help for logLikHmm() and for sp().

# Note that the order of the states has gotten swapped; 3-v.1[[1]]
# is identical to v.2[[1]]; for other k = 2, ..., 20, 3-v.1[[k]]
# is much more similar to v.2[[k]] than is v.1[[k]].

fit.let <- hmm(y.let,K=2,verb=TRUE,keep.y=TRUE))
v.3 <- viterbi(model=fit.let)
rownames(R) <- letters[1:5]
v.4 <- viterbi(y.let,tpm=P,Rho=R)
}
}
\keyword{models}
