\name{logLikHmm}
\alias{logLikHmm}
\title{
Log likelihood of a hidden Markov model
}
\description{
Calculate the log likelihood of a hidden Markov model with
discrete non-parametric observation distributions.
}
\usage{
logLikHmm(y, par)
}
\arguments{
  \item{y}{
  A sequence, or list of sequences, of observations from a hidden
  Markov model with discrete non-parametric observation distributions.
}
  \item{par}{
  A list with entries
  \itemize{
  \item \code{Rho}: The matrix whose columns consist of the
  probabilities specifying the distributions of the observations.
  \item \code{tpm}: The transition probability matrix of the Markov chain.
  \item (Optionally) \code{ispd}: The vector of probabilities specifying
  the initial state probability distribution.  If this component of \code{par}
  is not present then \code{ispd} is calculated as the stationary distribution
  determined by \code{tpm}.
  }
  This list may be an object returned by \code{\link{hmm}()}.
}
}
\details{
The observation values (the entries of \code{y}, or of the sequences
which constitute the entries of \code{y} if \code{y} is a list)
must be members of the set of row names of \code{Rho}.  More specifically,
if \code{Rho} has no row names then these observation values must all be
integers between \code{1} and \code{nrow(Rho)}.  If \code{Rho} \emph{has}
row names then these observations must be character strings and must be
members of \code{row.names(Rho)}.
}
\value{
The loglikehood of \code{y} given the parameter values specified in \code{par}.
}
\references{
See \code{\link{hmm}()} for references.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
  \url{http://www.stat.auckland.ac.nz/~rolf}
}

\seealso{
\code{\link{hmm}()}, \code{\link{pr}()}, \code{\link{sp}()}
}
\examples{
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
lll   <- sample(250:350,20,TRUE)
y.num <- sim.hmm(lll,P,R)
y.let <- sim.hmm(lll,P,R,yval=letters[1:5])
ll1 <- logLikHmm(y.num,list(tpm=P,Rho=R))
row.names(R) <- letters[1:5]
ll2 <- logLikHmm(y.let,list(tpm=P,Rho=R))
# Note: ll2 differs from ll1 since y.num and y.let were
# simulated from different random seeds.
ll3 <- logLikHmm(y.let,list(tpm=P,Rho=R,ispd=c(0.5,0.5)))
fit <- hmm(y.num,K=2)
ll4 <- logLikHmm(y.num,fit) # Use the fitted rather than the "true" parameters.
}
\keyword{ models }
