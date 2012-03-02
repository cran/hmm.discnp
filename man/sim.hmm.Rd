\name{sim.hmm}
\alias{sim.hmm}
\title{
Simulate discrete data from a hidden Markov model.
}
\description{
Simulates one or more replicates of discrete data
from a model such as is fitted by the function \code{hmm()}.
}
\usage{
sim.hmm(nsim,tpm,Rho,ispd=NULL,yval=NULL,verb=FALSE)
}
\arguments{
\item{nsim}{
Vector of the lengths of the sequences of observations
to be generated.
}
\item{tpm}{
The transition probability matrix for the underlying hidden Markov
chain(s).  Note that the rows of \code{tpm} must sum to 1.
Ignored if \code{ncol(Rho)==1}.
}
\item{Rho}{
A matrix specifying the probability of an observation taking on
one of a set of possible values, given the state of the underlying
hidden Markov chain.  Note that the columns of \code{Rho} must sum
to 1.  If \code{ncol(Rho)==1} the simulated data are i.i.d. from
the distribution specified by the single column of \code{Rho}.
}
\item{ispd}{A vector specifying the initial state probability
distribution of the chain.  If this is not specified it is
taken to be the stationary distribution of the chain, calculated
from \code{tpm}.
}
\item{yval}{Vector of numbers or character strings constituting
the possible values of the observations.  If not supplied it
is taken to equal \code{1:nrow(Rho)}.}
\item{verb}{Logical scalar.  If \code{TRUE} then the overall index
of the simulated value that has been reached is printed out every
1000 iterations.  Useful for reassurance when very \dQuote{large}
simulations are undertaken.}
}
\value{
If \code{length(nsim)==1} then the value returned is a vector of
length \code{nsim}.  If \code{length(nsim)>1} then the value returned
is a list of the same length as \code{nsim}, each component of which
is an independent vector of simulated observations.  The length
of component \code{i} of this list is equal to \code{nsim[i]}.
The values of the observations are entries of \code{yval}.
}

\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
  \url{http://www.math.unb.ca/~rolf}
}

\seealso{\code{\link{hmm}()}}
\examples{
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
lll   <- sample(250:350,20,TRUE)
y.num <- sim.hmm(lll,P,R)
y.let <- sim.hmm(lll,P,R,yval=letters[1:5])
}
\keyword{models}
