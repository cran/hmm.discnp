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
sim.hmm(nsim, tpm, Rho, nrep=1)
}
\arguments{
\item{nsim}{
The length of each sequence of observations.
}
\item{tpm}{
The transition probability matrix for the underlying hidden Markov
chain(s).  Note that the rows of \code{tpm} must sum to 1.
}
\item{Rho}{
A matrix specifying the probability of an observation taking on one
of a set of possible values, given the state of the underlying hidden
Markov chain.  Note that the columns of \code{Rho} must sum to 1.
}
\item{nrep}{
The number of independent replicate sequences of observations.
}}
\value{
A matrix with \code{nsim} rows and \code{nrep} columns; if
\code{nrep}==1, then the returned value is a vector of length
\code{nsim}.  Each column of the matrix forms an independent
replicate of observations from a hidden Markov model.  The values of
those observations are integers between 1 and \code{nrow(Rho)}.
}
\seealso{\code{\link{hmm}()}}
\examples{
P     <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R     <- matrix(c(0.33,0.24,0.05,0.24,0.23,
                  0.23,0.33,0.25,0.06,0.04),
                  5,2,byrow=TRUE)
y.sim <- sim.hmm(300,P,R,20)
}
\keyword{models}
