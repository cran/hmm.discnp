\name{hmm}
\alias{hmm}
\title{
Fit a hidden Markov model to discrete data.
}
\description{
Uses the EM algorithm to perform a maximum likelihood fit of a hidden
Markov model to discrete data where the observations come from one of
a number of finite discrete distributions, depending on the (hidden)
state of the Markov chain.  These distributions are specified
(non-parametrically) by a matrix \eqn{R = [\rho_{ij}]}{Rho =
[rho_ij]} where \eqn{\rho_{ij} = P(Y = y_i | S = j)}{rho_ij = P(Y =
y_i | S = j)}, \eqn{Y} being the observable random variable and
\eqn{S} being the hidden state.
}
\usage{
hmm(y, yval=NULL, par0=NULL, K=NULL, rand.start=NULL, mixture=FALSE,
    tolerance=1e-4, verbose=FALSE, itmax=200, crit='PCLL', data.name=NULL)
}
\arguments{
\item{y}{
A vector or matrix of discrete data; missing values are allowed.  If
\code{y} is a matrix, each column is interpreted as an independent
replicate of the observation sequence.
}
\item{yval}{
A vector of possible values for the data; it defaults to the sorted
unique values of \code{y}.  If any value of \code{y} does not match
some value of \code{yval}, it will be treated as a MISSING VALUE.
}
\item{par0}{
An optional (\emph{named}) list of starting values for the parameters
of the model, with components \code{tpm} (transition probability
matrix) and \code{Rho}.  The matrix \code{Rho} specifies the
probability that the observations take on each value in yval, given
the state of the hidden Markov chain.  The columns of \code{Rho}
correspond to states, the rows to the values of \code{yval}.

If \code{par0} is not specified, starting values are created by the
function \code{init.all()}.
}
\item{K}{
The number of states in the hidden Markov chain; if \code{par0} is
not specified \code{K} MUST be; if \code{par0} is specified, \code{K}
is ignored.

Note that \code{K=1} is acceptable; if \code{K} is 1 then all
observations are treated as being independent and the non-parametric
estimate of the distribution of the observations is calculated
in the obvious way.
}
\item{rand.start}{
A list consisting of two logical scalars which must be named
\code{tpm} and \code{Rho}, if \code{tpm} is TRUE then the function
init.all() chooses entries for then starting value of \code{tpm} at
random; likewise for \code{Rho}.  This argument defaults to
\code{list(tpm=FALSE,Rho=FALSE)}.
}
\item{mixture}{
A logical scalar; if TRUE then a mixture model (all rows of the
transition probability matrix are identical) is fitted rather than a
general hidden Markov model.
}
\item{tolerance}{
If the value of the quantity used for the stopping criterion is
less than tolerance then the EM algorithm is considered to have
converged.
}
\item{verbose}{
A logical scalar determining whether to print out details of the
progress of the EM algorithm.
}
\item{itmax}{
If the convergence criterion has not been met by the time \code{itmax}
EM steps have been performed, a warning message is printed out,
and the function stops.  A value is returned by the function
anyway, with the logical component "converged" set to FALSE.
}
\item{crit}{
The name of the stopping criterion, which must be one of "PCLL"
(percent change in log-likelihood; the default), "L2" (L-2 norm,
i.e.  square root of sum of squares of change in coefficients), or
"Linf" (L-infinity norm, i.e.  maximum absolute value of change in
coefficients).
}
\item{data.name}{
An identifying tag for the fit; if omitted, it defaults to the
name of data set \code{y} as determined by \code{deparse(substitute(y))}.
}
}
\value{
A list with components:

\item{Rho}{
The fitted value of the probability matrix \code{Rho} specifying the
distributions of the observations.
}
\item{tpm}{
The fitted value of the transition probabilty matrix \code{tpm}.
}
\item{ispd}{
The fitted initial state probability distribution, assumed to
be the (unique) stationary distribution for the chain, and thereby
determined by the transition probability matrix \code{tpm}.
}
\item{log.like}{
The final value of the log likelihood, as calculated through
recursion.
}
\item{converged}{
A logical scalar saying whether the algorithm satisfied the
convergence criterion before the maximum of itmax EM steps was
exceeded.
}
\item{nstep}{
The number of EM steps performed by the algorithm.
}
\item{data.name}{
An identifying tag, specified as an argument, or determined from
the name of the argument y by deparse(substitute(y)).
}}

\section{Note}{
If \code{K=1} then \code{tpm}, \code{ispd}, \code{converged}, and
\code{nstep} are all set equal to \code{NA} in the list returned
by this function.
}

\details{
The hard work is done by a Fortran subroutine "recurse" (actually
coded in Ratfor) which is dynamically loaded.
}

\references{
Rabiner, L. R., "A tutorial on hidden Markov models and
selected applications in speech recognition," Proc. IEEE vol. 77,
pp. 257 -- 286, 1989.

Zucchini, W. and Guttorp, P., "A hidden Markov model
for space-time precipitation," Water Resources Research vol.
27, pp. 1917-1923, 1991.

MacDonald, I. L., and Zucchini, W., "Hidden Markov
and Other Models for Discrete-valued Time Series, Chapman \& Hall,
London, 1997.

Liu, Limin, "Hidden Markov Models for Precipitation in a Region
of Atlantic Canada", Master's Report, University of New Brunswick,
1997.
}

\section{Warning}{
The ordering of the (hidden) states can be arbitrary.  What
the estimation procedure decides to call ``state 1'' may
not be what \emph{you} think of as being state number 1.
The ordering of the states will be affected by the starting
values used.
}

\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}
  \url{http://www.math.unb.ca/~rolf}
}

\seealso{\code{\link{sim.hmm}()}, \code{\link{mps}()},
         \code{\link{viterbi}()}}
\examples{
# See the help for sim.hmm() for how to generate y.sim.
\dontrun{
try <- hmm(y.sim,K=2,verb=TRUE)
try.mix <- hmm(y.sim,K=2,verb=TRUE,mixture=TRUE)
}
}
\keyword{models}
