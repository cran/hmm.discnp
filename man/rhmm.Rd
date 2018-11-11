\name{rhmm}
\alias{rhmm}
\alias{rhmm.default}
\alias{rhmm.hmm.discnp}
\title{
Simulate discrete data from a non-parametric hidden Markov model.
}
\description{
Simulates one or more replicates of discrete data
from a model such as is fitted by the function \code{hmm()}.
}
\usage{
rhmm(model,...,nsim,verbose=FALSE)
\method{rhmm}{default}(model, ..., nsim=1, verbose=FALSE, ylengths,
                       nafrac=NULL, fep=NULL, tpm, Rho, ispd=NULL, yval=NULL,
                       drop=TRUE)
\method{rhmm}{hmm.discnp}(model, ..., nsim=1, verbose=FALSE, inMiss=TRUE,
                          fep=NULL, drop=TRUE)
}
\arguments{
\item{model}{
   An object of class \code{hmm.discnp}.  This will have the form
   of a list specifying a hidden Markov model with discrete emissions
   and emission probabilities specified non-parametrically, i.e.
   by means of a table or tables.  Usually this will be an object
   returned by \code{\link{hmm}()}.  This argument is ignored
   by the default method.
}

\item{...}{
   Not used.
}

\item{nsim}{
    Integer scalar; the number of data sets to be simulated.
}

\item{verbose}{Logical scalar.  If \code{TRUE} then the overall index
of the simulated value that has been reached is printed out every
1000 iterations.  Useful for reassurance when very \dQuote{large}
simulations are undertaken.
}

\item{ylengths}{
    Integer values vector specify the lengths (or number of rows
    in the bivariate setting) of the individual observation
    sequences constituting a data set.
}
\item{nafrac}{
    See \code{\link{misstify}()} for an explanation of this argument.
    If specified a fraction \code{nafrac[[j]]} of column \code{j}
    of the data will be randomly set equal to \code{NA}.
}

\item{fep}{
    \dQuote{First entry present}.  See \code{\link{misstify}()}
    for an explanation of this argument.
}

\item{tpm}{
The transition probability matrix for the underlying hidden
Markov chain(s).  Note that the rows of \code{tpm} must sum to 1.
Ignored if \code{ncol(Rho)==1}.  Ignored by the \code{hmm.discnp}
method and extracted from \code{model}.
}

\item{Rho}{
An object specifying the probability distribution of the
observations, given the state of the underlying hidden Markov chain.
(I.e. the \dQuote{emission} probabilities.)  See \code{\link{hmm}()}.
Note that \code{Rho} can be such that the number of states is 1, in
which case the simulated data are i.i.d. from the single distribution
specified by \code{Rho}.  Ignored by the \code{hmm.discnp} method
and extracted from \code{model}.
}

\item{ispd}{A vector specifying the initial state probability
distribution of the chain.  If this is not specified it is taken
to be the stationary distribution of the chain, calculated from
\code{tpm}.  Ignored by the \code{hmm.discnp} method and extracted
from \code{model}.
}

\item{yval}{Vector of possible values of the observations, or (in
the bivariate setting) a list of two such vectors.  If not supplied
it is formed from the appropriates dimension names of \code{Rho}
or of the entries of \code{Rho} if \code{Rho} is a list.  Ignored by
the \code{hmm.discnp} method.
}
\item{drop}{
    Logical scalar; if \code{nsim} is 1 and if \code{drop} is
    \code{TRUE} then the list to be returned by this function
    is replaced by its first and only entry. (See below.)
}
\item{inMiss}{
    Logical scalar; if \code{TRUE} then missing values will
    be randomly inserted into the data in the fraction
    \code{nafrac} determined from \code{object}.
}
}

\value{
If \code{nsim>1} or \code{drop} is \code{FALSE} then the value
returned is a list of length \code{nsim}.  Each entry of this
list is in turn a list of the same length as \code{ylengths},
each component of which is an independent  vector or matrix of
simulated observations.  The length or number of rows of component
\code{i} of this list is equal to \code{ylengths[i]}.  The values
of the observations are entries of \code{yval} or of \emph{its}
entries when \code{yval} is a list.

If \code{nsim=1} and \code{drop} is \code{TRUE} then the (\dQuote{outer})
list described above is replaced by its first and only entry
}

\section{Note}{
You may find it useful to avail yourself of the function
\code{\link{nafracCalc}()} to determine the fraction of missing
values in a given existing (presumably \dQuote{real}) data set.
}

\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{\code{\link{hmm}()} \code{\link{nafracCalc}()}
  \code{\link{misstify}()}
}

\examples{
# To do: one or more bivariate examples.
\dontrun{
    y <- list(linLandFlows$deciles,ftLiardFlows$deciles)
    fit <- hmm(y,K=3)
    simX <- rhmm(fit)
}
}
\keyword{datagen}
