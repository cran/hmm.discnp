\name{fitted.hmm.discnp}
\Rdversion{1.1}
\alias{fitted.hmm.discnp}
\title{
  Fitted values of a discrete non-parametric hidden Markov
  model.
}
\description{
  Calculates the fitted values of a discrete non-parametric hidden
  Markov model as the conditional expectations of the observations,
  given the entire observation sequence (and the estimated parameters
  of the model).
}
\usage{
\method{fitted}{hmm.discnp}(object, ...)
}
\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{\dots}{
  Not used.
}
}
\details{
  This function is essentially a wrapper for \code{\link{sp}()}.
  This latter function gets called with \code{mean=TRUE}.
  The fitted values are the conditional expected values of the
  observations, given the entire observation sequence.  Note that the
  \dQuote{conditional expected value} concept makes sense \emph{only}
  if the observations are \emph{numeric}. If this is not the case,
  an error is thrown.  The observation sequence must be present in
  \code{object} (which will be the case if \code{object} was returned
  by \code{hmm()}).  If it is not present an error is thrown.

  However, if such an error is thrown, do not despair!  You
  \emph{do not} have to start from scratch when fitting your
  model with \code{keep.y==TRUE}.  If \code{fit} is your fitted
  model that you obtained \emph{without} setting \code{keep.y==TRUE},
  then you can just re-fit the model using \code{fit} as the
  starting values: \preformatted{
      fit2 <- hmm(<whatever>,par0=fit,keep.y=TRUE)
  }
  This will of course converge instantaneously.  You could also do:
  \preformatted{
      fit2 <- update(fit,data=<whatever>,keep.y=TRUE)
  }
}
\value{
  If the observations (stored in \code{object}) consist of a single
  sequence, then the object returned is a matrix.  If the model
  is univariate (see \code{\link{hmm}()}) then matrix has a single
  column constituting the corresponding sequence of fitted values.
  If the model is bivariate (either independent or dependent) then
  the martrix has two columns corresponding respectively to the two
  variables in the bivariate model.  If the observations consist
  of a list of observation sequences, then the object returned is a
  list of corresponding matrices whose columns are the fitted values.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}
\seealso{ \code{\link{sp}()} \code{link{predict.hmm.discnp}()} }
\examples{
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
lll <- sample(250:350,20,TRUE)
y   <- rhmm(ylengths=lll,nsim=1,drop=TRUE,tpm=P,Rho=R)
fit <- hmm(y,K=2,verb=TRUE,keep.y=TRUE,itmax=10)
fv  <- fitted(fit)
}

\keyword{models}
