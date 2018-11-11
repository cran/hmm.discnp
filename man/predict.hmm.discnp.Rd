\name{predict.hmm.discnp}
\Rdversion{1.1}
\alias{predict.hmm.discnp}
\title{
  Predicted values of a discrete non-parametric hidden Markov
  model.
}
\description{
  Calculates the predicted values of a discrete non-parametric hidden
  Markov model as the conditional expectations of the observations,
  given the entire observation sequence (and the [estimated] parameters
  of the model).
}
\usage{
\method{predict}{hmm.discnp}(object, y=NULL, ...)
}
\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{y}{A data structure from which the fitted model \code{object}
  \emph{could} have been fitted.
}
  \item{\dots}{
  Not used.
}
}
\details{
  This function is essentially the same as
  \code{\link{fitted.hmm.discnp}()}, except that it allows
  the calculation of fitted/predicted values for a data object
  possibly different from that to which the model was fitted.
  Like \code{fitted.hmm.discnp()} this function is essentially a
  wrapper for \code{\link{sp}()}.  Note that the predicted values
  are calculated as \dQuote{conditional expected values} which makes
  sense \emph{only} if the observations are \emph{numeric}. If this
  is not the case, an error is thrown.  If the observation list
  \code{y} is not supplied then it is taken to be the \code{y}
  component of \code{object}.  If there is no such component
  (i.e. if \code{object} is the value returned from \code{hmm()}
  with the argument \code{keep.y} set equal to \code{FALSE}) then
  an error is thrown.
}
\value{
  If the observations (stored in \code{y}) consist of a single
  sequence, then the object returned is a matrix.  If the model
  is univariate (see \code{\link{hmm}()}) then matrix has a single
  column constituting the corresponding sequence of fitted values.
  If the model is bivariate (either independent or dependent) then
  the matrix has two columns corresponding respectively to the two
  variables in the bivariate model.  If the observations consist
  of a list of observation sequences, then the object returned is a
  list of corresponding matrices whose columns are the fitted values.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}
\seealso{ \code{\link{sp}()} \code{link{fitted.hmm.discnp}()} }
\examples{
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
ll1 <- sample(250:350,20,TRUE)
y1  <- rhmm(ylengths=ll1,nsim=1,tpm=P,Rho=R,drop=TRUE)
fit <- hmm(y1,K=2,verb=TRUE,keep.y=TRUE,itmax=10)
fv  <- fitted(fit)
set.seed(176)
ll2 <- sample(250:350,20,TRUE)
y2  <- rhmm(ylengths=ll2,nsim=1,tpm=P,Rho=R,drop=TRUE)
pv  <- predict(fit,y=y2)
}

\keyword{models}
