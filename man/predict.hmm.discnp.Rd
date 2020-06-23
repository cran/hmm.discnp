\name{predict.hmm.discnp}
\alias{predict.hmm.discnp}
\title{
  Predicted values of a discrete non-parametric hidden Markov
  model.
}
\description{
  Calculates predicted values given a specification of a discrete
  non-parametric hidden Markov model.  The specification may be
  provided in the form of a \code{hmm.discnp} object as returned
  by \code{hmm()} or in the form of \dQuote{components} of such
  a model:  the data \code{y}, the transition probability matrix
  \code{tpm}, the emission probabilities \code{Rho}, etc.  If the
  data are numeric then these predicted values are the conditional
  expectations of the observations, given the entire observation
  sequence (and the --- possibly estimated --- parameters of the
  model).  If the data are categorical (whence \dQuote{expectations}
  make no sense) the \dQuote{predicted values} are taken to be the
  probabilities of each of the possible values of the observations,
  at each time point.
}
\usage{
\method{predict}{hmm.discnp}(object, y = NULL, tpm=NULL, Rho=NULL,
                             ispd=NULL, X=NULL,addIntercept=NULL,
                             warn=TRUE, drop=TRUE, ...)

}
\arguments{
  \item{object}{
  If not \code{NULL}, an object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{y}{A data structure to which the fitted model \code{object}
  \emph{could} have been fitted.  If \code{y} is NULL, an attempt
  is made to extract \code{y} from \code{model}.
}
  \item{tpm,Rho,ispd,X,addIntercept,warn}{
  See the help for \code{\link{sp}()}.
}
  \item{drop}{Logical scalar.  See the help for
  \code{\link{fitted.hmm.discnp}()}.
}
  \item{\dots}{
  Not used.
}
}
\details{
  This function is essentially the same as
  \code{\link{fitted.hmm.discnp}()}.  The main difference is
  that it allows the calculation of fitted/predicted values for a
  data object \code{y} possibly different from that to which the
  model was fitted.  Note that if both the argument \code{y} and
  \code{object[["y"]]} are present, the \dQuote{argument} value takes
  precedence.  This function also allows the model to be specfied
  in terms of individual components rather than as a fitted model
  of class \code{"hmm.discnp"}.  These components, (\code{tpm},
  \code{Rho}, \code{ispd}, \code{X}, \code{addIntercept}) if
  supplied, \emph{take precedence} over the corresponding components
  of \code{object}.  The opposite applies with \code{sp()}.  The
  function \code{fitted.hmm.discnp()} makes use \emph{only}
  of the components of \code{object}.
}
\value{
  See the help for \code{\link{fitted.hmm.discnp}()}.
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
ll1  <- sample(250:350,20,TRUE)
y1   <- rhmm(ylengths=ll1,nsim=1,tpm=P,Rho=R,drop=TRUE)
fit  <- hmm(y1,K=2,verb=TRUE,keep.y=TRUE,itmax=10)
fv   <- fitted(fit)
set.seed(176)
ll2  <- sample(250:350,20,TRUE)
y2   <- rhmm(ylengths=ll2,nsim=1,tpm=P,Rho=R,drop=TRUE)
pv   <- predict(fit,y=y2)
yval <- letters[1:5]
set.seed(171)
y3   <- rhmm(ylengths=ll2,yval=yval,nsim=1,tpm=P,Rho=R,drop=TRUE)
fit3 <- hmm(y3,K=2,verb=TRUE,keep.y=TRUE,itmax=10)
pv3  <- predict(fit3) # Same as fitted(fit3).
}

\keyword{models}
