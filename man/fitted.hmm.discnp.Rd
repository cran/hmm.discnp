\name{fitted.hmm.discnp}
\alias{fitted.hmm.discnp}
\title{
  Fitted values of a discrete non-parametric hidden Markov
  model.
}
\description{
  Calculates the fitted values of a discrete non-parametric hidden
  Markov model.  If the data are numeric these are the conditional
  expectations of the observations, given the entire observation
  sequence (and the estimated parameters of the model).  If the data
  are categorical (whence \dQuote{expectations} make no sense) the
  \dQuote{fitted values} are taken to be the probabilities of
  each of the possible values of the observations, at each time
  point.
}
\usage{
\method{fitted}{hmm.discnp}(object, warn=TRUE, drop=TRUE, ...)
}
\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by
  \code{\link{hmm}()}.
}
  \item{warn}{
  Logical scalar.  See the help for \code{\link{sp}()}.
}
  \item{drop}{
  Logical scalar.  If there is a single sequence of observations
  (i.e. if \code{object[["y"]]} consists of a matrix or a list of
  length 1) and if \code{drop} is \code{TRUE} then the returned value
  is a single entity (matrix, list of two matrices, or 3-dimensional
  array, depending on circumstances.  Otherwise the returned value
  is a list of such entities, one for each observation sequence.
  }
  \item{\dots}{
  Not used.
}
}
\details{
  The observation sequence(s) must be present in \code{object} (which
  will be the case if \code{object} was returned by \code{hmm()}
  and if the argument \code{keep.y} was set to \code{TRUE}).
  If it is not present an error is thrown.

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
  If the observations consist of a single sequence and if
  \code{drop} is \code{TRUE} then the returned value consists
  of a single object (matrix, list of two matrices, or 3-dimensional
  array, depending on circumstances; see below).  Otherwise the
  returned value is a list of such objects, one for each observation
  sequence.

  If the observations are numeric then the object corresponding
  to each observation sequence is a matrix.  If the model is
  univariate (see \code{\link{hmm}()}) then matrix has a single
  column constituting the sequence of fitted values corresponding to
  the observations in the given sequence.  The number of rows is the
  number of observations and the entry in row \code{t} is the fitted
  value (conditional expection) corresponding to the observation made
  at time \code{t}. If the model is bivariate (either independent
  or dependent) then the matrix has two columns corresponding
  respectively to the two variables in the bivariate model.

  If the observations are categorical then the nature of the object
  returned changes substantially depending on whether the data are
  univariate, bivariate independent or bivariate dependent.  (See
  \code{\link{hmm}()}.

  In the unvariate case the object corresponding to each sequence
  is a matrix, the number of rows of which is the number of
  observations and the number of columns of which is the number of
  unique \emph{possible} values of the observations. The entry of
  this matrix in row \code{t} and column \code{j} is the conditional
  probability that an emission, at time \code{t}, is equal to
  \eqn{u_i}{u[i]} where \eqn{u_1, \dots, u_m}{u[1], \dots, u[m]}
  are the unique possible values.

  In the bivariate independent case the object is a \emph{list}
  of two matrices, each of which is of the same nature as that
  produced in the univariate case, corresponding respectively to
  the first and second of the two variables.

  In the bivariate dependent case the object is a 3-dimensional
  array of dimension \eqn{m_1 \times m_2 \times n}{m1 x m2 x n}
  where \eqn{m_1}{m1} is the number of unique possible values
  of the first variable, \eqn{m_2}{m2} is the number of unique
  possible values of the second variable, and \eqn{n} is the number
  of observations.  The \code{(i,j,t)-th} entry of this array is
  the conditional probability that an emission, at time \code{t},
  is equal to \eqn{(u_i,v_j)}{(u[i],v[j])} where the \eqn{u_i}{u[i]}
  are the unique possible values of the first variable and the
  \eqn{v_j}{v[j]} are the unique possible values of the second
  variable.

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
