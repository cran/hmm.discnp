\name{cnvrtRho}
\alias{cnvrtRho}
\title{
    Convert Rho between forms.
}
\description{
    Converts a matrix specification of the emission probabilities
    (in which the probabilities are simply the entries of the matrix)
    to a data frame specification (in which the probabilities are
    a logistic-style function of the parameters) or vice versa.
}
\usage{
cnvrtRho(Rho)
}
\arguments{
  \item{Rho}{A specification of the emission probabilities of a
  discrete valued hidden Markov model.  It may be either a
  matrix of these probabilities, in which case it is converted
  to a three column data frame, or it may be a three column
  data frame, in which case it is converted to a matrix
  of probabilities.  See \code{\link{hmm}()} for more details
  about the structure of \code{Rho}, in either form.
}
}
\details{
  The data frame specification of \code{Rho} allows
  for predictor variables \code{x}.  If \code{Rho} is of the
  data frame form, and is designed to allow for predictor
  variables, then it will have more than three columns and
  cannot be converted to the matrix form.  In such
  cases \code{cnvrtRho} will throw an error.
}
\value{
   A data frame if the argument \code{Rho} is a matrix,
    or a matrix if the argument \code{Rho} is a data
   frame.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{
    \code{\link{hmm}()}
}
\examples{
Yval <- LETTERS[1:10]
Tpm  <- matrix(c(0.75,0.25,0.25,0.75),ncol=2,byrow=TRUE)
Rho  <- cbind(c(rep(1,5),rep(0,5)),c(rep(0,5),rep(1,5)))/5
rownames(Rho) <- Yval
newRho <- cnvrtRho(Rho)
oldRho <- cnvrtRho(newRho)
}
\keyword{ utilities }
