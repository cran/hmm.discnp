\name{nafracCalc}
\alias{nafracCalc}
\title{
  Calculate fractions of missing values.
}
\description{
  Calculate the fraction (univariate case) or fractions (bivariate
  case) of missing values in the data or in each component of the
  data.
}
\usage{
nafracCalc(y)
}
\arguments{
  \item{y}{
   A vector or a one or two column matrix of discrete data or a list of such
   vectors or matrices.
}
}
\value{
  A scalar (between 0 and 1) if the data are univariate or a pair of
  such scalars if the data are bivariate.  The values are equal to the
  ratios of the total count of missing values in the appropriate column
  to the total number of observations.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}
\seealso{
\code{\link{rhmm}()} \code{\link{misstify}()}
}
\examples{
xxx <- with(SydColDisc,split(y,f=list(locn,depth)))
nafracCalc(xxx) # 0.7185199

}

\keyword{utilities}
