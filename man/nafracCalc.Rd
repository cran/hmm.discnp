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
nafracCalc(y,drop=TRUE)
}
\arguments{
  \item{y}{
   A vector or a one or two column matrix of discrete data or a
   list of such vectors or matrices, or a \emph{list} of such lists
   (an object of class \code{"multipleHmmDataSets"} such as might
   be produced by \code{\link{rhmm}()}).
}
   \item{drop}{
   Logical scalar. If \code{y} is of class
   \code{"multipleHmmDataSets"} but actually consists of a
   single data set, and if \code{drop} is \code{TRUE}, then the
   returned value is not a list but rather the single component
   that such a list \dQuote{would have had} were \code{drop} equal
   to \code{FALSE}.  This argument is ignored if \code{y} is not of
   class \code{"multipleHmmDataSets"} or has length greater than 1.
}
}
\value{
  If \code{y} is \emph{not} of class \code{"multipleHmmDataSets"},
  then the returned value is a scalar (between 0 and 1) if the data
  are univariate or a pair (2-vector) of such scalars if the data
  are bivariate.  The values are equal to the ratios of the total
  count of missing values in the appropriate column to the total
  number of observations.

  If \code{y} \emph{is} of class \code{"multipleHmmDataSets"},
  and if \code{y} has length greater than 1 or \code{drop} is
  \code{FALSE}, then the returned value is a \emph{list} of such
  scalars or 2-vectors, each corresponding to one of the data sets
  constituting \code{y}.  If \code{y} has length equal to 1 and
  \code{drop} is \code{TRUE}, then the returned value is the same
  as if code{y} were not of class \code{"multipleHmmDataSets"}.
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
