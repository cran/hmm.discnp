\name{anova.hmm.discnp}
\alias{anova.hmm.discnp}
\title{
    Anova for hmm.discnp models
}
\description{
    Performs a likelihood ratio test to compare two discrete non-parametric
    hidden Markov models.
}
\usage{
\method{anova}{hmm.discnp}(object, ...)
}
\arguments{
  \item{object}{
  An object of class \dQuote{hmm.discnp} as returned by the
  function \code{\link{hmm}()}.
}
  \item{\dots}{
  A second object of class \dQuote{hmm.discnp}.  There must be
  only \emph{one} such object.
}
}

\value{
  A list with entries
    \item{stat}{The likelihood ratio statistic.}
    \item{df}{The degrees of freedom.}
    \item{pvalue}{The p-value of the test.}

  This list has an attribute \dQuote{details} which is a vector
  consisting of the first and second log likelihoods and the associated
  numbers of parameters, in order of these numbers of parameters.
  (See \bold{Warning}.)
}

\section{Warning}{
   Hidden Markov models can be numerically delicate and the fitting
   algorithm can converge to a local maximum of the likelihood surface
   which is not the global maximum.  Thus it is entirely possible
   for the log likelihood of the model with the greater number of
   parameters to be \emph{smaller} than that of the model with the
   lesser number of parameters.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{
   \code{\link{hmm}()}
}
\examples{
xxx  <- with(SydColDisc,split(y,f=list(locn,depth)))
fit1 <- hmm(xxx,K=1,itmax=10)
fit2 <- hmm(xxx,K=2,itmax=10)
anova(fit1,fit2)
}

\keyword{models}
\keyword{methods}
