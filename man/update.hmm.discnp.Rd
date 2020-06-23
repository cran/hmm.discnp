\name{update.hmm.discnp}
\alias{update.hmm.discnp}
\title{
   Update a fitted \code{hmm.discnp} model.
}
\description{
   An \code{update()} method for objects of class \code{hmm.discnp}.
}
\usage{
\method{update}{hmm.discnp}(object,..., data, Kplus1=FALSE,
        tpm2=NULL, verbose=FALSE, method=NULL, optimiser=NULL,
        stationary=NULL, mixture=NULL, cis=NULL, tolerance=NULL,
        itmax=NULL, crit=NULL, X=NULL, addIntercept=NULL)
}

\arguments{
  \item{object}{
  An object of class \code{hmm.discnp} as returned by \code{\link{hmm}()}.
}
  \item{\dots}{
  Not used.
}
  \item{data}{
  The data set to which the (updated) model is to be fitted.  See the
  description of the \code{y} argument of \code{\link{hmm}()} for more
  detail.
}
  \item{Kplus1}{
  Logical scalar.  Should the number of states be incremented by 1?
  If so then \code{tpm} (the transition probability matrix) is re-formed
  by \code{rbind()}-ing on a row all of whose entries are \code{1/K}
  (where \code{K} is the \dQuote{old} number of states) and then
  \code{cbind()}-ing on a column of zeroes.  The emission probability
  matrix \code{Rho} is reformed by \code{cbind()}-ing on a column all
  of whose entries are \code{1/m} where \code{m} is the number of discrete
  values of the emissions.

  Note that the intial likelihood of the \dQuote{new} model with
  \code{K+1} states will (should?) be exactly the same as that of
  of the \dQuote{old} fitted \code{K}-state model.

  The \code{Kplus1} argument is provided mainly so as to provide a
  set of starting values for the fitting process which will guarantee
  the log likelihood of a \code{K+1}-state model will be at least as
  large as that of a \code{K}-state model fitted to the same data set.

  Experience indicates that when \code{Kplus1=TRUE} is used, the
  fitting process does not \dQuote{move very far} from the maximum
  log likelihood found for the \code{K}-state model.  It is then
  advisable to try (many) random starting values so as to (try to)
  find the \dQuote{true} maximum for the \code{K+1}-state model.
}
  \item{tpm2}{
  The transtion probability matrix to use when updating a model
  fitted with \code{K=1} and \code{Kplus1=TRUE}. This argument
  is ignored otherwise.  The default value of this argument is
  \code{matrix(0.5,2,2)}.  The value of \code{tpm2} makes no
  difference to the \emph{initial} value of the likelihood
  of the \code{K=2} model (which will be identical to the
  likelihood of the fitted \code{K=1} model that is being
  updated).  Any two-by-two transition probability matrix
  \dQuote{will do}.  However the value of \code{tpm2} could
  conceivably have an impact on the final likelihood of the
  \code{K=2} model to which the fitting procedure converges.
  This is particularly true if the \code{method} is (or is
  switched to) \code{"LM"}.
}
  \item{verbose}{
  See the help for \code{hmm()}.
}
  \item{method}{
  See the help for \code{hmm()}.
}
  \item{optimiser}{
  See the help for \code{hmm()}.
}
  \item{stationary}{
  See the help for \code{hmm()}.
}
  \item{mixture}{
  See the help for \code{hmm()}.
}
  \item{cis}{
  See the help for \code{hmm()}.
}
  \item{tolerance}{
  See the help for \code{hmm()}.
}
  \item{itmax}{
  See the help for \code{hmm()}.
}
  \item{crit}{
  See the help for \code{hmm()}.
}
  \item{X}{
  See the help for \code{hmm()}.
}
  \item{addIntercept}{
  See the help for \code{hmm()}.
}
}
\details{
Except for argument \code{X}, any arguments that are left \code{NULL}
have their values supplied from the \code{args} component of \code{object}.
}
\value{
An object of class \code{hmm.discnp} with an additional component
\code{init.log.like} which is the initial log likelihood
calculated at the starting values of the parameters (which may
be modified from the parameters returned in the object being
updated, if \code{Kplus1} is \code{TRUE}).  The calculation is
done by the function \code{logLikHmm()}.  Barring the strange and
unforeseen, \code{init.log.like} should be (reassuringly) equal
to \code{object$log.like}.  See \code{\link{hmm}()} for details of
the other components of the returned value.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{
\code{\link{hmm}()} \code{\link{rhmm.hmm.discnp}()}
}

\examples{
set.seed(294)
fit  <- hmm(WoodPeweeSong,K=2,rand.start=list(tpm=TRUE,Rho=TRUE),itmax=10)
xxx  <- rhmm(fit,nsim=1)
sfit <- update(fit,data=xxx,itmax=10)
yyy  <- with(SydColDisc,split(y,f=list(locn,depth)))
f1   <- hmm(yyy,K=1)
f2   <- update(f1,data=yyy,Kplus1=TRUE) # Big improvement, but ...
\dontrun{
g2   <- hmm(yyy,K=2) # Substantially better than f2. 
}
}
\keyword{ methods }
\keyword{ models }
