\name{logLikHmm}
\alias{logLikHmm}
\title{
Log likelihood of a hidden Markov model
}
\description{
Calculate the log likelihood of a hidden Markov model with
discrete non-parametric observation distributions.
}
\usage{
logLikHmm(y, model=NULL, tpm=NULL, ispd=NULL, Rho=NULL, X=NULL,
          addIntercept=NULL, warn=TRUE)
}
\arguments{
  \item{y}{
  A vector, or list of vectors, or a two column matrix or a a list
  of such matrics, whose entries consist of observations from a
  hidden Markov model with discrete non-parametric observation
  distributions.
}
  \item{model}{An object specifying a hidden Markov model, usually
  returned by \code{\link{hmm}()}.
  }
  \item{tpm}{The transition probability matrix of the Markov chain.
  Ignored (and extracted from \code{model}) if \code{model}
  is non-\code{NULL}.
  }
  \item{ispd}{The vector of probabilities specifying the initial
  state probability distribution, or a matrix each of whose columns
  is a trivial (\dQuote{delta function}) vector specifying the
  \dQuote{most probable} initial state for each observation sequence.
  If \code{ispd} is missing then \code{ispd} is calculated as the
  stationary distribution determined by \code{tpm}.  Ignored (and
  extracted from \code{model}) if \code{model} is non-\code{NULL}.
  }
  \item{Rho}{An object specifying the \dQuote{emission} probabilities
  of the observations.  (See the \bold{Details} in the help for
  \code{\link{hmm}()}.)  Ignored (and extracted from \code{model})
  if \code{model} is non-\code{NULL}.
  }
  \item{X}{An optional \emph{numeric} matrix, or a list of such
  matrices, of \emph{predictors}.  The use of such predictors
  is (currently, at least) applicable only in the univariate
  emissions setting, and then only if \code{Rho} is supplied in
  the \dQuote{newstyle} format.  If \code{X} is a list it must
  be of the same length as \code{y} and all entries of this list
  must have the same number of columns.  The number of rows of each
  entry must be equal to the length of the corresponding entry of
  \code{y}. If \code{X} is a matrix then \code{y} should be a vector
  or one-column matrix (or a list with a single entry equal to such).
  }
  \item{addIntercept}{
  Logical scalar.  See the documentation of \code{\link{hmm}()}.
  If this argument is not specified, and if \code{model} is
  \code{NULL} then an error is thrown.
  }
  \item{warn}{Logical scalar; should a warning be issued if \code{Rho}
  hasn't got relevant dimension names?  (Note that if this is so,
  then the corresponding dimension names are formed from the sorted
  unique values of \code{y} or of the appropriate column(s) of \code{y}.
  And if \emph{this} is so, then the user should be sure that the
  ordering of the entries of \code{Rho} corresponds properly to the
  the sorted unique values of \code{y}.)  This argument is passed
  to the utility function \code{check.yval()} which actually issues
  the warning if \code{warn=TRUE}.
  }
}

\details{
If \code{y} is not provided the function simply returns the
\code{log.like} component of \code{model} (which could be
\code{NULL} if \code{model} was not produced by \code{hmm()}.

The observation values (the entries of the vector or matrix \code{y},
or of the vectors or matrices which constitute the entries of
\code{y} if \code{y} is a list) must be consistent with the
appropriate dimension names of \code{Rho} or of its entries when
\code{Rho} is a list.  More specifically, if \code{Rho} has dimension
names (or its entries have dimension names) then the observation
values must all be found as entries of the appropriate dimension
name vector.  If a vector of dimension names is \code{NULL} then
the corresponding dimension must be equal to the number of unique
observations of the appropriate variate.  integers between \code{1}
and \code{nrow(Rho)}.
}
\value{
The loglikehood of \code{y} given the parameter values specified
in \code{par}.
}
\references{
See \code{\link{hmm}()} for references.
}
\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{
\code{\link{hmm}()}, \code{\link{pr}()}, \code{\link{sp}()}
}
\examples{
# TO DO:  One or more bivariate examples.
P <- matrix(c(0.7,0.3,0.1,0.9),2,2,byrow=TRUE)
R <- matrix(c(0.5,0,0.1,0.1,0.3,
              0.1,0.1,0,0.3,0.5),5,2)
set.seed(42)
lll   <- sample(250:350,20,TRUE)
set.seed(909)
y.num <- rhmm(ylengths=lll,nsim=1,tpm=P,Rho=R,drop=TRUE)
set.seed(909)
y.let <- rhmm(ylengths=lll,nsim=1,tpm=P,Rho=R,yval=letters[1:5],drop=TRUE)
row.names(R) <- 1:5
ll1   <- logLikHmm(y.num,tpm=P,Rho=R)
row.names(R) <- letters[1:5]
ll2   <- logLikHmm(y.let,tpm=P,Rho=R)
ll3   <- logLikHmm(y.let,tpm=P,Rho=R,ispd=c(0.5,0.5))
fit   <- hmm(y.num,K=2,itmax=10)
ll4   <- logLikHmm(y.num,fit) # Use the fitted rather than the "true" parameters.
}
\keyword{ models }
