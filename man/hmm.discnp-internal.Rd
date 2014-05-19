\name{hmm.discnp-internal}
\alias{ffun}
\alias{init.all}
\alias{check.yval}
\alias{recurse}
\alias{revise.ispd}
\alias{revise.rho}
\alias{revise.tpm}
\alias{charList}
\title{Internal hmm.discnp functions.}
\description{
  Internal hmm.discnp functions.
}
\usage{
ffun(y,Rho)
init.all(nval,K,rand.start,mixture=FALSE)
check.yval(y,Rho,warn=TRUE)
recurse(fy,tpm,ispd,lns)
revise.ispd(tpm=NULL,gamma=NULL,lns=NULL,cis=TRUE)
revise.rho(y,gamma,yval)
revise.tpm(xi,mixture)
charList(y)
}
\details{
  These functions are auxilliary and are not intended to be called by
  the user.
}
\keyword{internal}
