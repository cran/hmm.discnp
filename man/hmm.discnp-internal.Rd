\name{hmm.discnp-internal}
\alias{ffun}
\alias{init.all}
\alias{recurse}
\alias{revise.ispd}
\alias{revise.rho}
\alias{revise.tpm}
\title{Internal hmm.discnp functions.}
\description{
  Internal hmm.discnp functions.
}
\usage{
ffun(y,Rho)
init.all(nval,K,rand.start,mixture=F)
recurse(fy,tpm,nc,epsilon)
revise.ispd(tpm)
revise.rho(y,gamma,m)
revise.tpm(xi,mixture)
}
\details{
  These functions are auxilliary and are not intended to be called by
  the user.
}
\keyword{internal}
