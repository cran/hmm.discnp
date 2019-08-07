\name{hmm.discnp-internal}
\alias{[.multipleHmmDataSets}
\alias{check.yval}
\alias{checkConstColumns}
\alias{checkyXoK}
\alias{derivf}
\alias{derivp}
\alias{derivpi}
\alias{ell}
\alias{expForm2p}
\alias{ffun}
\alias{fix.tpm}
\alias{forgethgl}
\alias{orgethgl}
\alias{rgethgl}
\alias{get.hgl}
\alias{get.gl}
\alias{get.hgl}
\alias{get.l}
\alias{getIspd}
\alias{getRho}
\alias{getTpm}
\alias{hmmBD}
\alias{hmmBI}
\alias{hmmNumOpt}
\alias{hmmLM}
\alias{hmmSD}
\alias{hmmUV}
\alias{init.all}
\alias{lmstep}
\alias{lse}
\alias{msRho}
\alias{makeDat}
\alias{paramExtract}
\alias{p2expForm}
\alias{phi2rho}
\alias{predictEngineHmmDiscnp}
\alias{recurse}
\alias{reparam}
\alias{revise.ispd}
\alias{revise.rho}
\alias{revise.tpm}
\alias{rho2phi}
\alias{steepest}
\alias{tidyList}
\alias{zeta2p}
\alias{zeta2Rho}
\title{Internal hmm.discnp functions.}
\description{
  Internal hmm.discnp functions.
}
\usage{
\method{[}{multipleHmmDataSets}(x,i)
check.yval(yval, Rho, type, warn=TRUE)
checkConstColumns(y,prednames)
checkyXoK(y, X)
derivf(theta, K)
derivp(theta, K)
derivpi(ispd, tpm, npar, dp)
ell(phi, G)
expForm2p(x)
ffun(Dat, Rho, type)
fix.tpm(x, K)
forgethgl(fy, y, ymiss, tpm, ispd, d1pi, d2pi, npar, d1p,
          d2p, m, d1f, d2f)
orgethgl(fy, y, ymiss, tpm, xispd, d1pi, d2pi, npar, d1p,
          d2p, m, d1f, d2f)
rgethgl(fy, y, ymiss, tpm, xispd, d1pi, d2pi, npar, d1p,
          d2p, m, d1f, d2f)
get.gl(theta, K, y)
get.hgl(theta, K, y, method=c("fortran","oraw","raw"))
get.l(theta, K, y)
getIspd(pars, K)
getRho(pars, K, rhovals, stationary, prednames)
getTpm(pars, K, stationary)
hmmBD(y, yval=NULL, par0=NULL, K=NULL, rand.start=NULL, stationary=cis,
      mixture=FALSE, cis=TRUE, tolerance, digits=NULL, verbose=FALSE,
      itmax=200, crit, bicm)
hmmBI(y, yval=NULL, par0=NULL, K=NULL, rand.start=NULL, stationary=cis,
      mixture=FALSE, cis=TRUE, tolerance, digits=NULL, verbose=FALSE,
      itmax=200, crit, bicm)
hmmNumOpt(Dat, par0, stationary, verbose, itmax, bicm, rhovals, npar, optimiser,
          optimMethod, hessian=FALSE, ...)
hmmLM(y, par0, itmax=200, crit, lmc=10, tolerance,
      bicm, rhovals, hglmethod, digits=NULL, verbose=FALSE) 
hmmSD(y, par0, itmax=200, crit, tolerance,
      bicm, rhovals, hglmethod, digits=NULL, verbose=FALSE) 
hmmUV(y, yval=NULL, par0=NULL, K=NULL, rand.start=NULL,
      method=c("EM","bf","LM","SD"), hglmethod=c("fortran","oraw","raw"),
      optimiser=c("nlm","optim"), optimMethod=NULL, stationary=cis,
      mixture=FALSE, cis=TRUE, tolerance, digits=NULL, verbose=FALSE,
      itmax=200, crit, bicm, newstyle=TRUE, X=NULL, addIntercept=TRUE,
      lmc=10, hessian=FALSE, ...)
init.all(nval, K, rand.start, mixture, indep, newstyle, yval, prednames)
lmstep(theta, K, y, lmc, hglmethod)
lse(z)
msRho(Rho0, G)
makeDat(y, X, addIntercept)
paramExtract(Rho, newstyle)
p2expForm(x)
phi2rho(phi, ijk)
predictEngineHmmDiscnp(stateProbs,Rho,numb,drop)
recurse(fy, tpm, ispd, lns)
reparam(pars, stationary)
revise.ispd(tpm=NULL, gamma=NULL, lns=NULL, cis=TRUE)
revise.rho(Dat, gamma, stnms, type)
revise.tpm(xi, mixture)
rho2phi(Rho)
steepest(K, y, theta)
tidyList(y, rp=c("response","predictor"),addIntercept)
}
\details{
  These functions are auxiliary and are not intended to be called
  by the user.
}
\keyword{internal}
