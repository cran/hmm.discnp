\name{hmm.discnp-internal}
\alias{[.multipleHmmDataSets}
\alias{check.yval}
\alias{checkConstColumns}
\alias{checkStartVal}
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
\alias{simference}
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
checkStartVal(par0,K,indep,yval,rand.start,mixture,prednames)
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
get.hgl(theta, K, y, hglmethod)
get.l(theta, K, y)
getIspd(pars, K)
getRho(pars, K, rhovals, stationary, prednames)
getTpm(pars, K, stationary)
hmmBD(y, par0, K, stationary,
      mixture, cis, tolerance, digits, verbose,
      itmax, crit, bicm)
hmmBI(y, par0, K, stationary,
      mixture, cis, tolerance, digits, verbose,
      itmax, crit, bicm)
hmmNumOpt(Dat, par0, stationary, verbose, itmax, bicm, rhovals, npar, optimiser,
          optimMethod, hessian=FALSE, ...)
hmmLM(y, par0, itmax=200, crit, lmc=10, tolerance,
      bicm, rhovals, hglmethod, digits=NULL, verbose=FALSE) 
hmmSD(y, par0, itmax=200, crit, tolerance,
      bicm, rhovals, hglmethod, digits=NULL, verbose=FALSE) 
hmmUV(y, par0, K, method,
      hglmethod, optimiser, optimMethod, stationary,
      mixture, cis, tolerance, digits, verbose,
      itmax, crit, bicm, X, addIntercept, lmc, hessian, ...)
init.all(K, rand.start, mixture, indep, yval, prednames)
lmstep(theta, K, y, lmc, hglmethod)
lse(z)
msRho(Rho0, G)
makeDat(y, X, addIntercept)
paramExtract(Rho)
p2expForm(x)
phi2rho(phi, ijk)
predictEngineHmmDiscnp(stateProbs,Rho,numb,drop)
recurse(fy, tpm, ispd, lns)
reparam(object, expForm=TRUE, stationary=NULL)
revise.ispd(tpm=NULL, gamma=NULL, lns=NULL, cis=TRUE)
revise.rho(Dat, gamma, type)
revise.tpm(xi, mixture)
rho2phi(Rho)
simference(object, expForm=TRUE, seed=NULL,nsim=100,verbose=TRUE)
steepest(K, y, theta)
tidyList(y, rp=c("response","predictor"),addIntercept=NULL, yval=NULL)
}
\details{
  These functions are auxiliary and are not intended to be called
  by the user.
}
\keyword{internal}
