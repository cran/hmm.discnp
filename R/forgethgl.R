forgethgl <- function(fy,y,ymiss,tpm,ispd,d1pi,d2pi,
                       npar,d1p,d2p,m,d1f,d2f) {
#
# Function forgethgl --- get Hessian, gradient, and log likelihood.
# for one observation sequence, Fortran method.
#

ny    <- length(y)
K     <- nrow(tpm)
xxx   <- .Fortran(
                "gethgl",
                NAOK=TRUE,
                fy=as.double(fy),
                y=y,
                ymiss=as.integer(ymiss),
                tpm=as.double(tpm),
                xispd=as.double(ispd),
                d1pi=as.double(d1pi),
                d2pi=as.double(d2pi),
                kstate=as.integer(K),
                n=as.integer(ny),
                npar=as.integer(npar),
                d1p=as.double(d1p),
                d2p=as.double(d2p),
                m=as.integer(m),
                d1f=as.double(d1f),
                d2f=as.double(d2f),
                alpha=double(K),
                alphw=double(K),
                a=double(K*npar),
                b=double(K*npar*npar),
                aw=double(K*npar),
                bw=double(K*npar*npar),
                xlc=double(ny),
                ll=double(1),
                grad=double(npar),
                hess=double(npar*npar),
                PACKAGE="hmm.discnp"
       )
ll    <- xxx$ll
grad  <- xxx$grad
hess  <- matrix(xxx$hess,npar,npar)
list(ll=ll,grad=grad,hess=hess)
}
