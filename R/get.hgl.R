get.hgl <- function(theta,K,y,method=c("fortran","oraw","raw")) {
#
# Function get.hgl --- get Hessian, gradient, and log likelihood.
#
meth <- match.arg(method)
npar <- length(theta)
tpm  <- getTpm(theta,K,stationary=TRUE)
ispd <- revise.ispd(tpm)
rrr  <- attr(y,"lvls")
Rho  <- getRho(theta,K,rhovals=rrr,stationary=TRUE,
               prednames="Intercept")
m    <- nrow(Rho)
dp   <- derivp(theta,K)
d1p  <- dp$d1p
d2p  <- dp$d2p
dpi  <- derivpi(ispd,tpm,npar,dp)
d1pi <- dpi$d1pi
d2pi <- dpi$d2pi
dfun <- derivf(theta,K)
d1f  <- dfun$d1f
d2f  <- dfun$d2f

# Run through the list "y":
ndat  <- length(y)
lll   <- numeric(ndat)
glist <- vector("list",ndat)
hlist <- vector("list",ndat)
ky    <- 0
fy    <- ffun(y,Rho,type=1)
j2    <- 0
for(yl in y) {
    ylv          <- as.integer(yl[,1])
    ymiss        <- is.na(ylv)
    ny           <- length(ylv)
    j1           <- j2 + 1
    j2           <- j2 + ny
    fyl          <- fy[,j1:j2]
    workfun      <- switch(meth,fortran=forgethgl,oraw=orgethgl,raw=rgethgl)
    xxx          <- workfun(fyl,ylv,ymiss,tpm,ispd,d1pi,d2pi,npar,d1p,d2p,m,d1f,d2f)
    ky           <- ky + 1
    lll[ky]      <- xxx$ll
    glist[[ky]]  <- xxx$grad
    hlist[[ky]]  <- xxx$hess
} 
ll    <- sum(lll)
gpre  <- matrix(unlist(glist),nrow=npar,ncol=ndat)
grad  <- apply(gpre,1,sum)
hpre  <- array(unlist(hlist),c(npar,npar,ndat))
hess  <- apply(hpre,c(1,2),sum)
list(ll=ll,grad=grad,hess=hess)
}
