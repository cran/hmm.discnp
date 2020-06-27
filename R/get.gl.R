get.gl <- function(theta,K,y) {
#
# Function get.gl --- get gradient and log likelihood.
#
npar <- length(theta)
tpm  <- getTpm(theta,K,stationary=TRUE)
ispd <- revise.ispd(tpm)
rrr  <- attr(y,"lvls")
Rho  <- getRho(theta,K,rhovals=rrr,stationary=TRUE,
               prednames="Intercept")
m    <- length(rrr)
dp   <- derivp(theta,K)
d1p  <- dp$d1p
dpi  <- derivpi(ispd,tpm,npar,dp)
d1pi <- dpi$d1pi
dfun <- derivf(theta,K)
d1f  <- dfun$d1f

# Run through the list "y":
ndat  <- length(y)
alist <- vector("list",ndat)
xll   <- numeric(ndat)
ky    <- 0     
fy    <- ffun(y,Rho,type=1)
j2    <- 0
for(yl in y) {
    ylv   <- yl[,1]
    ny    <- length(ylv)
    ymiss <- is.na(ylv)
    j1    <- j2+1
    j2    <- j2+ny
    xxx   <- .Fortran(
		"getgl",
		NAOK=TRUE,
		fy=as.double(fy[,j1:j2]),
		y=as.integer(ylv),
                ymiss=as.integer(ymiss),
		tpm=as.double(tpm),
		xispd=as.double(ispd),
		d1pi=as.double(d1pi),
		kstate=as.integer(K),
		n=as.integer(ny),
		npar=as.integer(npar),
		d1p=as.double(d1p),
                m=as.integer(m),
		d1f=as.double(d1f),
		alpha=double(K),
		alphw=double(K),
		a=double(K*npar),
		aw=double(K*npar),
		xlc=double(ny),
		PACKAGE="hmm.discnp"
		)
    ky <- ky + 1
    alist[[ky]] <- matrix(xxx$a,K,npar)/xxx$xlc[ny]
    xll[ky]  <- sum(log(xxx$xlc))
} 
a     <- array(unlist(alist),c(K,npar,ndat))
grad  <- apply(a,2,sum)
ll    <- sum(xll)
list(ll=ll,grad=grad)
}
