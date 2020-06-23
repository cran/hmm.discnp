get.l <- function(theta,K,y) {
#
# Function get.l --- get log likelihood.
#
tpm  <- getTpm(theta,K,stationary=TRUE)
ispd <- revise.ispd(tpm)
rrr  <- attr(y,"lvls")
Rho  <- getRho(theta,K,rhovals=rrr,stationary=TRUE,
               prednames="Intercept")

# Run through the list "y":
ndat <- length(y)
xll  <- numeric(ndat)
ky   <- 0
fy   <- ffun(y,Rho,type=1)
j2   <- 0
for(yl in y) {
    ylv <- yl[,1]
    ny  <- length(ylv)
    j1  <- j2 + 1
    j2  <- j2 + ny
    xxx <- .Fortran(
		"getl",
		fy=as.double(fy[,j1:j2]),
		tpm=as.double(tpm),
		xispd=as.double(ispd),
		kstate=as.integer(K),
		n=as.integer(ny),
		alpha=double(K),
		alphw=double(K),
		xlc=double(ny),
		PACKAGE="hmm.discnp"
		)
    ky <- ky + 1
    xll[ky] <- sum(log(xxx$xlc))
} 
sum(xll)
}
