mps <- function(y,object=NULL,tpm,Rho,ispd,yval=NULL) {
#
# Function mps: most probable states.
#

y  <- tidyup(y,yval)
nc <- ncol(y)

if(!is.null(object)) {
	tpm  <- object$tpm
	Rho  <- object$Rho
	ispd <- object$ispd
}
fy   <- ffun(y,Rho)
rp   <- recurse(fy, tpm, ispd, nc)
rslt <- apply(rp$gamma, 2, which.max)
if(nc==1) rslt else matrix(rslt,ncol=nc)
}
