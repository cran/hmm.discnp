recurse <- function(fy,tpm,ispd,lns)
{
#
# Function recurse to calculate the ``recursive probabilities'',
# given the parameters theta, and the observations y.
#

# Set a bunch of constants:
K  <- nrow(tpm)
K2 <- K*K
L  <- ncol(fy)
M  <- K*L
N  <- K*M - K2
nreps <- length(lns)
epsilon <- sqrt(.Machine$double.eps)

# Recursive probabilities:

        rp <- .Fortran(
                'recurse',
		fy=as.double(fy),
                xispd=as.double(ispd),
                tpm=as.double(tpm),
		nreps=as.integer(nreps),
                epsilon=as.double(epsilon),
		lns=as.integer(lns),
                nstate=as.integer(K),
                wrk=double(K2),
		xlc=double(L),
                alpha=double(M),
                beta=double(M),
                gamma=double(M),
                xi=double(N),
		xisum=double(K2),
                PACKAGE="hmm.discnp"
        )

list(gamma=matrix(rp$gamma,nrow=K),xi=matrix(rp$xisum,nrow=K),llc=rp$xlc)
}
