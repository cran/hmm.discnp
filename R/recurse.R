recurse <- function(fy,tpm,ispd,nc)
{
#
# Function recurse to revise the ``recursive probabilities'',
# given the parameters theta, and the observations y.
#

# Set a bunch of constants:
K  <- nrow(tpm)
K2 <- K*K
L  <- ncol(fy)
M  <- K*L
N  <- K*M - K2
n  <- L/nc
epsilon <- 10*.Machine$double.eps

# Recursive probabilities:

        rp <- .Fortran(
                'recurse',
		fy=as.double(fy),
                xispd=as.double(ispd),
                tpm=as.double(tpm),
		nc=as.integer(nc),
                epsilon=as.double(epsilon),
		n=as.integer(n),
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
