recurse <- function(fy,tpm,nc,epsilon)
{
#
# Function recurse to revise the ``recursive probabilities'',
# given the parameters theta, and the observations y.
#

# Get the initial state probability distribution.
ispd <- revise.ispd(tpm)

# Set a bunch of constants:
K  <- nrow(tpm)
K2 <- K*K
L  <- ncol(fy)
M  <- K*L
N  <- K*M - K2
n  <- L/nc

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
