derivpi <- function(ispd,tpm,npar,dp) {
K     <- length(ispd)
npro  <- K*(K-1)
d1p   <- dp$d1p
d2p   <- dp$d2p
A     <- diag(K)-t(tpm)
A[K,] <- rep(1,K)
d1pi  <- matrix(0,K,npar)
d2pi  <- array(0,c(K,npar,npar))
for(j in 1:npro) {
	v        <- t(d1p[,,j])%*%ispd
	v[K]     <- 0
	d1pi[,j] <- solve(A,v)
}

for(j in 1:npro) {
	for(k in 1:npro) {
		u          <- t(d1p[,,j])%*%d1pi[,k]
		v          <- t(d2p[,,j,k])%*%ispd
		w          <- t(d1p[,,k])%*%d1pi[,j]
		z          <- u + v + w
		z[K]       <- 0
		d2pi[,j,k] <- solve(A,z)
	}
}
list(d1pi=d1pi,d2pi=d2pi)
}
