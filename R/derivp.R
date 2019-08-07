derivp <- function(theta,K) {

npar <- length(theta)
d1p  <- array(0,c(K,K,npar))
d2p  <- array(0,c(K,K,npar,npar))

# d1p[i,j,k] = the derivative of p_ij = P[i,j] with respect
# to theta_k = where theta_k = zeta_{i1,j1}, k = (j1-1)*K+i1.
# These derivatives are identically zero unless i1 = i since
# p_ij does not involve zeta_{i1,j1} unless i1 = i.

# d2p[i,j,k,l] = the second derivative of p_ij with respect to
# zeta_{i1,j1} and zeta_{i2,j2}, k = (j1-1)*K + i1, and
# l = (j2-1)*K + i2.  These derivatives are identically zero
# unless i1 = i and i2 = i.
 
zeta <- theta[1:(K*(K-1))]
M    <- cbind(matrix(zeta,nrow=K),0)
M    <- M-apply(M,1,max)
E    <- exp(M)
Id   <- diag(K)
den  <- apply(E,1,sum)
Km1  <- K-1
for(i in 1:K) {
    for(j in 1:K) {
        for(j1 in 1:Km1) {
            k <- i + K*(j1-1)
            d1p[i,j,k] <- E[i,j]*(Id[j,j1]*den[i] - E[i,j1])/den[i]^2
            for(j2 in 1:Km1) {
                l <- i + K*(j2-1)
                a <- den[i]*Id[j,j1]*Id[j,j2]
                b <- E[i,j2]*Id[j,j1]+E[i,j1]*(Id[j,j2]+Id[j1,j2])
                c <- 2*E[i,j1]*E[i,j2]/den[i]
                d2p[i,j,k,l] <- E[i,j]*(a-b+c)/den[i]^2
            }
        }
    }
}
list(d1p=d1p,d2p=d2p)
}
