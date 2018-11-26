derivf <- function(theta,K) {
npar <- length(theta)
m    <- 2 - K + npar/K
npr  <- K*(K-1)
d1f  <- array(0,c(m,K,npar))
d2f  <- array(0,c(m,K,npar,npar))
Id   <- diag(max(K,m))
phi  <- theta[(npr+1):npar]
M    <- rbind(matrix(phi,ncol=K),0)
M    <- t(t(M) - apply(M,2,max))
E    <- exp(M)
den  <- apply(E,2,sum)
mm1  <- m - 1
for(i in 1:m) {
    for(j in 1:K) {
        for(k in 1:mm1) {
            h <- (j-1)*mm1 + k
# d1f[i,j,npr+h] = the derivative of rho_ij w.r.t. theta_{npr+h} (the (npr+h)-th
# parameter) which is equal to phi_{k,j} where h = (j-1)*(m-1) + k.
            d1f[i,j,npr+h] <- E[i,j]*(Id[i,k]*den[j] - E[k,j])/den[j]^2
# d2f[i,j,npr+h,npr+n] = the second derivative of rho_ij w.r.t. theta_{npr+h}
# and theta_{npr+n}, i.e. w.r.t. phi_{k,j} where h = (j-1)*(m-1) + k, and
# phi_{l,j} where n = (j-1)*(m-1) + l.
            for(l in 1:mm1) {
                n <- (j-1)*mm1 + l
                a <- den[j]*Id[i,k]*Id[i,l]
                b <- E[l,j]*(Id[i,k]+Id[k,l])+E[k,j]*Id[i,l]
                c <- 2*E[k,j]*E[l,j]/den[j]
                d2f[i,j,npr+h,npr+n] <- E[i,j]*(a-b+c)/den[j]^2
            }
        }
    }
}
list(d1f=d1f,d2f=d2f)
}
