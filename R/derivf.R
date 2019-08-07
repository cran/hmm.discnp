derivf <- function(theta,K) {
npar <- length(theta)
m    <- 2 - K + npar/K
npr  <- K*(K-1)
d1f  <- array(0,c(m,K,npar))
d2f  <- array(0,c(m,K,npar,npar))

# d1f[i,j,k] = the derivative of rho_ij w.r.t. theta_k
# where theta_k = phi_{i1,j1}, where k = (j1-1)*(m-1) + i1 + npr
# These derivatives are identically zero unless j1 = j since
# rho_ij does not involve phi_{i1,j1} unless j1 = j.

# d2f[i,j,k,l] = the second derivative of rho_ij w.r.t. theta_k
# and theta_l, where theta_k = phi_{i1,j1} and theta_l = phi_{i2,j2},
# k = (j1-1)*(m-1) + i1 + npr and l = (j2-1)*(m-1) + i2 + npr
# These derivatives are identically zero unless j1 = j and j2 = j.

Id   <- diag(m)
phi  <- theta[(npr+1):npar]
M    <- rbind(matrix(phi,ncol=K),0)
M    <- t(t(M) - apply(M,2,max))
E    <- exp(M)
den  <- apply(E,2,sum)
mm1  <- m - 1
for(i in 1:m) {
    for(j in 1:K) {
        for(i1 in 1:mm1) {
            k <- (j-1)*mm1 + i1 + npr
            d1f[i,j,k] <- E[i,j]*(Id[i,i1]*den[j] - E[i1,j])/den[j]^2
            for(i2 in 1:mm1) {
                l <- (j-1)*mm1 + i2 + npr
                a <- den[j]*Id[i,i1]*Id[i,i2]
                #b <- E[i2,j]*(Id[i,i1]+Id[i1,i2])+E[i1,j]*Id[i1,i2]
                b <- E[i1,j]*(Id[i,i2]+Id[i1,i2])+E[i2,j]*Id[i,i1]
                c <- 2*E[i1,j]*E[i2,j]/den[j]
                d2f[i,j,k,l] <- E[i,j]*(a-b+c)/den[j]^2
            }
        }
    }
}
list(d1f=d1f,d2f=d2f)
}
