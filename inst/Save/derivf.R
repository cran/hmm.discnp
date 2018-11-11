derivf <- function(Rho,npar) {
K    <- ncol(Rho)
m    <- nrow(Rho)
nst  <- K*(K-1)
d1f  <- array(0,c(m,K,npar))
d2f  <- array(0,c(m,K,npar,npar))
Id   <- diag(m)
den  <- apply(Rho,2,sum)
for(i in 1:m) {
    for(j in 1:K) {
        for(k in 1:(m-1)) {
            h <- (j-1)*(m-1) + k
# d1f[i,j,nst+h] = the derivative of rho_ij w.r.t. theta_{nst+h} (the (nst+h)-th
# parameter) which is equal to phi_{k,j} where h = (j-1)*(m-1) + k.
            d1f[i,j,nst+h] <- Rho[i,j]*(Id[i,k]*den[j] - Rho[k,j])/den[j]^2
# d2f[i,j,h,n] = the second derivative of rho_ij w.r.t. theta_{nst+h} and
# theta_{nst+n}, i.e. w.r.t. phi_{i,k} where h = (k-1)*m + i, and
# phi_{i,n} where n = (l-1)*m + i.
            for(l in 1:(K-1)) {
                n  <- (j-1)*(m-1) + l
                a <- den[j]*Id[i,k]*Id[j,l]
                b <- Rho[l,j]*(Id[i,k]+Id[k,l])+Rho[k,j]*Id[i,l]
                c <- 2*Rho[k,j]*Rho[l,j]/den[j]
                d2f[i,j,nst+h,nst+n] <- Rho[i,j]*(a-b+c)/den[j]^2
            }
        }
    }
}
list(d1f=d1f,d2f=d2f)
}
