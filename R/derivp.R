derivp <- function(theta,K) {

npar <- length(theta)
d1p  <- array(0,c(K,K,npar))
d2p  <- array(0,c(K,K,npar,npar))

zeta <- theta[1:(K*(K-1))]
E    <- exp(cbind(matrix(zeta,nrow=K),0))
Id   <- diag(K)
den  <- apply(E,1,sum)
for(i in 1:K) {
    for(j in 1:K) {
        for(k in 1:(K-1)) {
            m <- (K-1)*(i-1) + k
            d1p[i,j,m] <- E[i,j]*(Id[j,k]*den[i] - E[i,k])/den[i]^2
            for(l in 1:(K-1)) {
                n <- (K-1)*(i-1) + l 
                a <- den[i]*Id[j,k]*Id[j,l]
                b <- E[i,l]*Id[j,k]+E[i,k]*(Id[j,l]+Id[k,l])
                c <- 2*E[i,k]*E[i,l]/den[i]
                d2p[i,j,m,n] <- E[i,j]*(a-b+c)/den[i]^2
            }
        }
    }
}
list(d1p=d1p,d2p=d2p)
}
