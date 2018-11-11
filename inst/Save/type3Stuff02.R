#
# Auxiliary stuff needed when type==3.
#
p2zeta <- function (x) {
    z  <- log(x)
    nok <- z==-Inf
    if(any(nok)) {
        btm <- min(z[!nok])
        z[nok] <- min(btm,-300)
    }
    z <- z - z[length(z)]
    return(z)
}

zeta2p <- function(x){
    m  <- max(x)
    xr <- exp(x-m)
    xr/sum(xr)
}

rho2zeta <- function(Rho) {
K <- dim(Rho)[3]
zeta <- vector("list",K)
for(k in 1:K) {
    Lk <- as.vector(Rho[,,k])
    zeta[[k]] <- p2zeta(Lk)
}
unlist(zeta)
}

zeta2Rho <- function(zeta,ijk){
m      <- ijk[1]
n      <- ijk[2]
K      <- ijk[3]
mn     <- m*n
Rho    <- array(0,dim=c(m,n,K))
istop  <- 0
for(k in 1:K) {
    istart <- 1+istop
    istop  <- istart+mn-1
    rrr   <- zeta2p(zeta[istart:istop])
    Rho[,,k] <- matrix(rrr,nrow=m,ncol=n)
}
Rho
}

msRho <- function(Rho0,G) {
#
# Perform the M-step for Rho numerically.
#
    zeta0 <- rho2zeta(Rho0)
    zeta  <- try(optim(zeta0,ell,method="BFGS",control=list(fnscale=-1),G=G)$par)
## Get rid of this conditional browser when (???) things seem to be working.
#  if(inherits(zeta,"try-error")) browser()
    if(inherits(zeta,"try-error"))  browser()
    ijk    <- dim(G) - c(1,1,0)
    zeta2Rho(zeta,ijk)
}

lse <- function(z){
    z <- as.vector(z)
    m <- max(z)
    z <- z-m
    log(sum(exp(z))) + m
}

ell <- function(zeta,G) {
#
# Expected log likelihood.
#
    ijk <- dim(G)
    m   <- ijk[1]
    n   <- ijk[2]
    K   <- ijk[3]
    zeta <- array(zeta,dim=c(m-1,n-1,K))
    Zk   <- apply(zeta,3,lse)
    Zik  <- apply(zeta,c(1,3),lse)
    Zjk  <- apply(zeta,c(2,3),lse)
# Part 1:
    A  <- sum(G[-m,-n,,drop=FALSE]*zeta)
    B  <- apply(G[-m,-n,,drop=FALSE],3,sum)
    P1 <- A - sum(B*Zk)
# Part 2:
    A  <- G[-m,n,]*Zik
    B  <- apply(G[-m,n,,drop=FALSE],3,sum)
    P2 <- sum(A) - sum(B*Zk)
# Part 3:
    A  <- G[m,-n,]*Zjk
    B  <- apply(G[m,-n,,drop=FALSE],3,sum)
    P3 <- sum(A) - sum(B*Zk)
# Aw' done!
    P1+P2+P3
}
