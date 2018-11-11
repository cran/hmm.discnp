msRho <- function(Rho0,G) {
#
# Perform the M-step for Rho numerically.
#
    phi0 <- rho2phi(Rho0)
    phi  <- try(optim(phi0,ell,method="BFGS",control=list(fnscale=-1),G=G)$par)
## Get rid of this conditional browser when (???) things seem to be working.
    if(inherits(phi,"try-error")) {
        if(interactive()) browser() else stop("Problem with optim().\n")
    }
    ijk <- dim(G) - c(1,1,0)
    phi2rho(phi,ijk)
}

rho2phi <- function(Rho) {
K <- dim(Rho)[3]
phi <- vector("list",K)
for(k in 1:K) {
    Lk <- as.vector(Rho[,,k])
    phi[[k]] <- p2expForm(Lk)
}
unlist(phi)
}

phi2rho <- function(phi,ijk){
m      <- ijk[1]
n      <- ijk[2]
K      <- ijk[3]
mn     <- m*n
Rho    <- array(0,dim=c(m,n,K))
istop  <- 0
for(k in 1:K) {
    istart <- 1+istop
    istop  <- istart+mn-1
    rrr   <- expForm2p(phi[istart:istop])
    Rho[,,k] <- matrix(rrr,nrow=m,ncol=n)
}
Rho
}

ell <- function(phi,G) {
#
# Expected log likelihood.
#
    ijk <- dim(G)
    m   <- ijk[1]
    n   <- ijk[2]
    K   <- ijk[3]
    phi <- array(phi,dim=c(m-1,n-1,K))
    Zk   <- apply(phi,3,lse)
    Zik  <- apply(phi,c(1,3),lse)
    Zjk  <- apply(phi,c(2,3),lse)
# Part 1:
    A  <- sum(G[-m,-n,,drop=FALSE]*phi)
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

lse <- function(z){
    z <- as.vector(z)
    m <- max(z)
    z <- z-m
    log(sum(exp(z))) + m
}
