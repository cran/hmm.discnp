rgethgl <- function(fy,y,ymiss,tpm,xispd,d1pi,d2pi,npar,d1p,d2p,m,d1f,d2f) {
#
# Function rgethgl --- get Hessian, gradient, and log likelihood.
# for one observation sequence, (new) raw R method.
#

# Create blank arrays.
kstate <- nrow(tpm)
n      <- length(y)
alpha  <- numeric(kstate)
xlc    <- numeric(n)
a      <- matrix(NA,kstate,npar)
b      <- array(NA,dim=c(kstate,npar,npar))

alpha  <- xispd*fy[,1]
sxlc   <- sum(alpha)
xlc[1] <- sxlc
alpha  <- alpha/sxlc
d1fx   <- if(is.na(y[1])) 0 else d1f[y[1],,]
a      <- xispd*d1fx + fy[,1]*d1pi
for(j in 1:kstate) {
    d1fx   <- if(is.na(y[1])) rep(0,npar) else d1f[y[1],j,]
    d2fx   <- if(is.na(y[1])) matrix(0,npar,npar) else d2f[y[1],j,,]
    b[j,,] <- (xispd[j]*d2fx + as.matrix(d1fx)%*%d1pi[j,]
               + as.matrix(d1pi[j,])%*%d1fx + fy[j,1]*d2pi[j,,])
}

strangeProd <- function(A,B) {
# Due to Chuck Berry and Jeff Newmiller 05/08/2018
  nca <- ncol(A)
  ncb <- ncol(B)
  j.index <- rep(seq.int(nca),times=ncb)
  k.index <- rep(seq.int(ncb),each=nca)
  array(A[,j.index]*B[,k.index],c(nrow(A),nca,ncb))
}


if(n>1) {
for(kt in 2:n) {
# Do the b's:
    prevb <- b
    for(j in 1:kstate) {
        Part1  <- apply(alpha*d2p[,j,,],c(2,3),sum)
        Part2  <- apply(strangeProd(a,d1p[,j,]),c(2,3),sum)/sxlc
        Part3  <- apply(strangeProd(d1p[,j,],a),c(2,3),sum)/sxlc
        Part4  <- apply(tpm[,j]*prevb,c(2,3),sum)/sxlc
        Part5  <- alpha%*%d1p[,j,]
        Part6  <- tpm[,j]%*%a/sxlc
        Part7  <- alpha%*%tpm[,j]
        Mult1  <- fy[j,kt]
        Mult2  <- if(is.na(y[kt])) rep(0,npar) else d1f[y[kt],j,]
        Mult3  <- if(is.na(y[kt])) matrix(0,npar,npar) else d2f[y[kt],j,,]
        b[j,,] <- (Mult1*(Part1 + Part2 + Part3 + Part4) +
                  as.matrix(Mult2)%*%(Part5 + Part6) +
                  t(Part5 + Part6)%*%Mult2 +
                  Mult3*Part7[1,1])
    }

# Do the a's:
    preva <- a
    for(j in 1:kstate) {
        Part1 <- alpha%*%d1p[,j,]
        Part2 <- tpm[,j]%*%preva/sxlc
        Part3 <- alpha%*%tpm[,j]
        Mult1 <- fy[j,kt]
        Mult2 <- if(is.na(y[kt])) rep(0,npar) else d1f[y[kt],j,]
        a[j,] <- Mult1*(Part1 + Part2) + Mult2*Part3[1,1]
    }

# Do the alpha's:
    alpha   <- as.vector(fy[,kt]*alpha%*%tpm)
    sxlc    <- sum(alpha)
    alpha   <- alpha/sxlc
    xlc[kt] <- sxlc
    }
}

# Finish off:
ll   <- sum(log(xlc))
grad <- apply(a,2,sum)/sxlc
hess <- apply(b,c(2,3),sum)/sxlc - as.matrix(grad)%*%grad

list(ll=ll,grad=grad,hess=hess)
}
