revise.rho <- function(Dat,gamma,stnms,type) {

lvls <- attr(Dat,"lvls")

if(type==1) {
    DF  <- do.call(rbind,Dat)
    K   <- nrow(gamma)
    xxx <- vector("list",K)
    for(k in 1:K) {
        mnfit    <- multinom(y ~ 0 + ., data=DF,weights=gamma[k,],trace=FALSE)
#
# The nnet::multinom() function uses the convention that the *first*
# exponent is 0, rather than the last, whereas the convention that is
# used elsewhere in this package is that the *last* exponent is set to 0.
# This messes up the procedures applied by the reparam() and getRho()
# functions.
# We can rectify the situation by subtracting the last row of
# the matrix produced from all rows (including itself), making
# the last row equal to the zero vector.
# In the simple intercept-only setting this could be done by:
#            x <- c(0,coef(mnfit))
#            xxx[[k]] <- x - x[length(x)]
# It is only slightly more complicated in the general case.
# Note that coef(mnfit) is always a matrix, even if there is
# only a single column.
        M <- rbind(0,coef(mnfit))
        xxx[[k]] <- t(t(M) - M[nrow(M),])
    }
    CCC <- do.call(rbind,xxx)
    rownames(CCC) <- NULL
    colnames(CCC) <- attr(Dat,"prednames")
    Rho <- data.frame(y=factor(rep(lvls,K),levels=lvls),
                      state=factor(rep(1:K,each=length(lvls))),CCC)
    return(Rho)
}

if(type==2) {
    yf  <- factor(unlist(lapply(Dat,function(x){x$y})),levels=lvls)
    Rho <- apply(gamma,1,function(x,y){xtabs(x ~ y)},y=yf)
    Rho <- t(t(Rho)/apply(Rho,2,sum))
    rownames(Rho) <- lvls
    return(Rho)
}

if(type==3) {
    Rho <- vector("list",2)
    for(j in 1:2) {
	yj   <- factor(unlist(lapply(Dat,function(x){x[,j]})),levels=lvls[[j]])
        Rhoj <- apply(gamma,1,function(x,y){xtabs(x ~ y)},y=yj)
	Rho[[j]] <- t(t(Rhoj)/apply(Rhoj,2,sum))
        rownames(Rho[[j]]) <- lvls[[j]]
    }
    return(Rho)
}

if(type==4) {
# Here Rho must be a 3 dimensional array.  We shall set it up
# so that the third dimension ("layers") corresponds to "state".
# Each layer is a matrix whose (i,j)th entry is the estimated
# probability that X = x_i and Y = y_j where X and Y are the
# two variables that are emitted and where x_i and y_j are the
# possible values of X and Y respectfully.

aNA  <- any(is.na(unlist(Dat)))
ym   <- do.call(rbind,Dat)
Rho0 <- array(0,dim=c(sapply(lvls,length),length(stnms)))
G    <- array(0,dim=dim(Rho0)+c(1,1,0))
m    <- dim(G)[1]    
n    <- dim(G)[2]    
K    <- dim(G)[3]
X    <- factor(ym[,1],levels=c(lvls[[1]],NA),exclude=NULL)
Y    <- factor(ym[,2],levels=c(lvls[[2]],NA),exclude=NULL)
for(k in 1:K) {
    G[,,k] <- xtabs(gamma[k,] ~ X + Y,exclude=NULL)
    Rho0[,,k] <- G[-m,-n,k]
    Rho0[,,k] <- Rho0[,,k]/sum(Rho0[,,k])
}
if(aNA) {
# Call upon optim() to find the optimum value of Rho,
# Using Rho0 as starting values.
    Rho <- msRho(Rho0,G)
} else {
    Rho <- Rho0
}
dimnames(Rho) <- c(lvls,list(stnms))
return(Rho)
}
}
