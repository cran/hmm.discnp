ffun <- local({

f1 <- function(Y,Rho) {
# Since f1 is being called, "newstyle" is TRUE so Rho is a data
# frame.  Note that Y is a data frame whose first column is y
# (a factor) and whose other columns are the predictors (possibly
# only the "Intercept", or "constant" predictor).
#
# Bozhemoi!!! Should the following *ever* be tested!!!
#
    y      <- Y$y
    Pred   <- Y[,-1,drop=FALSE]
    state  <- levels(Rho$state)
    rslt   <- vector("list",length(state))
    names(rslt) <- state
    indmat <- cbind(y,seq(along=y))
# Note that we are cbinding a *factor* to a numeric vector, which
# coerces the factor to numeric mode so that the resulting entries
# are integers, with i corresponding to the i-th level of the factor.
# That makes the indexing of "Tmp" by "indmat" in the following
# work correctly.
    for(k in state) {
        B   <- as.matrix(Rho[Rho$state==k,-(1:2)])
        Tmp <- B%*%t(Pred)
        Tmp <- t(t(Tmp) - apply(Tmp,2,max))
        Tmp <- exp(Tmp)
        den <- apply(Tmp,2,sum)
        num <- Tmp[indmat]
        prb <- num/den
        prb[is.na(prb)] <- 1
        rslt[[k]] <- prb
    }
    return(do.call(cbind,rslt))
}

f2 <- function(Y,Rho) {
# Since f2 is being called newstyle is FALSE here, Rho is a matrix.
# The argument Y is a data frame whose first column is a factor.
# Any other columns of Y are ignored here.
    rslt <- lapply(1:ncol(Rho),function(k,yf,Rho){
                                        v <- Rho[cbind(yf,k)]
                                        v[is.na(v)] <- 1
                                        return(v)
                                },yf=Y[,1],Rho=Rho)
    return(do.call(cbind,rslt))
}

f3 <- function(Y,Rho) {
# Here Y is a two-column data frame each column of which
# is a factor.
    f2(Y[,1,drop=FALSE],Rho[[1]]) * f2(Y[,2,drop=FALSE],Rho[[2]])
}

f4 <- function(Y,Rho) {
# Here Y is a two-column data frame each column of which is
# a factor.  We convert Y to a two-column matrix (of mode *character*).
Y <- as.matrix(Y)
lll <- apply(Y,1,is.na)
if(!any(lll)) {
    xxx <- lapply(1:nrow(Y),function(i,Rho,y){Rho[y[i,1],y[i,2],]},
                  Rho=Rho,y=Y)
} else {
    xxx <- lapply(1:nrow(Y),function(j,lll,Rho,y){
        u <- lll[,j]
        v <- y[j,]
        nat <- switch(1+sum(u),1,1+which(u),4)
        switch(nat,
            Rho[v[1],v[2],],
            apply(Rho[,v[2],,drop=FALSE],3,sum),
            apply(Rho[v[1],,,drop=FALSE],3,sum),
            rep(1,dim(Rho)[3])
        )
    },lll=lll,Rho=Rho,y=Y)
}
matrix(unlist(xxx),byrow=TRUE,ncol=dim(Rho)[3])
}

f  <- list(f1,f2,f3,f4)

function(Dat,Rho,type) {
#
# Function ffun to calculate
#        f(y) = Pr(Y=y | the model parameters and predictors)
# for each entry of each member of the list Dat, for each value of
# the state k, k = 1, ..., K.  The members of the list Dat may be
# factors (univariate, newstyle=FALSE), data frames whose columns
# are factors (bivariate, dependent) or matrices (bivariate dependent).
# In the bivariate setting, the rows of the data frames or matrices
# constitute the observations.  The returned result, fy, is a matrix
# whose ***rows*** correspond to the states and whose ***columns***
# correspond to observations.
#
# The "type" argument:
# type = 1 <--> univariate, newstyle = TRUE  -- f1
# type = 2 <--> univariate, newstyle = FALSE -- f2
# type = 3 <--> bivariate, independent       -- f3
# type = 4 <--> bivariate, dependent         -- f4

fy <- lapply(Dat,f[[type]],Rho=Rho)
if(any(sapply(fy,is.null))) {
   if(interactive) browser() else stop("Null fy component(s).\n")
}
fy <- do.call(rbind,fy)
t(fy)
}
})
