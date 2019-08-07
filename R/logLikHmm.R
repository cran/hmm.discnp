logLikHmm <- function(y,model=NULL,tpm=NULL,ispd=NULL,Rho=NULL,
                      X=NULL,addIntercept=NULL,warn=TRUE) {
#
# Function logLikHmm.  To calculate the log likelihood of a sequence,
# or collection (list) of sequences, of observations which come
# from a hidden Markov model with discrete non-parametric observation
# distributions.  These distributions are specified by:
#
# Univariate case.
#
#   * a matrix Rho where
#     P(Y = y_i | S = k) = Rho[i,k], OR
#
#   * a data frame with columns "y", "state", and further
#     columns of coefficients corresponding to the numerical
#     predictors given in X.  See the "Details" in the help
#     for hmm().
#
# Bivariate independent case.
#
#   * a pair of matrices Rho[[1]] and Rho[[2]] where
#     P(Y1 = y_i & Y2 = y_j | S = k) = Rho[[1]][i,k] * Rho[[2]][j,k]
#
# Bivariate dependent case.
#
#   * a 3-dimensional array Rho where
#     P(Y1 = y_i & Y2 = y_j | S = k) = Rho[i,j,k]
#
# In the foregoing Y, Y1, Y2 are observable random variables,
# y_i (resp. y_j) is the i-th (resp. j-th) possible value of such
# variables, and S is the hidden state.
#

# If y is not provided, just extract the log likelihood from "model".
if(missing(y)) {
   if(!is.null(model)) return(model$log.like)
   else stop("At least one of \"y\" or \"model\" must be supplied.\n")
}

# Get the parameters and addIntercept.
if(!is.null(model)) {
    Rho  <- model$Rho
    tpm  <- model$tpm
    ispd <- model$ispd
    addIntercept <- model$args$addIntercept
}

if(is.null(Rho)) stop("\"Rho\" not supplied.\n")

# Set the type:
type <- if(is.data.frame(Rho)) 1
        else if(is.matrix(Rho)) 2
        else if(is.list(Rho) & !is.data.frame(Rho)) 3
        else if(is.array(Rho)) 4
if(is.null(type)) stop("\"Rho\" is not of an appropriate form.\n")

if(inherits(y,"madeDat")) {
    Dat <- y
} else {
    y <- tidyList(y,rp="response")

# Make sure that the entries of the vectors in y correspond
# to the appropriate dimension names of Rho.
    Rho <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)

# If we are using predictors, tidy them up.
if(type==1 & !is.null(X)) {
    X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
    checkyXoK(y,X)
}

# Form the data list.
Dat <- makeDat(y,X)
}

# If K=1 do the triv thing:
K <- switch(EXPR=type,length(levels(Rho$state)),ncol(Rho),
            ncol(Rho[[1]]),dim(Rho)[3])
if(K==1) return(sum(log(ffun(Dat,Rho,type))))

# K is not equal to 1; need tpm and ispd.
if(is.null(tpm)) stop("Transition probability matrix not supplied.\n")
if(is.null(ispd)) {
    ispd <- revise.ispd(tpm=tpm)
}

lns <- sapply(Dat,nrow)
fy  <- ffun(Dat,Rho,type)
rp  <- try(recurse(fy,tpm,ispd,lns))
if(inherits(rp,"try-error")) {
    if(interactive()) browser() else stop("Problem with recurse().\n")
}
sum(log(rp$llc))
}
