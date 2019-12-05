sp <- function (y, model = NULL, tpm=NULL, Rho=NULL, ispd=NULL, X=NULL,
                addIntercept=NULL, warn=TRUE, drop=TRUE) {
#
# State probabilities.
#
    if (!is.null(model)) {
        if(missing(y)) y <- model[["y"]]
        tpm <- model$tpm
        Rho <- model$Rho
        ispd <- model$ispd
        addIntercept <- model$args$addIntercept
        X <- model[["X"]]
    }
    if(missing(y) | is.null(y))
            stop("No observations supplied.\n")
    if(is.null(Rho)) stop("\"Rho\" not supplied.\n")
    if(is.null(tpm)) stop("Transition probability matrix not supplied.\n")
    if(is.null(ispd)) ispd <- revise.ispd(tpm)

# Set the type:
if(inherits(Rho,"data.frame")) {
    type <- 1
} else if(inherits(Rho,"list")) {
    type <- 3
} else if(inherits(Rho,c("matrix","array"))) {
    if(length(dim(Rho))==2) type <- 2
    else if(length(dim(Rho))==3) type <- 4
    else stop("Object \"Rho\" can be of dimension 2 or 3 only.\n")
} else {
    stop("Object \"Rho\" has an incorrect class.\n")
}

# Tidy up y and check on compatibility of y and Rho.
    y <- tidyList(y)
    Rho  <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)

# If we are using predictors, tidy them up.
    if(is.data.frame(Rho)) {
        if(ncol(Rho) > 3) {
            if(is.null(X))
                stop("Predictors \"X\" are needed and were not supplied.\n")
            X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
            checkyXoK(y,X)
        }
    }

# Form the data list.
    Dat  <- makeDat(y,X)

# Calculate the gammas.
    lns  <- sapply(y,nrow)
    fy   <- ffun(Dat,Rho,type)
    rp   <- recurse(fy, tpm, ispd, lns)
    prbs <- rp$gamma

    nseq <- length(lns)
    if (nseq == 1) {
	if(drop) return(prbs) else return(list(prbs))
    }
    xxx <- vector("list",nseq)
    istop <- 0
    for(i in 1:nseq) {
        istart <- istop+1
        istop  <- istop + lns[i]
        xxx[[i]] <- prbs[,istart:istop]
    }
    xxx
}
