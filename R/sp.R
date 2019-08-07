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
    type <- switch(class(Rho),data.frame=1, matrix=2, list=3, array=4, NULL)
    if(is.null(type)) stop("Argument \"Rho\" is not of an appropriate form.\n")

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
