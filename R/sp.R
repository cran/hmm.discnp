sp <- function (y, model = NULL, tpm=NULL, Rho=NULL, ispd=NULL, X=NULL,
                addIntercept=NULL,means=FALSE, warn=TRUE) {
# State probabilities.
#
# To do:  Something about supplying newstyle and X. (06/04/2017)
    if (!is.null(model)) {
        tpm <- model$tpm
        Rho <- model$Rho
        ispd <- model$ispd
        addIntercept <- model$args$addIntercept
    }
    if(is.null(Rho)) stop("\"Rho\" not supplied.\n")
    if(is.null(tpm)) stop("Transition probability matrix not supplied.\n")
# Set the type:
    type <- if(is.data.frame(Rho)) 1
            else if(is.matrix(Rho)) 2
            else if(is.list(Rho) & !is.data.frame(Rho)) 3
            else if(is.array(Rho)) 4
    if(is.null(type)) stop("Argument \"Rho\" is not of an appropriate form.\n")
    newstyle <- type==1

    if(is.null(ispd)) ispd <- revise.ispd(tpm)
    if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
    }
    y    <- tidyList(y)
    Rho  <- check.yval(y,Rho,type,warn=warn)

# If we are using predictors, tidy them up.
    if(type==1 & !is.null(X)) {
        X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
        checkyXoK(y,X)
    }

# Form the data list.
    Dat  <- makeDat(y,X)

# Calculate the gammas.
    lns  <- sapply(y,nrow)
    fy   <- ffun(Dat, Rho,type)
    rp   <- recurse(fy, tpm, ispd, lns)
    prbs <- rp$gamma

    if (means) {
        switch(type,
# Univariate, Rho a data frame.
        {
            yval <- as.numeric(levels(Rho$y))
            if (any(is.na(yval))) 
                stop("Non-numeric y-values; means make no sense.\n")
            Roe <- cnvrtRho(Rho)
            cmns <- apply(yval * Roe, 2, sum)
            mns  <- apply(cmns * prbs, 2, sum)
        },

# Univariate, Rho a matrix.
        {
            yval <- as.numeric(row.names(Rho))
            if (any(is.na(yval))) 
                stop("Non-numeric y-values; means make no sense.\n")
            cmns <- apply(yval * Rho, 2, sum)
            mns  <- apply(cmns * prbs, 2, sum)
        },

# Bivariate independent:
        {
            mns <- vector("list", 2)
            for (j in 1:2) {
                yval <- as.numeric(row.names(Rho[[j]]))
                if (any(is.na(yval))) 
                  stop("Non-numeric y-values; means make no sense.\n")
                cmns <- apply(yval * Rho[[j]], 2, sum)
                mns[[j]] <- apply(cmns * prbs, 2, sum)
            }
        },
# Bivariate dependent:
        {
            yval <- vector("list", 2)
            yval[[1]] <- as.numeric(dimnames(Rho)[[1]])
            if (any(is.na(yval[[1]]))) 
                stop("Non-numeric y-values in first variable; means make no sense.\n")
            yval[[2]] <- as.numeric(dimnames(Rho)[[2]])
            if (any(is.na(yval[[2]]))) 
                stop("Non-numeric y-values in second variable; means make no sense.\n")
            mns <- vector("list", 2)
            for (j in 1:2) {
                RT <- apply(Rho, c(j, 3), sum)
                cmns <- apply(yval * RT, 2, sum)
                mns[[j]] <- apply(cmns * prbs, 2, sum)
            }
       })
    }
    nseq <- length(lns)
    if (nseq == 1) {
	if(means) return(list(probs=prbs,means=mns))
	return(prbs)
    }
    xxx <- vector("list",nseq)
    if(means) yyy <- vector("list",nseq)
    istop <- 0
    for(i in 1:nseq) {
        istart <- istop+1
        istop  <- istop + lns[i]
        xxx[[i]] <- prbs[,istart:istop]
        if(means) {
            if(is.list(mns)) {
                yyy[[i]] <- list(mns[[1]][istart:istop],mns[[2]][istart:istop])
            } else yyy[[i]] <- mns[istart:istop]
        }
    }
if(means) return(list(probs=xxx,means=yyy))
xxx
}
