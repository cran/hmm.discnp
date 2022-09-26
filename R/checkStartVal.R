checkStartVal <- function(par0,K,indep,yval,rand.start,
                          mixture,prednames) {
#
bivar <- inherits(yval,"list")
if(is.null(par0)) {
    if(bivar & is.null(indep)) {
        stop("Neither \"indep\" nor \"par0$Rho\" have been supplied.\n")
    }
    if(is.null(K)) {
        stop("One of \"par0\" and \"K\" must be specified.\n")
    } else if(K==1) {
        par0 <- NA
    } else {
        par0 <- init.all(K,rand.start,mixture,indep,yval,prednames)
    }
} else {

# Check on par0$tpm.
    tpm <- par0$tpm
    if(!inherits(tpm,"matrix") || nrow(tpm) != ncol(tpm)) {
        whinge <- paste0("Argument \"par0\" if supplied must",
                         " have a component\n",
                         "  \"tpm\" which is a square matrix.\n")
        stop(whinge)
    }
    if(is.null(K)) {
        K <- nrow(tpm)
    } else {
        if(K != nrow(tpm)) {
            stop("The values of \"K\" and \"par0$tpm\" are inconsistent.\n")
        }
    }
# Check on par0$Rho.
    if(bivar) {
        if(inherits(par0$Rho,"data.frame")) {
            stop("In the bivariate setting \"par0$Rho\" cannot be a data frame.\n")
        } else if(inherits(par0$Rho,"list")) {
            if(is.null(indep)) {
                indep <- TRUE
            } else if(!indep) {
                 stop(paste("The value of \"indep\" and that of \"par0$Rho\"\n",
                            "are inconsistent.\n",sep=""))
            }
        } else {
            if(is.null(indep)) {
                indep <- FALSE
            } else if(indep) {
                 stop(paste("The value of \"indep\" and that of \"par0$Rho\"\n",
                            "are inconsistent.\n",sep=""))
            }
        }
        if(indep) {
            for(i in 1:2) {
                if(length(yval[[i]]) == nrow(par0$Rho[[i]])) {
                    rownames(par0$Rho[[i]]) <- yval[[i]]
                } else {
                    ordnl <- if(i==1) "first" else "second"
                    whinge <- paste0("The length of the ",ordnl," component ",
                                     " of \"yval\" is incompatible\n",
                                     " with the row dimension of the ",ordnl,
                                     " component of par0$Rho.\n")
                    stop(whinge)
                }
            }
            if(!identical(colnames(par0$Rho[[1]]),colnames(par0$Rho[[2]]))) {
                whinge <- paste0("Inconsistency in the column names of\n",
                                 "  the components of par0$Rho.\n")
                stop(whinge)
            }
        } else {

# Do the dependent case.
             yn <- sapply(yval,length)
             if(!isTRUE(all.equal(yn,dim(par0$Rho)[1:2]))) {
                stop("Lengths of \"yval\" and dim(par0$Rho) are incompatible.\n")
             }
             dimnames(par0$Rho)[1:2] <- yval
        }
    } else {

# Convert par0$Rho to a data frame if necessary.
        if(inherits(par0$Rho,"matrix")) {
            if(is.null(rownames(par0$Rho))) {
                if(length(yval) != nrow(par0$Rho)) {
                    whinge <- paste0("The specified \"Rho\" has no row",
                                     " names and\n","  the specified ",
                                     "\"yval\" is of the wrong length\n",
                                     "  to form such row names.\n")
                    stop(whinge)
                }
                rownames(par0$Rho) <- yval
            }
            newRho <- cnvrtRho(par0$Rho)
            ncx <- length(prednames)
            if(prednames[[1]]=="Intercept") {
                Cx <- matrix(0,nrow(newRho),ncx-1)
                colnames(Cx) <- prednames[-1]
                newRho <- cbind(newRho,as.data.frame(Cx))
            } else {
                Cx <- matrix(0,nrow(newRho),ncx)
                colnames(Cx) <- prednames
                newRho <- cbind(newRho[,-3],as.data.frame(Cx))
            }
            par0$Rho <- newRho
        }
        lvls <- levels(par0$Rho$y)
        if(!isTRUE(all.equal(sort(yval),sort(lvls)))) {
            whinge <- paste0("The levels of \"y\" in the specified \"Rho\"\n",
                             "  and \"yval\" are incompatible.\n")
            stop(whinge)
        }
    }
}

# Get the state names set properly; take these to be the
# common row and column names of tpm.
if(identical(par0,NA)) {
    stnms <- NULL
} else {
    if(inherits(par0$Rho,"data.frame")) {
        type <- 1
    } else if(inherits(par0$Rho,"list")) {
        type <- 2
    } else if(inherits(par0$Rho,"array")) {
        type <- 3
    } else {
        stop("This is impossible, but par0$Rho is of the wrong class.\n")
    }
    stnms <- rownames(par0$tpm)
    if(!identical(stnms,colnames(par0$tpm)))
        stop("Inconsistency in row and column names of par0$tpm.\n")
    if(is.null(stnms)) {
        stnms <- switch(EXPR=type,
                     levels(par0$Rho$state),
                     colnames(par0$Rho[[1]]),
                     dimnames(par0$Rho)[[3]]
                 )
    }
    if(is.null(stnms)) stnms <- as.character(1:K)
}

list(par0=par0,K=K,yval=yval,indep=indep,stnms=stnms)
}
