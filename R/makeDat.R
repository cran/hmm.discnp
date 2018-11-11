makeDat <- function(y,X,addIntercept) {
#
# Combines the list of response objects y and the list of predictor
# objects X into a single list of objects whose entries form the
# "complete" data, i.e. encompass both response and predictor.
#
# If the data are unvariate, then each entry of the result returned
# by this function is a data frame, the first column of which is
# a *factor* constituting the (univariate) observations and whose
# other columns consist of corresponding predictors.  (There may
# be only one such "predictor", the constant term or intercept,
# in which case the predictor is a column of 1's.)
#
# If the data are bivariate then each entry of the result returned
# by this function is a data frame whose first column is a *factor*
# constituting the first of the (bivariate) observations and
# whose second column is a *factor* constituting the second of the
# (bivariate) observations.  The data consist of responses only
# (no predictors).
#

# Make sure that X, if present, is of class "tidyList".
    if(!is.null(X)) {
        X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
    }

# Get the levels of the "y" factor (the putatively possible values
# of the emissions.
    lvls <- attr(y,"lvls")

# Do the right thing, according to parity.
    parity <- attr(y,"parity")
    if(parity=="univar") {
        n    <- length(y)
        rslt <- vector("list",n)
        for(i in seq(along=rslt)) {
            yi <- data.frame(y=factor(y[[i]],levels=lvls))
            Xi <- if(is.null(X)) data.frame(Intercept=1) else X[[i]]
            rslt[[i]] <- cbind(yi,Xi)
        }
        attr(rslt,"prednames") <- if(is.null(X)) "Intercept" else attr(X,"prednames")
        attr(rslt,"lvls") <- lvls
        class(rslt) <- c("madeDat","list")
        return(rslt)
    }
    if(parity=="bivar") {
        rslt <- lapply(y,function(x,lvls){xdf <- as.data.frame(x)
                               xdf[,1] <- factor(xdf[,1],levels=lvls[[1]])
                               xdf[,2] <- factor(xdf[,2],levels=lvls[[2]])
                               return(xdf)
                              },lvls=lvls)
        attr(rslt,"lvls") <- lvls
        class(rslt) <- c("madeDat","list")
        return(rslt)
    }
    stop("Something must be wrong with the parity attribute of \"y\".\n")
} 
