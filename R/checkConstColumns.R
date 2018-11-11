checkConstColumns <- function(y,prednames) {
ny    <- length(y)
eps   <- sqrt(.Machine$double.eps)
noCC  <- TRUE
oneOK <- TRUE
for(i in ny) {
    xxx <- apply(y[[i]],2,function(x){all(diff(x)==0)})
    if(sum(xxx)>1) {
        j <- which(xxx)
        whinge <- paste0("Predictor matrix ",i," has at least two\n",
                         "  constant columns.\n")
        stop(whinge)
    }
    xxx <- apply(y[[i]],2,function(x){all(abs(x-1) <= eps)})
    if(sum(xxx)) noCC <- FALSE
    if(!xxx[1]) oneOK <- FALSE
    if(!(noCC | oneOK)) {
        whinge <- paste0("Predictor matrix ",i," has a column of 1's\n",
                         "  which is not the first column.\n")
        stop(whinge)
    }
}
if(!(noCC | oneOK)) {
    whinge("Some predictor matrices have a first column of 1's\n",
           " and some don't.\n")
    stop(whinge)
}

if(oneOK) {
    if(!is.null(prednames) && prednames[1] != "Intercept")
        stop("The initial column of 1's must be named \"Intercept\".\n")
}
oneOK
}
