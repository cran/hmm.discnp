fitted.hmm.discnp <- function(object,...) {
y <- object[["y"]]
if(is.null(y)) stop("Observations \"y\" were not kept.\n")

# Check on numeracy.
num <- object$numeric
if(!num)
	stop(paste("Original observations were not numeric;\n",
                   "fitted values make no sense.\n"))
if(!is.list(y)) y <- list(y)
# Do *not* need to convert the y values (which are character)
# to numeric; it's the capacity of the dimnames of Rho to be
# interpreted as numeric that determines the validity of the
# calculations.

# Do the calculations.
Rho <- object$Rho
if(is.data.frame(Rho)) {
    if(ncol(Rho) > 3) {
        X <- object$X
        if(is.null(X)) stop("Predictors \"X\" are needed and were not kept.\n")
    } else X <- NULL
} else X <- NA
sp(y,object,X=X,means=TRUE)$means
}
