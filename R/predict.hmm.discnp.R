predict.hmm.discnp <- function(object,y=NULL,...) {
if(is.null(y)) y <- object$y
if(is.null(y)) stop("Observations \"y\" were not kept.\n")
if(!is.list(y)) y <- list(y)
if(!all(sapply(y,is.numeric)))
	stop(paste("Some values of \"y\" are not numeric;\n",
                   "predicted values make no sense.\n"))
sp(y,object,means=TRUE)$means
}
