fitted.hmm.discnp <- function(object,...) {
y <- object$y
if(is.null(y)) stop("Observations \"y\" were not kept.\n")
if(!is.numeric(y[[1]]))
	stop("Observations are not numeric; fitted values make no sense.\n")
sp(y,object,means=TRUE)$means
}
