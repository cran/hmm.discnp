tidyup <- function(y,yval=NULL) {
if(is.null(yval)) yval <- sort(unique(as.vector(y)))
nval <- length(yval)
if(is.matrix(y)) {
	nr <- nrow(y)
	nc <- ncol(y)
}
else {
	nr <- length(y)
	nc <- 1
}
matrix(match(y,yval),nr,nc)
}
