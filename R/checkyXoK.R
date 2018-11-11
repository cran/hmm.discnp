checkyXoK <- function(y,X) {
    lns <- sapply(y,nrow)
    nrX <- sapply(X,nrow)
    if(!isTRUE(all.equal(lns,nrX,check.attributes=FALSE)))
        stop(paste0("Mismatch between lengths of response vectors and\n",
                    "numbers of rows of predictor matrices.\n"))
    invisible()
}
