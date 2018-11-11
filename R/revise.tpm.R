revise.tpm <- function(xi,mixture) {
    if(mixture)  {
        matrix(apply(xi,2,sum)/sum(xi),byrow=TRUE,
               nrow=nrow(xi),ncol=ncol(xi))
    } else {
        den <- apply(xi,1,sum)
        if(any(den == 0))
            stop("At least one row of new tpm would be all zeroes.\n")
        xi/den
    }
}
