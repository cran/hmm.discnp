nafracCalc <- function(y,drop=TRUE){
    if(inherits(y,"multipleHmmDataSets")) {
        rslt <- lapply(y,nafracCalc)
        if(drop & length(y)==1) rslt <- rslt[[1]]
        return(rslt)
    }
    if(!is.list(y)) y <- list(y)
    cls <- unique(sapply(y,class))
    if(length(cls) > 1) stop("The list \"y\" has entries of different classes.\n")
    if(!cls %in% c("matrix","character","integer","numeric","factor")) {
        whinge <- paste0("The class of the entries of the list \"y\" is\n",
                         "\"",cls,"\" which is inappropriate.\n")
        stop(whinge)
    }
    aaa <- sapply(y,function(x){apply(as.matrix(x),2,function(z){sum(is.na(z))})})
    bbb <- sapply(y,function(x){apply(as.matrix(x),2,length)},simplify="array")
    num <- if(is.matrix(aaa)) apply(aaa,1,sum) else (sum(aaa))
    den <- if(is.matrix(bbb)) apply(bbb,1,sum) else (sum(bbb))
    num/den
}
