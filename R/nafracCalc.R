nafracCalc <- function(y,drop=TRUE){
    if(inherits(y,"multipleHmmDataSets")) {
        rslt <- lapply(y,nafracCalc)
        if(drop & length(y)==1) rslt <- rslt[[1]]
        return(rslt)
    }
    if(!is.list(y)) y <- list(y)

# Check that the entries of y all have the same class and that this
# class is appropriate.
    elsie <- unique(sapply(y,function(x){length(class(x))}))
    if(length(elsie) > 1)
        stop("The list \"y\" has entries with classes of different lengths.\n")
    m   <- matrix(unlist(lapply(y,class)),nrow=length(y),byrow=TRUE)
    ndm <- !duplicated(m)
    if(sum(ndm) > 1)
        stop("The list \"y\" has entries of different classes.\n")
    if(!inherits(y[[1]],
            c("matrix","character","integer","numeric","factor"))) {
        cls <- paste(m[ndm,],collapse=" ")
        whinge <- paste0("The class of the entries of the list \"y\" is\n",
                         "  \"",cls,"\" which is inappropriate.\n")
        stop(whinge)
    }
 
# Now count the missing values.
    aaa <- sapply(y,function(x){apply(as.matrix(x),2,function(z){sum(is.na(z))})})
    bbb <- sapply(y,function(x){apply(as.matrix(x),2,length)},simplify="array")
    num <- if(is.matrix(aaa)) apply(aaa,1,sum) else (sum(aaa))
    den <- if(is.matrix(bbb)) apply(bbb,1,sum) else (sum(bbb))
    num/den
}
