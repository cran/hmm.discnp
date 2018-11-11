nafracCalc <- function(y){
    if(is.null(y)) return(0)
    if(!is.list(y)) y <- list(y)
    aaa <- sapply(y,function(x){apply(as.matrix(x),2,function(z){sum(is.na(z))})})
    bbb <- sapply(y,function(x){apply(as.matrix(x),2,length)},simplify="array")
    num <- if(is.matrix(aaa)) apply(aaa,1,sum) else (sum(aaa))
    den <- if(is.matrix(bbb)) apply(bbb,1,sum) else (sum(bbb))
    num/den
}
