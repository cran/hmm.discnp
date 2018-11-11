misstify <- local({

fix1 <- function(x,nafrac,p2){
if(ncol(x)==1) return(x)
# p2 = prob. of both entries *non*-missing in bivariate case.
    pq <- nafrac/sum(nafrac)
    x1 <- x[1,]
    j  <- sample(1:2,1,prob=pq)
    ind <- rbinom(1,1,p2)
    if(ind==0) x1[j] <- NA
    x[1,] <- x1
    x
}

fix2 <- function(x,nafrac,present){
    pad <- if(present) 0 else NULL
    s   <- if(present) 1 else 0
    ina1 <- c(pad,rbinom(nrow(x)-s,1,nafrac[1]))
    x[ina1 == 1,1] <- NA
    if(ncol(x)==2) {
        ina2 <- c(pad,rbinom(nrow(x)-s,1,nafrac[2]))
        x[ina2 == 1,2] <- NA
    }
    x
}

function(X,nafrac,fep=NULL) {
    if(any(nafrac >= 1))
        stop("Makes no sense to have a missing-fraction equal to 1.\n")
    if(is.null(fep)) {
        fep <- list(present=TRUE)
    }
    bivar <- ncol(X[[1]]) == 2
    if(bivar) {
        if(length(nafrac) == 1) nafrac <- rep(nafrac,2)
        if(length(fep) == 1) {
            fep <- c(fep,list(p2=prod(1-nafrac)/(1 - prod(nafrac))))
        } else if(fep[[2]] < 0 | fep[[2]] > 1) {
            whinge <- paste0("Component \"p2\" of \"fep\" must be a probability.\n")
            stop(whinge)
        }
    }
    X <- lapply(X,fix2,nafrac=nafrac,present=fep[[1]])
    if(fep[[1]] & bivar && fep[[2]] < 1) {
        X <- lapply(X,fix1,nafrac=nafrac,p2=fep[[2]])
    }
    X
}})
