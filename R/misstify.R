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

function(y,nafrac,fep=NULL) {
    if(inherits(y,"multipleHmmDataSets")) {
        bivar <- ncol(y[[1]][[1]]) == 2
        ny <- length(y)
        if(ny==0)
            stop("Something is wrong. Argument \"y\" has zero length.\n")
        nafracv <- rep(nafrac,length=if(bivar) 2*ny else ny)
        nafracl <- lapply(1:ny,function(k,nafracv,bivar){
                                   if(bivar) nafracv[2*(k-1) + 1:2] else nafracv[k]
                               },nafracv=nafracv,bivar=bivar)
        xxx <- lapply(1:ny,function(k,y,nafrac,fep){
                                        misstify(y[[k]],nafracl[[k]],fep)
                                    },y=y,nafrac=nafrac,fep=fep)
        class(xxx) <- class(y)
        return(xxx)
    }
    nafrac <- rep(nafrac,length=ncol(y[[1]]))
    if(!all(nafrac >= 0 & nafrac < 1)) {
        whinge <- paste0("All entries of \"nafrac\" must be probabilities",
                         " strictly less than 1.\n")
        stop(whinge)
     }
    if(is.null(fep)) {
        fep <- list(present=TRUE)
    }
    bivar <- length(nafrac) == 2
    if(bivar) {
        if(length(fep) == 1) {
            fep <- c(fep,list(p2=prod(1-nafrac)/(1 - prod(nafrac))))
        } else if(fep[[2]] < 0 | fep[[2]] > 1) {
            whinge <- paste0("Component \"p2\" of \"fep\" must be a probability.\n")
            stop(whinge)
        }
    }
    y <- lapply(y,fix2,nafrac=nafrac,present=fep[[1]])
    if(fep[[1]] & bivar && fep[[2]] < 1) {
        y <- lapply(y,fix1,nafrac=nafrac,p2=fep[[2]])
    }
    y
}})
