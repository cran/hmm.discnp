pr <- function(s,y,model=NULL,tpm=NULL,Rho=NULL,ispd=NULL,warn=TRUE) {
# Probability of state sequence(s).
# To do:  Something about supplying X.

if(!missing(y) && inherits(y,"multipleHmmDataSets")) {
    if(!inherits(s,"list") | length(s) != length(y)) {
        whinge <- paste0("Argument \"y\" is of class \"multipleHmmDataSet\"\n",
                         "  and argument \"s\" is of an incompatible",
                         " form or length.\n")
        stop(whinge)
    }
    ndat <- length(s)
    rslt <- lapply(1:ndat,function(k,s,y,tpm,Rho,ispd,warn){
                              sk <- s[[k]]
                              yk <- y[[k]]
                              pr(sk,yk,tpm=tpm,Rho=Rho,ispd=ispd,warn=warn)
                          },s=s,y=y,tpm=tpm,Rho=Rho,ispd=ispd,warn=warn)
    return(rslt)
}

if(!is.null(model)) {
	tpm  <- model$tpm
	Rho  <- model$Rho
	ispd <- model$ispd
}
if(is.null(tpm) | is.null(Rho))
    stop("At least one of \"tpm\" and \"Rho\" was not supplied.\n")

# Convert Rho if necessary (just for consistency with the
# overall pattern; we'll actually convert it *back* again later).
if(inherits(Rho,"matrix")) Rho <- cnvrtRho(Rho)

# Set the type:
if(inherits(Rho,"data.frame")) {
    type <- 1
} else if(inherits(Rho,"list")) {
    type <- 2
} else if(inherits(Rho,"array")) {
    type <- 3
} else {
    stop("Object \"Rho\" has an incorrect class.\n")
}

if(missing(y)) {
	y <- if(!is.null(model)) model[["y"]] else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
y   <- tidyList(y)
y   <- makeDat(y,X=NULL)

if(is.null(ispd)) ispd <- revise.ispd(tpm)

if(!inherits(s,"list")) s <- list(s)
nseq <- length(s)
if(!length(y)%in%c(1,nseq))
	stop(paste("Mismatch between number of state sequences\n",
                   "and number of observation sequences.\n"))
rslt <- numeric(nseq)
Rho  <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)

# Set the names of the states.
stnms <- switch(type,
             levels(Rho$state),
             colnames(Rho[[1]]),
             dimnames(Rho)[[3]]
         )
if(is.null(stnms)) stnms <- as.character(1:nrow(tpm))
names(ispd) <- stnms
rownames(tpm) <- colnames(tpm) <- stnms

# Get the log likelihood vector; needed later.
fy  <- ffun(y,Rho,type) # To do: Need to allow predictors X here.
lns <- sapply(y,nrow)
rp  <- recurse(fy, tpm, ispd, lns)
ll  <- log(rp$llc)

# Possibly convert Rho back to matrix form.
if(type==1) Rho <- cnvrtRho(Rho)

# Rubber hits road.
for(i in 1:nseq) {
	ii <- if(length(y)==1) 1 else i
        si <- s[[i]]
        n  <- length(si)
        y1 <- y[[ii]][,1]
        if(n != length(y1))
            stop(paste("Mismatch between length of state sequence\n",
                       "and length of observation sequence.\n",sep=""))
        if(type > 1) y2 <- y[[ii]][,2]
        bit1 <- log(ispd[si[1]]) + sum(log(tpm[cbind(si[-n],si[-1])]))
        bit2 <- switch(type,
                       sum(log(Rho[cbind(y1,si)])),
                       sum(c(log(Rho[[1]][cbind(y1,si)]),log(Rho[[2]][cbind(y2,si)]))),
                       sum(log(Rho[y1,y2,si]))
                 )
	rslt[i] <- bit1 + bit2
}


if(nseq==1) return(exp(rslt-sum(ll))) # sum(ll)???

part2 <- numeric(nseq)
jstop <- 0
for(i in 1:nseq) {
        jstart <- jstop+1
        jstop  <- jstop + lns[i]
        part2[i] <- sum(ll[jstart:jstop])
}
exp(rslt-part2)
}
