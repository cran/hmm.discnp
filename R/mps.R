mps <- function(y,model=NULL,tpm=NULL,Rho=NULL,ispd=NULL,warn=TRUE){
#
# Function mps: most probable states.
#
# To do:  Something about supplying X. (10/01/2020)
#

if(inherits(y,"multipleHmmDataSets")) {
    rslt <- lapply(y,function(x,model,tpm,Rho,ispd,warn){
                              mps(x,model,tpm,Rho,ispd,warn)
                     },model=model,tpm=tpm,Rho=Rho,
                       ispd=ispd,warn=warn)
    return(rslt)
}

if(!is.null(model)) {
	tpm  <- model$tpm
	Rho  <- model$Rho
	ispd <- model$ispd
}
if(is.null(tpm) | is.null(Rho))
    stop("At least one of \"tpm\" and \"Rho\" was not supplied.\n")


# Convert Rho if necessary.
if(inherits(Rho,"matrix")) Rho <- cnvrtRho(Rho)

stnms <- if(is.null(rownames(tpm))) as.character(1:nrow(tpm)) else rownames(tpm)
if(is.null(ispd)) ispd <- revise.ispd(tpm)
if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
y <- tidyList(y)
y <- makeDat(y,X=NULL)

# Set the type:
if(inherits(Rho,"data.frame")) {
    type <- 1
} else if(inherits(Rho,"list")) {
    type <- 2
} else if(inherits(Rho,c("array"))) {
    type <- 3
} else {
    stop("Object \"Rho\" has an incorrect class.\n")
}

Rho  <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)
lns  <- sapply(y,nrow)
nseq <- length(y)
fy   <- ffun(y,Rho,type)
rp   <- recurse(fy, tpm, ispd, lns)
xxx  <- apply(rp$gamma, 2, which.max)
if(!is.null(stnms)) xxx <- stnms[xxx]
if(nseq==1) return(xxx)
rslt  <- list()
jstop <- 0
for(i in 1:nseq) {
        jstart <- jstop+1
        jstop  <- jstop + lns[i]
	rslt[[i]] <- xxx[jstart:jstop]
}
rslt
}
