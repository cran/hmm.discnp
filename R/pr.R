pr <- function(s,y,model=NULL,tpm,Rho,ispd=NULL,warn=TRUE) {
# Probability of state sequence(s).
# To do:  Something about supplying X.

if(!is.null(model)) {
	tpm  <- model$tpm
	Rho  <- model$Rho
	ispd <- model$ispd
}

# If Rho is presented in the "logistic style" parameterisation,
# convert it to a matrix of probabilities.
if(class(Rho)=="data.frame") Rho <- cnvrtRho(Rho)

# Set the type:
type <- switch(class(Rho),data.frame=1, matrix= 2, list=3, array=4, NULL)
if(is.null(type)) stop("Argument \"Rho\" is not of an appropriate form.\n")
# Note: type = 1 can't happen here.

if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
if(is.null(ispd)) ispd <- revise.ispd(tpm)
y   <- tidyList(y)
y   <- makeDat(y,X=NULL)

if(!is.list(s)) s <- list(s)
nseq <- length(s)
if(!length(y)%in%c(1,nseq))
	stop(paste("Mismatch between number of state sequences\n",
                   "and number of observation sequences.\n"))
rslt <- numeric(nseq)
Rho  <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)
for(i in 1:nseq) {
	ii <- if(length(y)==1) 1 else i
        si <- s[[i]]
        n  <- length(si)
        y1 <- y[[ii]][,1]
        if(n != length(y1))
            stop(paste("Mismatch between length of state sequence\n",
                       "and length of observation sequence.\n",sep=""))
        if(type > 2) y2 <- y[[ii]][,2]
        bit1 <- log(ispd[si[1]]) + sum(log(tpm[cbind(si[-n],si[-1])]))
        bit2 <- switch(type,
                       NA,
                       sum(log(Rho[cbind(y1,si)])),
                       sum(c(log(Rho[[1]][cbind(y1,si)]),log(Rho[[2]][cbind(y2,si)]))),
                       sum(log(Rho[y1,y2,si]))
                 )
	rslt[i] <- bit1 + bit2
}

fy <- ffun(y,Rho,type)
# To do: Need to take account of predictors X here.
lns <- sapply(y,nrow)
rp  <- recurse(fy, tpm, ispd, lns)
ll  <- log(rp$llc)

if(nseq==1) return(exp(rslt-ll))

jstop <- 0
for(i in 1:nseq) {
        jstart <- jstop+1
        jstop  <- jstop + lns[i]
	rslt[i] <- exp(rslt[i]-sum(ll[jstart:jstop]))
}
rslt
}
