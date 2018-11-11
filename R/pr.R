pr <- function(s,y,model=NULL,tpm,Rho,ispd=NULL,warn=TRUE) {
# Probability of state sequence(s).
# To do:  Something about supplying newstyle and X. (06/04/2017)

if(!is.null(model)) {
	tpm <- model$tpm
	Rho <- model$Rho
	ispd <- model$ispd
}

# Set the type:
type <- if(is.matrix(Rho)) 1 else if(is.list(Rho)) 2 else if(is.array(Rho)) 3
if(is.null(type)) stop("Argument \"Rho\" is not of an appropriate form.\n")

if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
if(is.null(ispd)) ispd <- revise.ispd(tpm)
y   <- tidyList(y)
Rho <- check.yval(y,Rho,type,warn=warn)

if(!is.list(s)) s <- list(s)
nseq <- length(s)
if(!length(y)%in%c(1,nseq))
	stop(paste("Mismatch between number of state sequences\n",
                   "and number of observation sequences.\n"))
rslt <- numeric(nseq)
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

fy  <- ffun(y, Rho,type) # Argument "type" was previously not
                         # supplied (06/02/2017).
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
