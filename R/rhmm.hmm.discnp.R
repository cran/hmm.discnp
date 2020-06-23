rhmm.hmm.discnp <- function(model,...,nsim=1,verbose=FALSE,
                            inMiss=TRUE,fep=NULL,drop=TRUE,forceNumeric=TRUE) {
#
# Function rhmm.hmm to simulate data from a fitted hidden Markov
# model with discrete emissions, the probabilities of which are
# specified non-parametrically (i.e. by means of a table or tables).
#

if(!inherits(model,"hmm.discnp"))
    stop("Argument \"model\" is not of class \"hmm.discnp\".\n")

tpm  <- model$tpm
Rho  <- model$Rho
ispd <- model$ispd
yval <- if(is.null(model$y)) NULL else attr(model$y,"uval")
nafrac <- if(inMiss) model$nafrac else NULL

ylengths <- model$ylengths
rhmm.default(ylengths=ylengths,nsim=nsim,verbose=verbose,
             nafrac=nafrac,fep=fep,tpm=tpm,Rho=Rho,ispd=ispd,
             yval=yval,drop=drop,forceNumeric=forceNumeric)
}
