predict.hmm.discnp <- function(object, y=NULL, tpm=NULL, Rho=NULL,
                               ispd=NULL, X=NULL, addIntercept=NULL,
                               warn=TRUE, drop=TRUE, ...) {
if(missing(object)) {
    if(is.null(y)) stop("No observations supplied.\n")
    if(is.null(tpm)) stop("Transition probabilities not supplied.\n")
    if(is.null(Rho)) stop("Emission probabilities not supplied.\n")
} else {
    if(is.null(tpm))          tpm          <- object$tpm
    if(is.null(Rho))          Rho          <- object$Rho
    if(is.null(ispd))         ispd         <- object$ispd
    if(is.null(X))            X            <- object$X
    if(is.null(addIntercept)) addIntercept <- object$addIntercept
    if(is.null(y)) {
       if(is.null(object$y)) {
           stop("No observations supplied.\n")
       }
       y    <- object$y
       numb <- object$numeric
    } else {
       y    <- tidyList(y)
       numb <- attr(y,"numeric")
    }
}
model      <- list(y=y,tpm=tpm,Rho=Rho,ispd=ispd,X=X,addIntercept=addIntercept)
stateProbs <- sp(model=model,warn=warn,drop=drop)
predictEngineHmmDiscnp(stateProbs,Rho,numb,drop=drop)
}
