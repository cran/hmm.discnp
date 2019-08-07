fitted.hmm.discnp <- function(object,warn=TRUE,drop=TRUE,...) {
y <- object[["y"]]
if(is.null(y)) stop("Observations \"y\" were not kept.\n")

# Check on numeracy.  If the observations are numeric, then
# the fitted values are the expected values.  If the observations
# are categorical then the "fitted values" consist of ... ???
numb <- object$numeric

# Get the state probabilities.
stateProbs <- sp(model=object,warn=warn)

# Dig out Rho.
Rho <- object$Rho

predictEngineHmmDiscnp(stateProbs,Rho,numb,drop)
}
