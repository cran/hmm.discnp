rhmm <- function(model,...,nsim,verbose=FALSE){
#
# Generic function rhmm (r for random) to simulate data from a
# hidden Markov model with discrete emissions the probabilities
# of which are specified non-parametrically (i.e. by means of a
# table or tables).
#
UseMethod("rhmm")
}
