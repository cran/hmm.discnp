ffun <- function(y,Rho)
{
#
# Function ffun.  To calculate f(y) = f(Y=y | the model parameters)
# for each value of the state k, k = 1, ..., K.  The returned result,
# fy, is a matrix whose ***rows*** correspond to the states and whose
# ***columns*** correspond to the observations y.
#

K <- ncol(Rho)
fy <- Rho[y,1:K]
fy[is.na(fy)] <- 1
t(fy)
}
