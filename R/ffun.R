ffun <- function(y,Rho)
{
#
# Function ffun to calculate f(x) = Pr(Y=x | the model parameters)
# for each entry of each vector in the list y, for each value of the
# state k, k = 1, ..., K.  The returned result, fy, is a matrix whose
# ***rows*** correspond to the states and whose ***columns***
# correspond to the observations y.
#

fy <- lapply(y,function(x,Rho){Rho[x,1:ncol(Rho)]},Rho=Rho)
fy <- do.call(rbind,fy)
fy[is.na(fy)] <- 1
t(fy)
}
