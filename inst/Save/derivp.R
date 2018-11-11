derivp <- function(Zeta) {
# Zeta is an m x K matrix of parameters zeta_ij determining
# a corresponding matrix of probabilities P = [p_ij] such
# that either (1) the *rows* sum to 1 or (2) the *columns* sum to 1.
#
# In the case (1) the p_ij are expressed in terms of the
# zeta_ij as p_ij = exp(zeta_ij)/(1 + Delta_i) where
# Delta_i = exp(zeta_i1) + exp(zeta_i2) + ... + exp(zeta_{i,K-1}.
# The last entry of each row of Zeta (i.e. zeta_{i,K} for the
# i-th row) must be zero.
#
# In the case (2) the p_ij are expressed in terms of the
# zeta_ij as p_ij = exp(zeta_ij)/(1 + Delta_j) where
# Delta_j = exp(zeta_1j) + exp(zeta_2j) + ... + exp(zeta_{m-1,j}.
# The last entry of each column of Zeta (i.e. zeta_{m,j} for the
# j-th column) must be zero.
# 
# The arrays returned are of (respectively) the first and second derivatives
# of the p_ij with respect to phi_k, k = 1, ..., npar = m*(K-1) in
# the first case and (m-1)*K in the second case.  The phi_k are the
# entries of Zeta "strung" out column by column (with redundancies
# removed).  I.e. in the first case phi is as.vetor(Zeta[,-K]) and
# in the second case phi is as.vector(Zeta[-m,])).

m    <- nrow(Zeta)
K    <- ncol(Zeta)

# Check whether last column of Zeta is all zeroes, or the
# last row is all zeroes.
if(all(Zeta[,K]==0) {
   case <- 1
} else if(all(Zeta[m,]==0)) {
   case <- 2
} else {
    stop("Either the last column or the last row of \"Zeta\" must be all zeroes.\n")
}

switch(EXPR=case,
    {npar <- m*(K-1)
    d1    <- array(0,c(m,K,npar))
    d2    <- array(0,c(m,K,npar,npar))
    E     <- exp(Zeta)
    delta <- diag(K)
    Dp1   <- apply(E,1,sum) # "Delta_i + 1"

    for(i in 1:m) {
        for(j in 1:K) {
            for(ell in 1:npar)
                jp <- ceiling(ell/K)
                ip <- ell - K*(jp-1)
                d1[i,j,ell] <- ... # dp_ij/dphi_ell
# Urrrrkkkkk.             
                ki <- (n-1)*(i-1) + k
                s  <- delta[j,k]*Dp1[i] - E[i,k]
                d1[i,j,ki] <- E[i,j]*s/Dp1[i]^2
                for(l in 1:(n-1)) {
                    li <- (n-1)*(i-1) + l 
                    ds <- delta[j,k]*E[i,l] - delta[k,l]*E[i,k]
                    a  <- (ds + delta[j,l]*s)*Dp1[i]
                    b  <- 2*E[i,l]*s
                    d2[i,j,ki,li] <- E[i,j]*(a-b)/Dp1[i]^3
            }
        }
    }
}
list(d1=d1,d2=d2)
}
