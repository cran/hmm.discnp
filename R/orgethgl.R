orgethgl <- function(fy,y,ymiss,tpm,xispd,d1pi,d2pi,npar,d1p,d2p,m,d1f,d2f) {
#
# Function orgethgl --- get Hessian, gradient, and log likelihood.
# for one observation sequence, original (for-loop-based) raw R method.
#

# Set zero.
zero = 0

# Create blank arrays.
kstate = nrow(tpm)
n      = length(y)
alpha  = numeric(kstate)
alphw  = numeric(kstate)
xlc    = numeric(n)
a      = matrix(NA,kstate,npar)
aw     = matrix(NA,kstate,npar)
b      = array(NA,dim=c(kstate,npar,npar))
bw     = array(NA,dim=c(kstate,npar,npar))

# Initialize; i.e. do the t = 1 case:
sxlc = zero
for(j in 1:kstate) {
    alpha[j] = xispd[j]*fy[j,1]
    sxlc = sxlc + alpha[j]
    for(k1 in 1:npar) {
        d1fx1 = if(is.na(y[1])) 0 else d1f[y[1],j,k1]
        a[j,k1] = xispd[j]*d1fx1 + fy[j,1]*d1pi[j,k1]
        for(k2 in 1:npar) {
                        d1fx2 = if(is.na(y[1])) 0 else d1f[y[1],j,k2]
                        d2fx  = if(is.na(y[1])) 0 else d2f[y[1],j,k1,k2]
                        b[j,k1,k2] = (xispd[j]*d2fx +
                                      d1pi[j,k1]*d1fx2  +
                                      d1pi[j,k2]*d1fx1 +
                                      fy[j,1]*d2pi[j,k1,k2])
        }
    }
}
xlc[1] = sxlc
for(j in 1:kstate) {
    alpha[j] = alpha[j]/sxlc
}

if(n>1) {
for(kt in 2:n) {
# Do the b's:
    for(j in 1:kstate) {
        for(k1 in 1:npar) {
            d1fx1 = if(is.na(y[kt])) 0 else d1f[y[kt],j,k1]
            for(k2 in 1:npar) {
                d1fx2 = if(is.na(y[kt])) 0 else d1f[y[kt],j,k2]
                d2fx  = if(is.na(y[kt])) 0 else d2f[y[kt],j,k1,k2]
                vvv   = zero
                xxx   = zero
                yy1   = zero
                yy2   = zero
                zz1   = zero
                zz2   = zero
                www   = zero
                for(i in 1:kstate) {
                    vvv = vvv+alpha[i]*d2p[i,j,k1,k2]
                    xxx = (xxx + a[i,k1]*d1p[i,j,k2] + a[i,k2]*d1p[i,j,k1] +
                                                       b[i,k1,k2]*tpm[i,j])
                    yy1 = yy1 + alpha[i]*d1p[i,j,k2]
                    yy2 = yy2 + a[i,k2]*tpm[i,j]
                    zz1 = zz1 + alpha[i]*d1p[i,j,k1]
                    zz2 = zz2 + a[i,k1]*tpm[i,j]
                    www = www + alpha[i]*tpm[i,j]
                }
                vvv = fy[j,kt]*vvv
                xxx = fy[j,kt]*xxx/sxlc
                yyy = d1fx1*(yy1 + yy2/sxlc)
                zzz = d1fx2*(zz1 + zz2/sxlc)
                www = d2fx*www
                bw[j,k1,k2] = vvv + xxx + yyy + zzz + www
            }
        }
    }
    for(j in 1:kstate) {
        for(k1 in 1:npar) {
            for(k2 in 1:npar) {
                b[j,k1,k2] = bw[j,k1,k2]
            }
        }
    }

# Do the a's:
    for(j in 1:kstate) {
        for(k in 1:npar) {
            d1fx = if(is.na(y[kt])) 0 else d1f[y[kt],j,k]
            xxx = zero
            yyy = zero
            zzz = zero
            for(i in 1:kstate) {
                xxx = xxx + alpha[i]*d1p[i,j,k]
                yyy = yyy + a[i,k]*tpm[i,j]
                zzz = zzz + alpha[i]*tpm[i,j]
            }
            aw[j,k] = fy[j,kt]*(xxx + yyy/sxlc) + d1fx*zzz
        }
    }
    for(j in 1:kstate) {
        for(k in 1:npar) {
            a[j,k] = aw[j,k]
        }
    }

# Do the alpha's:
    sxlc = zero
    for(j in 1:kstate) {
        alphw[j] = zero
        for(i in 1:kstate) {
            alphw[j] = alphw[j] + alpha[i]*tpm[i,j]
        }
        alphw[j] = fy[j,kt]*alphw[j]
        sxlc = sxlc + alphw[j]
    }
    xlc[kt] = sxlc
    for(j in 1:kstate) {
        alpha[j] = alphw[j]/sxlc
    }
}

# Finish off:
ll   <- sum(log(xlc))
grad <- apply(a,2,sum)/sxlc
hess <- matrix(NA,npar,npar)
for(k1 in 1:npar) {
    for(k2 in 1:npar) {
        xxx = zero
        yyy = zero
        zzz = zero
        for(i in 1:kstate) {
            xxx = xxx + b[i,k1,k2]
            yyy = yyy + a[i,k1]
            zzz = zzz + a[i,k2]
        }
        hess[k1,k2] = (xxx - yyy*zzz/sxlc)/sxlc
    }
}
}


list(ll=ll,grad=grad,hess=hess)
}
