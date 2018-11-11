subroutine gethgl(fy,y,ymiss,tpm,xispd,d1pi,d2pi,kstate,n,
                  npar,d1p,d2p,m,d1f,d2f,alpha,alphw,a,b,aw,bw,
                  xlc,ll,grad,hess)
implicit double precision(a-h,o-z)

double precision ll
integer   y(n)
integer   ymiss(n)
dimension fy(kstate,n)
dimension tpm(kstate,kstate), xispd(kstate)
dimension d1pi(kstate,npar), d2pi(kstate,npar,npar)
dimension d1p(kstate,kstate,npar), d2p(kstate,kstate,npar,npar)
dimension d1f(m,kstate,npar), d2f(m,kstate,npar,npar)
dimension alpha(kstate), alphw(kstate)
dimension a(kstate,npar), b(kstate,npar,npar)
dimension aw(kstate,npar), bw(kstate,npar,npar)
dimension xlc(n), grad(npar), hess(npar,npar)
#

# Set zero.
kt = 1
zero = 0.d0

# Initialize; i.e. do the t = 1 case:
sxlc = zero
do j = 1,kstate {
    alpha(j) = xispd(j)*fy(j,1)
    sxlc = sxlc + alpha(j)
    do k1 = 1,npar {
        if(ymiss(1) == 1) {
            d1fx1 = 0
        } else {
            d1fx1 = d1f(y(1),j,k1)
        }
        a(j,k1) = xispd(j)*d1fx1 + fy(j,1)*d1pi(j,k1)
        do k2 = 1,npar {
            if(ymiss(1) == 1) {
                d1fx2 = 0
            } else {
                d1fx2 = d1f(y(1),j,k2)
            }
            if(ymiss(1) == 1) {
                d2fx = 0
            } else {
                d2fx = d2f(y(1),j,k1,k2)
            }
            b(j,k1,k2) = (xispd(j)*d2fx + d1pi(j,k1)*d1fx2 +
                                          d1pi(j,k2)*d1fx1 +
                                          fy(j,1)*d2pi(j,k1,k2))
        }
    }
}
xlc(1) = sxlc

do j = 1,kstate {
    alpha(j) = alpha(j)/sxlc
}

if(n>1) {
do kt = 2,n {
# Do the b's:
    do j = 1,kstate {
        do k1 = 1,npar {
            if(ymiss(kt) == 1) {
                d1fx1 = 0
            } else {
                d1fx1 = d1f(y(kt),j,k1)
            }
            do k2 = 1,npar {
                if(ymiss(kt) == 1) {
                    d1fx2 = 0
                    d2fx  = 0
                } else {
                    d1fx2 = d1f(y(kt),j,k2)
                    d2fx  = d2f(y(kt),j,k1,k2)
                }
                vvv = zero
                xxx = zero
                yy1 = zero
                yy2 = zero
                zz1 = zero
                zz2 = zero
                www = zero
                do i = 1,kstate {
                    vvv = vvv+alpha(i)*d2p(i,j,k1,k2)
                    xxx = (xxx + a(i,k1)*d1p(i,j,k2) + a(i,k2)*d1p(i,j,k1) +
                                                       b(i,k1,k2)*tpm(i,j))
                    yy1 = yy1 + alpha(i)*d1p(i,j,k2)
                    yy2 = yy2 + a(i,k2)*tpm(i,j)
                    zz1 = zz1 + alpha(i)*d1p(i,j,k1)
                    zz2 = zz2 + a(i,k1)*tpm(i,j)
                    www = www + alpha(i)*tpm(i,j)
                }
                vvv = fy(j,kt)*vvv
                xxx = fy(j,kt)*xxx/sxlc
                yyy = d1fx1*(yy1 + yy2/sxlc)
                zzz = d1fx2*(zz1 + zz2/sxlc)
                www = d2fx*www
                bw(j,k1,k2) = vvv + xxx + yyy + zzz + www
            }
        }
    }

    do j = 1,kstate {
        do k1 = 1,npar {
            do k2 = 1,npar {
                b(j,k1,k2) = bw(j,k1,k2)
            }
        }
    }

# Do the a's:
    do j = 1,kstate {
        do k = 1,npar {
            if(ymiss(kt) == 1) {
                d1fx = 0
            } else {
                d1fx = d1f(y(kt),j,k)
            }
            xxx = zero
            yyy = zero
            zzz = zero
            do i = 1, kstate {
                xxx = xxx + alpha(i)*d1p(i,j,k)
                yyy = yyy + a(i,k)*tpm(i,j)
                zzz = zzz + alpha(i)*tpm(i,j)
            }
            aw(j,k) = fy(j,kt)*(xxx + yyy/sxlc) + d1fx*zzz
        }
    }
    do j = 1,kstate {
        do k = 1,npar {
            a(j,k) = aw(j,k)
        }
    }

# Do the alpha's:
    sxlc = zero
    do j = 1,kstate {
        alphw(j) = zero
        do i = 1,kstate {
            alphw(j) = alphw(j) + alpha(i)*tpm(i,j)
        }
        alphw(j) = fy(j,kt)*alphw(j)
        sxlc = sxlc + alphw(j)
    }
    xlc(kt) = sxlc
    do j = 1,kstate {
        alpha(j) = alphw(j)/sxlc
    }
}
}

# Finish off:

# Log likelihood.
ll = zero
do kt = 1,n {
   ll = ll + log(xlc(kt))
}

# Gradient.
do j = 1,npar {
   xxx = zero
   do i = 1,kstate {
       xxx = xxx + a(i,j)
   }
   grad(j) = xxx/sxlc
}

# Hessian.
do j1 = 1,npar {
    do j2 = 1,npar {
        xxx = zero
        yyy = zero
        zzz = zero
        do i = 1,kstate {
            xxx = xxx + b(i,j1,j2)
        }
        hess(j1,j2) = xxx/sxlc - grad(j1)*grad(j2)
    }
}

return
end
