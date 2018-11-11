subroutine getgl(fy,y,ymiss,tpm,xispd,d1pi,kstate,n,
                 npar,d1p,m,d1f,alpha,alphw,a,aw,xlc)
implicit double precision(a-h,o-z)

integer   y(n)
integer   ymiss(n)
dimension fy(kstate,n)
dimension tpm(kstate,kstate), xispd(kstate)
dimension d1pi(kstate,npar), d1p(kstate,kstate,npar)
dimension d1f(m,kstate,npar), alpha(kstate), alphw(kstate)
dimension a(kstate,npar), aw(kstate,npar), xlc(n)

# Set zero.
zero = 0.d0

# Initialize; i.e. do the t = 1 case:
sxlc = zero
do j = 1,kstate {
	alpha(j) = xispd(j)*fy(j,1)
	sxlc = sxlc + alpha(j)
	do k = 1,npar {
                if(ymiss(1) == 1) {
                    d1fx = zero
                } else {
                    d1fx = d1f(y(1),j,k)
                }
		a(j,k) = xispd(j)*d1fx + fy(j,1)*d1pi(j,k)
	}
}
xlc(1) = sxlc

do j = 1,kstate {
	alpha(j) = alpha(j)/sxlc
}

if(n>1) {
do kt = 2,n {
	kstart = 1+(kt-1)*kstate
# Do the a's:
	do j = 1,kstate {
		do k = 1,npar {
                    if(ymiss(kt) == 1) {
                        d1fx = zero
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

return
end
