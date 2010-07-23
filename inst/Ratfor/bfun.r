subroutine bfun(fy,xispd,tpm,epsilon,n,nstate,wrk,beta)
implicit double precision(a-h,o-z)
dimension wrk(1), xispd(1)
dimension fy(nstate,1), tpm(nstate,1), beta(nstate,1)

# Set some constants.
one  = 1.d0
zero = 0.d0

# Set the last beta's.

do j = 1,nstate {
	beta(j,n) = one # Doesn't matter that the beta_n's are
}                       # not rescaled.

# Run through the remaining n-1 of the betas (recursing!), backwards.
do ktb = 2,n {
	kt  = n - ktb + 1
	ktp = kt + 1
	tsum = zero
	do i = 1,nstate {
		wrk(i) = zero
		do j = 1,nstate {
			wrk(i) = wrk(i) + tpm(i,j)*beta(j,ktp)*fy(j,ktp)
		}
		tsum = tsum + wrk(i)
	}
	if(tsum < epsilon) {
		do j = 1,nstate {
			beta(j,kt) = one/nstate
		}
	}
	else {
		do j = 1,nstate {
			beta(j,kt) = wrk(j)/tsum
		}
	}
}

return
end