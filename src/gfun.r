subroutine gfun(alpha,beta,epsilon,n,nstate,wrk,gamma)
implicit double precision(a-h,o-z)
dimension alpha(nstate,1), beta(nstate,1), gamma(nstate,1)
dimension wrk(1)

zero = 0.d0
ook  = 1.d0/dble(nstate)

do kt = 1,n {
	tsum = zero
	do i = 1,nstate {
		wrk(i) = alpha(i,kt)*beta(i,kt)
		tsum = tsum + wrk(i)
	}
	if(tsum<epsilon) {
		do i = 1,nstate {
			gamma(i,kt) = ook
		}
	}
	else {
		do i = 1,nstate {
			gamma(i,kt) = wrk(i)/tsum
		}
	}
}

return
end
