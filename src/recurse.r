subroutine recurse(fy,xispd,tpm,nreps,epsilon,n,nstate,wrk,xlc,
                   alpha,beta,gamma,xi,xisum)

implicit double precision(a-h,o-z)
dimension xispd(1), xlc(1)
dimension tpm(nstate,1), wrk(nstate,1)
dimension fy(nstate,1), alpha(nstate,1), beta(nstate,1), gamma(nstate,1)
dimension xi(nstate,nstate,1), xisum(nstate,1)

# Set zero:
zero = 0.d0

# Run through the replicates.
do k = 1,nreps {
	kstart = (k-1)*n + 1

# Update the alpha's.
	call afun(fy(1,kstart),xispd,tpm,epsilon,n,nstate,wrk,
                  xlc(kstart),alpha(1,kstart))

# Update the beta's.
	call bfun(fy(1,kstart),xispd,tpm,epsilon,n,nstate,wrk,
                  beta(1,kstart))

# Update the gamma's.
	call gfun(alpha(1,kstart),beta(1,kstart),epsilon,n,
                  nstate,wrk,gamma(1,kstart))

# Update the xi's.
	call xfun(alpha(1,kstart),beta(1,kstart),fy(1,kstart),
                  tpm,epsilon,n,nstate,wrk,xi(1,1,kstart-k+1))
}

# Sum up the xi's.
kstop = nreps*(n-1)
do i = 1,nstate {
	do j = 1,nstate {
		xisum(i,j) = zero
		do k = 1,kstop {
			xisum(i,j) = xisum(i,j) + xi(i,j,k)
		}
	}
}

return
end
