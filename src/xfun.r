subroutine xfun(alpha,beta,fy,tpm,epsilon,n,nstate,wrk,xi)
implicit double precision(a-h,o-z)
dimension alpha(nstate,1), beta(nstate,1), fy(nstate,1)
dimension tpm(nstate,1), wrk(nstate,1), xi(nstate,nstate,1)

one  = 1.d0
zero = 0.d0
dns2 = dble(nstate*nstate)

do ktp = 2,n {
	kt = ktp - 1
	tsum = zero
	do i = 1,nstate {
		do j = 1, nstate {
			wrk(i,j) = alpha(i,kt)*fy(j,ktp)*beta(j,ktp)*tpm(i,j)
			tsum = tsum + wrk(i,j)
		}
	}
	if(tsum<epsilon) {
		do i = 1,nstate {
			do j = 1,nstate {
				xi(i,j,kt) = one/dns2
			}
		}
	}
	else {
		do i = 1,nstate {
			do j = 1,nstate {
				xi(i,j,kt) = wrk(i,j)/tsum
			}
		}
	}
}

return
end