subroutine getl(fy,tpm,xispd,kstate,n,alpha,alphw,xlc)
implicit double precision(a-h,o-z)

dimension fy(kstate,n)
dimension tpm(kstate,kstate), xispd(kstate)
dimension alpha(kstate), alphw(kstate)
dimension xlc(n)

# Set zero.
zero = 0.d0

# Initialize; i.e. do the t = 1 case:
sxlc = zero
do j = 1,kstate {
	alpha(j) = xispd(j)*fy(j,1)
	sxlc = sxlc + alpha(j)
}
xlc(1) = sxlc

do j = 1,kstate {
	alpha(j) = alpha(j)/sxlc
}

if(n>1) {
do kt = 2,n {
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
