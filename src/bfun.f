      subroutine bfun(fy,xispd,tpm,epsilon,n,nstate,wrk,beta)
      implicit double precision(a-h,o-z)
      dimension wrk(1), xispd(1)
      dimension fy(nstate,1), tpm(nstate,1), beta(nstate,1)
      one = 1.d0
      zero = 0.d0
      do 23000 j = 1,nstate 
      beta(j,n) = one
23000 continue
      do 23002 ktb = 2,n 
      kt = n - ktb + 1
      ktp = kt + 1
      tsum = zero
      do 23004 i = 1,nstate 
      wrk(i) = zero
      do 23006 j = 1,nstate 
      wrk(i) = wrk(i) + tpm(i,j)*beta(j,ktp)*fy(j,ktp)
23006 continue
      tsum = tsum + wrk(i)
23004 continue
      if(.not.(tsum .lt. epsilon))goto 23008
      do 23010 j = 1,nstate 
      beta(j,kt) = one/nstate
23010 continue
      goto 23009
23008 continue
      do 23012 j = 1,nstate 
      beta(j,kt) = wrk(j)/tsum
23012 continue
23009 continue
23002 continue
      return
      end
