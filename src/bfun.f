C Output from Public domain Ratfor, version 1.0
      subroutine bfun(fy,xispd,tpm,epsilon,n,nstate,wrk,beta)
      implicit double precision(a-h,o-z)
      dimension wrk(1), xispd(1)
      dimension fy(nstate,1), tpm(nstate,1), beta(nstate,1)
      one = 1.d0
      zero = 0.d0
      do23000 j = 1,nstate 
      beta(j,n) = one
23000 continue
23001 continue
      do23002 ktb = 2,n 
      kt = n - ktb + 1
      ktp = kt + 1
      tsum = zero
      do23004 i = 1,nstate 
      wrk(i) = zero
      do23006 j = 1,nstate 
      wrk(i) = wrk(i) + tpm(i,j)*beta(j,ktp)*fy(j,ktp)
23006 continue
23007 continue
      tsum = tsum + wrk(i)
23004 continue
23005 continue
      if(tsum .lt. epsilon)then
      do23010 j = 1,nstate 
      beta(j,kt) = one/nstate
23010 continue
23011 continue
      else
      do23012 j = 1,nstate 
      beta(j,kt) = wrk(j)/tsum
23012 continue
23013 continue
      endif
23002 continue
23003 continue
      return
      end
