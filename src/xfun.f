C Output from Public domain Ratfor, version 1.0
      subroutine xfun(alpha,beta,fy,tpm,epsilon,n,nstate,nm1,wrk,xi)
      implicit double precision(a-h,o-z)
      dimension alpha(nstate,n), beta(nstate,n), fy(nstate,n)
      dimension tpm(nstate,nstate), wrk(nstate,nstate), xi(nstate,nstate
     *,nm1)
      one = 1.d0
      zero = 0.d0
      dns2 = dble(nstate*nstate)
      do23000 ktp = 2,n 
      kt = ktp - 1
      tsum = zero
      do23002 i = 1,nstate 
      do23004 j = 1, nstate 
      wrk(i,j) = alpha(i,kt)*fy(j,ktp)*beta(j,ktp)*tpm(i,j)
      tsum = tsum + wrk(i,j)
23004 continue
23005 continue
23002 continue
23003 continue
      if(tsum.lt.epsilon)then
      do23008 i = 1,nstate 
      do23010 j = 1,nstate 
      xi(i,j,kt) = one/dns2
23010 continue
23011 continue
23008 continue
23009 continue
      else
      do23012 i = 1,nstate 
      do23014 j = 1,nstate 
      xi(i,j,kt) = wrk(i,j)/tsum
23014 continue
23015 continue
23012 continue
23013 continue
      endif
23000 continue
23001 continue
      return
      end
