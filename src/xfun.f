      subroutine xfun(alpha,beta,fy,tpm,epsilon,n,nstate,wrk,xi)
      implicit double precision(a-h,o-z)
      dimension alpha(nstate,1), beta(nstate,1), fy(nstate,1)
      dimension tpm(nstate,1), wrk(nstate,1), xi(nstate,nstate,1)
      one = 1.d0
      zero = 0.d0
      dns2 = dble(nstate*nstate)
      do 23000 ktp = 2,n 
      kt = ktp - 1
      tsum = zero
      do 23002 i = 1,nstate 
      do 23004 j = 1, nstate 
      wrk(i,j) = alpha(i,kt)*fy(j,ktp)*beta(j,ktp)*tpm(i,j)
      tsum = tsum + wrk(i,j)
23004 continue
23002 continue
      if(.not.(tsum.lt.epsilon))goto 23006
      do 23008 i = 1,nstate 
      do 23010 j = 1,nstate 
      xi(i,j,kt) = one/dns2
23010 continue
23008 continue
      goto 23007
23006 continue
      do 23012 i = 1,nstate 
      do 23014 j = 1,nstate 
      xi(i,j,kt) = wrk(i,j)/tsum
23014 continue
23012 continue
23007 continue
23000 continue
      return
      end
