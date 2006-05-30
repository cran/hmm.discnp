C Output from Public domain Ratfor, version 1.0
      subroutine gfun(alpha,beta,epsilon,n,nstate,wrk,gamma)
      implicit double precision(a-h,o-z)
      dimension alpha(nstate,1), beta(nstate,1), gamma(nstate,1)
      dimension wrk(1)
      zero = 0.d0
      ook = 1.d0/dble(nstate)
      do23000 kt = 1,n 
      tsum = zero
      do23002 i = 1,nstate 
      wrk(i) = alpha(i,kt)*beta(i,kt)
      tsum = tsum + wrk(i)
23002 continue
23003 continue
      if(tsum.lt.epsilon)then
      do23006 i = 1,nstate 
      gamma(i,kt) = ook
23006 continue
23007 continue
      else
      do23008 i = 1,nstate 
      gamma(i,kt) = wrk(i)/tsum
23008 continue
23009 continue
      endif
23000 continue
23001 continue
      return
      end
