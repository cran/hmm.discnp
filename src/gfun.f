      subroutine gfun(alpha,beta,epsilon,n,nstate,wrk,gamma)
      implicit double precision(a-h,o-z)
      dimension alpha(nstate,1), beta(nstate,1), gamma(nstate,1)
      dimension wrk(1)
      zero = 0.d0
      ook = 1.d0/dble(nstate)
      do 23000 kt = 1,n 
      tsum = zero
      do 23002 i = 1,nstate 
      wrk(i) = alpha(i,kt)*beta(i,kt)
      tsum = tsum + wrk(i)
23002 continue
      if(.not.(tsum.lt.epsilon))goto 23004
      do 23006 i = 1,nstate 
      gamma(i,kt) = ook
23006 continue
      goto 23005
23004 continue
      do 23008 i = 1,nstate 
      gamma(i,kt) = wrk(i)/tsum
23008 continue
23005 continue
23000 continue
      return
      end
