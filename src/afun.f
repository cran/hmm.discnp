      subroutine afun(fy,xispd,tpm,epsilon,n,nstate,wrk,xlc,alpha)
      implicit double precision(a-h,o-z)
      dimension wrk(1), xispd(1), xlc(1)
      dimension fy(nstate,1), tpm(nstate,1), alpha(nstate,1)
      one = 1.d0
      zero = 0.d0
      dummy = -one
      tsum = zero
      do 23000 j = 1,nstate 
      wrk(j) = fy(j,1)*xispd(j)
      tsum = tsum + wrk(j)
23000 continue
      if(.not.(tsum .lt. epsilon))goto 23002
      xlc(1) = dummy
      do 23004 j = 1,nstate 
      alpha(j,1) = one/nstate
23004 continue
      goto 23003
23002 continue
      xlc(1) = tsum
      do 23006 j = 1,nstate 
      alpha(j,1) = wrk(j)/tsum
23006 continue
23003 continue
      do 23008 kt = 2,n 
      tsum = zero
      ktm = kt - 1
      do 23010 j = 1,nstate 
      wrk(j) = zero
      do 23012 i = 1,nstate 
      wrk(j) = wrk(j) + alpha(i,ktm)*tpm(i,j)
23012 continue
      wrk(j) = fy(j,kt)*wrk(j)
      tsum = tsum + wrk(j)
23010 continue
      if(.not.(tsum .lt. epsilon))goto 23014
      xlc(kt) = dummy
      do 23016 j = 1,nstate 
      alpha(j,kt) = one/nstate
23016 continue
      goto 23015
23014 continue
      xlc(kt) = tsum
      do 23018 j = 1,nstate 
      alpha(j,kt) = wrk(j)/tsum
23018 continue
23015 continue
23008 continue
      return
      end
