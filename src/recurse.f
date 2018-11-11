C Output from Public domain Ratfor, version 1.03
      subroutine recurse(fy,xispd,tpm,nreps,epsilon,lns,nstate,nis,cis, 
     *wrk,xlc,ntot,nxi,alpha,beta,gamma,xi,xisum)
      implicit double precision(a-h,o-z)
      integer cis
      dimension xispd(nstate,nis), xlc(ntot), lns(nreps)
      dimension tpm(nstate,nstate), wrk(nstate,nstate)
      dimension fy(nstate,ntot), alpha(nstate,ntot), beta(nstate,ntot)
      dimension gamma(nstate,ntot), xi(nstate,nstate,nxi), xisum(nstate,
     *nstate)
      zero = 0.d0
      one = 1.d0
      kstop = 0
      do23000 k = 1,nreps 
      n = lns(k)
      nm1 = n - 1
      kstart = 1 + kstop
      if(cis .gt. 0)then
      kis = 1
      else
      kis = k
      endif
      call afun(fy(1,kstart),xispd(1,kis),tpm,epsilon,n,nstate,wrk, xlc(
     *kstart),alpha(1,kstart))
      call bfun(fy(1,kstart),xispd(1,kis),tpm,epsilon,n,nstate,wrk, beta
     *(1,kstart))
      call gfun(alpha(1,kstart),beta(1,kstart),epsilon,n, nstate,wrk,gam
     *ma(1,kstart))
      call xfun(alpha(1,kstart),beta(1,kstart),fy(1,kstart), tpm,epsilon
     *,n,nstate,nm1,wrk,xi(1,1,kstart-k+1))
      kstop = kstop + lns(k)
23000 continue
23001 continue
      kstop = kstop - nreps
      do23004 i = 1,nstate 
      do23006 j = 1,nstate 
      xisum(i,j) = zero
      do23008 k = 1,kstop 
      xisum(i,j) = xisum(i,j) + xi(i,j,k)
23008 continue
23009 continue
23006 continue
23007 continue
23004 continue
23005 continue
      return
      end
