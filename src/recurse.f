C Output from Public domain Ratfor, version 1.0
      subroutine recurse(fy,xispd,tpm,nreps,epsilon,lns,nstate,wrk,xlc, 
     *ntot,nxi,alpha,beta,gamma,xi,xisum)
      implicit double precision(a-h,o-z)
      dimension xispd(nstate), xlc(ntot), lns(nreps)
      dimension tpm(nstate,nstate), wrk(nstate,nstate)
      dimension fy(nstate,ntot), alpha(nstate,ntot), beta(nstate,ntot)
      dimension gamma(nstate,ntot), xi(nstate,nstate,nxi), xisum(nstate,
     *nstate)
      zero = 0.d0
      kstop = 0
      do23000 k = 1,nreps 
      n = lns(k)
      nm1 = n - 1
      kstart = 1 + kstop
      call afun(fy(1,kstart),xispd,tpm,epsilon,n,nstate,wrk, xlc(kstart)
     *,alpha(1,kstart))
      call bfun(fy(1,kstart),xispd,tpm,epsilon,n,nstate,wrk, beta(1,ksta
     *rt))
      call gfun(alpha(1,kstart),beta(1,kstart),epsilon,n, nstate,wrk,gam
     *ma(1,kstart))
      call xfun(alpha(1,kstart),beta(1,kstart),fy(1,kstart), tpm,epsilon
     *,n,nstate,nm1,wrk,xi(1,1,kstart-k+1))
      kstop = kstop + lns(k)
23000 continue
23001 continue
      kstop = kstop - nreps
      do23002 i = 1,nstate 
      do23004 j = 1,nstate 
      xisum(i,j) = zero
      do23006 k = 1,kstop 
      xisum(i,j) = xisum(i,j) + xi(i,j,k)
23006 continue
23007 continue
23004 continue
23005 continue
23002 continue
23003 continue
      return
      end
