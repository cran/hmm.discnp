      subroutine recurse(fy,xispd,tpm,nreps,epsilon,n,nstate,wrk,xlc, 
     *alpha,beta,gamma,xi,xisum)
      implicit double precision(a-h,o-z)
      dimension xispd(1), xlc(1)
      dimension tpm(nstate,1), wrk(nstate,1)
      dimension fy(nstate,1), alpha(nstate,1), beta(nstate,1), gamma(
     *nstate,1)
      dimension xi(nstate,nstate,1), xisum(nstate,1)
      zero = 0.d0
      do 23000 k = 1,nreps 
      kstart = (k-1)*n + 1
      call afun(fy(1,kstart),xispd,tpm,epsilon,n,nstate,wrk, xlc(kstart)
     *,alpha(1,kstart))
      call bfun(fy(1,kstart),xispd,tpm,epsilon,n,nstate,wrk, beta(1,
     *kstart))
      call gfun(alpha(1,kstart),beta(1,kstart),epsilon,n, nstate,wrk,
     *gamma(1,kstart))
      call xfun(alpha(1,kstart),beta(1,kstart),fy(1,kstart), tpm,
     *epsilon,n,nstate,wrk,xi(1,1,kstart-k+1))
23000 continue
      kstop = nreps*(n-1)
      do 23002 i = 1,nstate 
      do 23004 j = 1,nstate 
      xisum(i,j) = zero
      do 23006 k = 1,kstop 
      xisum(i,j) = xisum(i,j) + xi(i,j,k)
23006 continue
23004 continue
23002 continue
      return
      end
