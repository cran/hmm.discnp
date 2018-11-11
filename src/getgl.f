C Output from Public domain Ratfor, version 1.03
      subroutine getgl(fy,y,ymiss,tpm,xispd,d1pi,kstate,n, npar,d1p,m,d1
     *f,alpha,alphw,a,aw,xlc)
      implicit double precision(a-h,o-z)
      integer y(n)
      integer ymiss(n)
      dimension fy(kstate,n)
      dimension tpm(kstate,kstate), xispd(kstate)
      dimension d1pi(kstate,npar), d1p(kstate,kstate,npar)
      dimension d1f(m,kstate,npar), alpha(kstate), alphw(kstate)
      dimension a(kstate,npar), aw(kstate,npar), xlc(n)
      zero = 0.d0
      sxlc = zero
      do23000 j = 1,kstate 
      alpha(j) = xispd(j)*fy(j,1)
      sxlc = sxlc + alpha(j)
      do23002 k = 1,npar 
      if(ymiss(1) .eq. 1)then
      d1fx = zero
      else
      d1fx = d1f(y(1),j,k)
      endif
      a(j,k) = xispd(j)*d1fx + fy(j,1)*d1pi(j,k)
23002 continue
23003 continue
23000 continue
23001 continue
      xlc(1) = sxlc
      do23006 j = 1,kstate 
      alpha(j) = alpha(j)/sxlc
23006 continue
23007 continue
      if(n.gt.1)then
      do23010 kt = 2,n 
      kstart = 1+(kt-1)*kstate
      do23012 j = 1,kstate 
      do23014 k = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx = zero
      else
      d1fx = d1f(y(kt),j,k)
      endif
      xxx = zero
      yyy = zero
      zzz = zero
      do23018 i = 1, kstate 
      xxx = xxx + alpha(i)*d1p(i,j,k)
      yyy = yyy + a(i,k)*tpm(i,j)
      zzz = zzz + alpha(i)*tpm(i,j)
23018 continue
23019 continue
      aw(j,k) = fy(j,kt)*(xxx + yyy/sxlc) + d1fx*zzz
23014 continue
23015 continue
23012 continue
23013 continue
      do23020 j = 1,kstate 
      do23022 k = 1,npar 
      a(j,k) = aw(j,k)
23022 continue
23023 continue
23020 continue
23021 continue
      sxlc = zero
      do23024 j = 1,kstate 
      alphw(j) = zero
      do23026 i = 1,kstate 
      alphw(j) = alphw(j) + alpha(i)*tpm(i,j)
23026 continue
23027 continue
      alphw(j) = fy(j,kt)*alphw(j)
      sxlc = sxlc + alphw(j)
23024 continue
23025 continue
      xlc(kt) = sxlc
      do23028 j = 1,kstate 
      alpha(j) = alphw(j)/sxlc
23028 continue
23029 continue
23010 continue
23011 continue
      endif
      return
      end
