C Output from Public domain Ratfor, version 1.03
      subroutine gethgl(fy,y,ymiss,tpm,xispd,d1pi,d2pi,kstate,n, npar,d1
     *p,d2p,m,d1f,d2f,alpha,alphw,a,b,aw,bw, xlc,ll,grad,hess)
      implicit double precision(a-h,o-z)
      double precision ll
      integer y(n)
      integer ymiss(n)
      dimension fy(kstate,n)
      dimension tpm(kstate,kstate), xispd(kstate)
      dimension d1pi(kstate,npar), d2pi(kstate,npar,npar)
      dimension d1p(kstate,kstate,npar), d2p(kstate,kstate,npar,npar)
      dimension d1f(m,kstate,npar), d2f(m,kstate,npar,npar)
      dimension alpha(kstate), alphw(kstate)
      dimension a(kstate,npar), b(kstate,npar,npar)
      dimension aw(kstate,npar), bw(kstate,npar,npar)
      dimension xlc(n), grad(npar), hess(npar,npar)
      kt = 1
      zero = 0.d0
      sxlc = zero
      do23000 j = 1,kstate 
      alpha(j) = xispd(j)*fy(j,1)
      sxlc = sxlc + alpha(j)
      do23002 k1 = 1,npar 
      if(ymiss(1) .eq. 1)then
      d1fx1 = 0
      else
      d1fx1 = d1f(y(1),j,k1)
      endif
      a(j,k1) = xispd(j)*d1fx1 + fy(j,1)*d1pi(j,k1)
      do23006 k2 = 1,npar 
      if(ymiss(1) .eq. 1)then
      d1fx2 = 0
      else
      d1fx2 = d1f(y(1),j,k2)
      endif
      if(ymiss(1) .eq. 1)then
      d2fx = 0
      else
      d2fx = d2f(y(1),j,k1,k2)
      endif
      b(j,k1,k2) = (xispd(j)*d2fx + d1pi(j,k1)*d1fx2 + d1pi(j,k2)*d1fx1 
     *+ fy(j,1)*d2pi(j,k1,k2))
23006 continue
23007 continue
23002 continue
23003 continue
23000 continue
23001 continue
      xlc(1) = sxlc
      do23012 j = 1,kstate 
      alpha(j) = alpha(j)/sxlc
23012 continue
23013 continue
      if(n.gt.1)then
      do23016 kt = 2,n 
      do23018 j = 1,kstate 
      do23020 k1 = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx1 = 0
      else
      d1fx1 = d1f(y(kt),j,k1)
      endif
      do23024 k2 = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx2 = 0
      d2fx = 0
      else
      d1fx2 = d1f(y(kt),j,k2)
      d2fx = d2f(y(kt),j,k1,k2)
      endif
      vvv = zero
      xxx = zero
      yy1 = zero
      yy2 = zero
      zz1 = zero
      zz2 = zero
      www = zero
      do23028 i = 1,kstate 
      vvv = vvv+alpha(i)*d2p(i,j,k1,k2)
      xxx = (xxx + a(i,k1)*d1p(i,j,k2) + a(i,k2)*d1p(i,j,k1) + b(i,k1,k2
     *)*tpm(i,j))
      yy1 = yy1 + alpha(i)*d1p(i,j,k2)
      yy2 = yy2 + a(i,k2)*tpm(i,j)
      zz1 = zz1 + alpha(i)*d1p(i,j,k1)
      zz2 = zz2 + a(i,k1)*tpm(i,j)
      www = www + alpha(i)*tpm(i,j)
23028 continue
23029 continue
      vvv = fy(j,kt)*vvv
      xxx = fy(j,kt)*xxx/sxlc
      yyy = d1fx1*(yy1 + yy2/sxlc)
      zzz = d1fx2*(zz1 + zz2/sxlc)
      www = d2fx*www
      bw(j,k1,k2) = vvv + xxx + yyy + zzz + www
23024 continue
23025 continue
23020 continue
23021 continue
23018 continue
23019 continue
      do23030 j = 1,kstate 
      do23032 k1 = 1,npar 
      do23034 k2 = 1,npar 
      b(j,k1,k2) = bw(j,k1,k2)
23034 continue
23035 continue
23032 continue
23033 continue
23030 continue
23031 continue
      do23036 j = 1,kstate 
      do23038 k = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx = 0
      else
      d1fx = d1f(y(kt),j,k)
      endif
      xxx = zero
      yyy = zero
      zzz = zero
      do23042 i = 1, kstate 
      xxx = xxx + alpha(i)*d1p(i,j,k)
      yyy = yyy + a(i,k)*tpm(i,j)
      zzz = zzz + alpha(i)*tpm(i,j)
23042 continue
23043 continue
      aw(j,k) = fy(j,kt)*(xxx + yyy/sxlc) + d1fx*zzz
23038 continue
23039 continue
23036 continue
23037 continue
      do23044 j = 1,kstate 
      do23046 k = 1,npar 
      a(j,k) = aw(j,k)
23046 continue
23047 continue
23044 continue
23045 continue
      sxlc = zero
      do23048 j = 1,kstate 
      alphw(j) = zero
      do23050 i = 1,kstate 
      alphw(j) = alphw(j) + alpha(i)*tpm(i,j)
23050 continue
23051 continue
      alphw(j) = fy(j,kt)*alphw(j)
      sxlc = sxlc + alphw(j)
23048 continue
23049 continue
      xlc(kt) = sxlc
      do23052 j = 1,kstate 
      alpha(j) = alphw(j)/sxlc
23052 continue
23053 continue
23016 continue
23017 continue
      endif
      ll = zero
      do23054 kt = 1,n 
      ll = ll + log(xlc(kt))
23054 continue
23055 continue
      do23056 j = 1,npar 
      xxx = zero
      do23058 i = 1,kstate 
      xxx = xxx + a(i,j)
23058 continue
23059 continue
      grad(j) = xxx/sxlc
23056 continue
23057 continue
      do23060 j1 = 1,npar 
      do23062 j2 = 1,npar 
      xxx = zero
      yyy = zero
      zzz = zero
      do23064 i = 1,kstate 
      xxx = xxx + b(i,j1,j2)
23064 continue
23065 continue
      hess(j1,j2) = xxx/sxlc - grad(j1)*grad(j2)
23062 continue
23063 continue
23060 continue
23061 continue
      return
      end
