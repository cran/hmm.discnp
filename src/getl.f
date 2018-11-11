C Output from Public domain Ratfor, version 1.03
      subroutine getl(fy,tpm,xispd,kstate,n,alpha,alphw,xlc)
      implicit double precision(a-h,o-z)
      dimension fy(kstate,n)
      dimension tpm(kstate,kstate), xispd(kstate)
      dimension alpha(kstate), alphw(kstate)
      dimension xlc(n)
      zero = 0.d0
      sxlc = zero
      do23000 j = 1,kstate 
      alpha(j) = xispd(j)*fy(j,1)
      sxlc = sxlc + alpha(j)
23000 continue
23001 continue
      xlc(1) = sxlc
      do23002 j = 1,kstate 
      alpha(j) = alpha(j)/sxlc
23002 continue
23003 continue
      if(n.gt.1)then
      do23006 kt = 2,n 
      sxlc = zero
      do23008 j = 1,kstate 
      alphw(j) = zero
      do23010 i = 1,kstate 
      alphw(j) = alphw(j) + alpha(i)*tpm(i,j)
23010 continue
23011 continue
      alphw(j) = fy(j,kt)*alphw(j)
      sxlc = sxlc + alphw(j)
23008 continue
23009 continue
      xlc(kt) = sxlc
      do23012 j = 1,kstate 
      alpha(j) = alphw(j)/sxlc
23012 continue
23013 continue
23006 continue
23007 continue
      endif
      return
      end
