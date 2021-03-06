      subroutine writevertices(vp,n,ioppbc,ioutpdb,iconntyp,icrot,crot)
      character*3 pbcres
      dimension vp(3,n),crot(3,3)
      common /pbcresname/ pbcres(10)
      real*8 c8(3)
c     print *,'WRITEVERTICES n=',n,' icrot=',icrot
c     print *,'crot=',crot
      if (iconntyp .eq. 1) return
      if (icrot .eq. 1) call rotate_c(vp,n,crot,vp,'WVERTICES',9)
      do ia=1,n
        do k=1,3
          c8(k)=vp(k,ia)
        end do
        call writepdbd(ioutpdb,c8,ia,1,'VERT',pbcres(ioppbc),'V',1.0,0.)
      end do
      return
      end
