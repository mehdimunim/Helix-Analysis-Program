      subroutine genimdist123dim(crm,cell,ncell1,ncell,ixyzexcld,
     -  ixyzincld,icmin,rmin2)
      dimension crm(3),cell(3,27)
c     print *,'GENIMDIST123 ncell,ncell1,ixyzexcld,ixyzincld=',
c    -  ncell,ncell1,ixyzexcld,ixyzincld
      if (ixyzincld+ixyzexcld .eq. 0) then
        call genimdist(crm,cell,ncell1,ncell,icmin,rmin2)
      else if (ixyzexcld .gt. 0) then
c       2D PBC
        rmin2=1000000.0
        do ic=ncell1,ncell
c         Compare only with cells in the 2D plane
          if (cell(ixyzexcld,ic) .eq. 0.0) then
            r2=dist2(crm,cell(1,ic))
            if (r2 .lt. rmin2) then
              rmin2=r2
              icmin=ic
            end if
          end if
        end do
      else
c       1D PBC
        rmin2=1000000.0
        do ic=ncell1,ncell
c         Compare only with cells along the 1D axis
          nz=0
          do k=1,3
            if (k .ne. ixyzincld) then
              if (cell(k,ic) .eq. 0.0) nz=nz+1
            end if
          end do
          if (nz .eq. 2) then
            r2=(crm(ixyzincld)-cell(ixyzincld,ic))**2
            if (r2 .lt. rmin2) then
              rmin2=r2
              icmin=ic
            end if
          end if
        end do
      end if
c     print *,'ICMIN=',icmin
c      write (77,1717) ncell1,ncell,crm,(cell(k,1),k=1,3),rmin2
c1717  format(' ncell1,ncell=',2i3,' crm=',3f10.5,' cell=',3f10.5,
c    -  ' rmin2=',f10.5)
      return
      end
