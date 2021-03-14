      subroutine genimdist(crm,cell,ncell1,ncell,icmin,rmin2)
      dimension crm(3),cell(3,27)
c     print *,'GENIMDIST ncell,ncell1=',ncell,ncell1
      rmin2=1000000.0
      icmin=0
      do ic=ncell1,ncell
        r2=dist2(crm,cell(1,ic))
        if (r2 .lt. rmin2) then
          rmin2=r2
          icmin=ic
        end if
      end do
      if (icmin .eq. 0) then
        write (6,1717) crm,((i,cell(k,i),k=1,3),i=ncell1,ncell)
1717    format(' ERROR: ic=0 for crm=',3f10.5,/,(' cell',i3,'=',3f10.5))
        stop
      end if
      return
      end
