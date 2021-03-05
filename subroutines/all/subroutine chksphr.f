      subroutine chksphr(c,ih,n,nnh,rmin,czero)
      dimension c(3,n),ih(n),czero(3)
c     Check if nothing is outside the sphere
      do ii=1,nnh
        i=ih(ii)
        rr=dist2(c(1,i),czero)
        if (rr-0.001 .gt. rmin**2) then
          print *,'PROGRAM ERROR: Atom ',i,' is outside the sphere'
          rr=sqrt(rr)
          print *,'  rmin=',rmin,' r=',rr
        end if
      end do
      return
      end
