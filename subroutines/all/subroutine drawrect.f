      subroutine drawrect(ips0,kc,lw,ix0,ix1,iy0,iy1,nrep)
c     Draw a rectangle
c      write (6,1000) kc,ix0,ix1,iy0,iy1
c1000  format(' DRAWRECT kc=',i1,' ix1,2=',2i5,' iy1,2=',2i5)
      ips=iabs(ips0)
      if (nrep .le. 1) then
        write (ips,3004) lw
        write (ips,3000)
        call rgbcolor(ips,kc)
        call pswrite(ips,ix0,iy0,'m',1)
        call pswrite(ips,ix0,iy1,'l',1)
        call pswrite(ips,ix1,iy1,'l',1)
        call pswrite(ips,ix1,iy0,'l',1)
        call pswrite(ips,ix0,iy0,'l',1)
        write(ips,3001)
      end if
      return
3000  format('np')
3001  format('sk')
3004  format(i5,' lw')
      end
