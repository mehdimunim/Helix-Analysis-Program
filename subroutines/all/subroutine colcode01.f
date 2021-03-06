      subroutine colcode01(ips,ixright,iydown,ncode,nrep)
c     Plot color code
      iyd0=iydown-2
      iyd1=iydown+3
      iyd2=iydown+8
      if (nrep .le. 1) then
        write (ips,3004) 12
        do icol=1,ncode
          call rgbcolor(ips,9)
          ix0=ixright+(icol-1)*50
          write (ips,3001) ix0,-iyd2
          write (ips,3012) ncode-icol,ncode
          write (ips,*) 'np'
          call rgbcolor(ips,icol)
          write (ips,3001) ix0+35,-iyd1
          write (ips,3002) ix0+45,-iyd1
          write (ips,*) 'sk'
        end do
      end if
      do icol=1,ncode
        call drawrect(-ips,9,1,ix0+35,ix0+45,-iyd2,-iyd0,nrep)
      end do
      return
3001  format(i4,i5,' m')
3002  format(i4,i5,' l')
3004  format(i5,' setlinewidth')
3012  format('(>',i1,'/',i1,':) show')
      end
