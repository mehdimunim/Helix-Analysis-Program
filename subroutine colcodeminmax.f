      subroutine colcodeminmax(ips,ixright,iydown,nrep,ncode,maxcode,
     -  rmin,rmax)
      character*6 limline
c     Plot color code
c     print *,'COLCODEMINMAX ips,ncode,iydown=',ips,ncode,iydown
      iyd0=iydown-2
      iyd1=iydown+3
      iyd2=iydown+8
      absmax=amax1(abs(rmin),abs(rmax))
      if (nrep .le. 1) then
        write (ips,3004) 12
        do icol=0,ncode+1
          call rgbcolor(ips,9)
          ix0=ixright+(icol)*50
          write (ips,3001) ix0,-iyd2
          if (icol .le. ncode) then
            write (ips,3012) '<=',icol,ncode
            icoldrw=icol
          else
            write (ips,3012) '>',ncode,ncode
            icoldrw=maxcode+1
          end if
          write (ips,*) 'np'
          call rgbcolor(ips,icoldrw)
          write (ips,3001) ix0+35,-iyd1
          write (ips,3002) ix0+45,-iyd1
          write (ips,*) 'sk'
        end do
        iyd3=iyd2+15
        call rgbcolor(ips,9)
        do icol=0,ncode+1
          ix0=ixright+(icol)*50
          write (ips,3001) ix0,-iyd3
          rlim=rmin+float(icol)*(rmax-rmin)/float(ncode)
          if (icol .gt. ncode) rlim=rmax
          if (absmax .lt. 1.0) then
            write (limline,3013) rlim
          else if (absmax .lt. 10.0) then
            write (limline,3014) rlim
          else
            write (limline,3015) rlim
          end if
          if (icol .le. ncode) then
            write (ips,3016) '<=',limline
            icoldrw=icol
          else
            write (ips,3016) '>',limline
            icoldrw=maxcode+1
          end if
        end do
      end if
      call rgbcolor(ips,9)
      do icol=0,ncode+1
        ix0=ixright+(icol)*50
        call drawrect(-ips,9,1,ix0+35,ix0+45,-iyd2,-iyd0,nrep)
      end do
      return
3001  format(i4,1x,i5,' m')
3002  format(i4,1x,i5,' l')
3004  format(i5,' setlinewidth')
3012  format('(',a,i1,'/',i1,':) show')
3013  format(f6.3)
3014  format(f6.2)
3015  format(f6.1)
3016  format('(',a,a,') show')
      end
