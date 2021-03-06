      subroutine calcperp(start,dir,from,orig,perpvec,itest)
      real*8 start,dir,from,orig,perpvec
      dimension start(3),dir(3),from(3),orig(3),perpvec(3)
c     For a line from start in the direction dir, calculate the normal to it
c     from the point from. The normal meets the line at orig and its direction
c     is perpvec
      real*8 ddot,dmag,dsx,xx,ddsx,perpvecfac
      dimension dsx(3),xx(3)
      call dvdif(from,start,dsx)
      ddsx=ddot(dir,dsx)
      do k=1,3
        perpvec(k)=dir(k)*ddsx-dsx(k)
      end do
c     call dcross(dsx,dir,xx)
c     call dcross(xx,dir,perpvec)
      call dvnorm(perpvec)
      perpvecfac=ddot(perpvec,dsx)
      do k=1,3
        orig(k)=from(k)-perpvecfac*perpvec(k)
      end do
      if (itest .gt. 0) then
        write (6,1000) 'Start     ',start
        write (6,1000) 'Dir       ',dir,' Magn=',dmag(dir)
        write (6,1000) 'From      ',from
        write (6,1000) 'Orig      ',orig
        write (6,1000) 'Perpvec   ',perpvec,' Magn=',dmag(dir)
        if (itest .gt. 1) then
          cc=ddot(dir,perpvec)
          print *,'Dir . Perpvec (->0)=',cc
          call dvdif(start,orig,xx)
          cc=ddot(dir,xx)/dmag(xx)
          print *,'Dir . (orig-start) (->1)=',cc
        end if
      end if
      return
1000  format(1x,a10,3f12.6,a,f12.6)
      end
