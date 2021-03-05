      subroutine dvdiv(x,y,z)
      real*8 x(3), y, z(3)
c     /*- z =  x / y -*/
      integer i
      if (y .lt. 1e-12) then
         write(*,*) 'dvdiv: can''t divide vector by zero!'
      else
         do i=1,3
            z(i)=x(i)/y
         end do
      end if
      end
