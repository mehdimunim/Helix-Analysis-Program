      subroutine rotabout(v,v0,t,axis)
c     rotate vector v about point v0 by t radians along the axis
c     specified.
      real*8 v(3),v0(3),t,temp(3)
      character*1 axis

      st=sin(t)
      ct=cos(t)
      if ((axis .eq. 'x') .or. (axis .eq. 'X')) then
         temp(1)=v(1)
         temp(2)=v0(2) + v(2)*ct - v0(2)*ct - v(3)*st + v0(3)*st
         temp(3)=v0(3) + v(3)*ct - v0(3)*ct + v(2)*st - v0(2)*st
      else
         if ((axis .eq. 'y') .or. (axis .eq. 'Y')) then
            temp(1)=v0(1) + v(1)*ct - v0(1)*ct + v(3)*st - v0(3)*st
            temp(2)=v(2)
            temp(3)=v0(3) + v(3)*ct - v0(3)*ct - v(1)*st + v0(1)*st
         else
            if ((axis .eq. 'z') .or. (axis .eq. 'Z')) then
               temp(1)=v0(1) + v(1)*ct - v0(1)*ct - v(2)*st + v0(2)*st
               temp(2)=v0(2) + v(2)*ct - v0(2)*ct + v(1)*st - v0(1)*st
               temp(3)=v(3)
            else
               return           !no rotation if wrong axis specified
            end if
         end if
      end if
      do i=1,3
         v(i)=temp(i)
      end do
      end
