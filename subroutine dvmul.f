      subroutine dvmul(x,y,z)
      real*8 x(3), y, z(3)
c     /*- z = x * y -*/
      z(1)=x(1)*y
      z(2)=x(2)*y
      z(3)=x(3)*y
      end
