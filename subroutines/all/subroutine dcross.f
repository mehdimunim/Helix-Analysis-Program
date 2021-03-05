      subroutine dcross(x,y,z)
      real*8 x(3), y(3), z(3)
c     /*- compute z = x cross y -*/
      z(1)= x(2)*y(3)-y(2)*x(3)
      z(2)=-x(1)*y(3)+y(1)*x(3)
      z(3)= x(1)*y(2)-y(1)*x(2)
      end
