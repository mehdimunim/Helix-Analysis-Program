real*8 function ddistsq(x,x0)
      real*8 x(3),x0(3)
c     /*- returns squared distance -*/
      ddistsq=(x(1)-x0(1))**2 + (x(2)-x0(2))**2 + (x(3)-x0(3))**2
      return
      end