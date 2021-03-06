      def dvproj(x,y,z):
      #real*8 x(3),y(3),z(3)
c/*-
c     compute z = the projection on x of y
c     projxy requires that x be a unit vector
c     formula: projxy = (y dot x)*x
c-*/
      #real*8 tmp(3),ddot
      call dvset(tmp,x)
      call dvnorm(tmp)
      call dvmul(tmp,ddot(tmp,y),z)
      return
      # end
