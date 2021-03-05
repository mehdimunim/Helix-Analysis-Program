      subroutine polar(x,y,z,r,theta,phi)
c     Return polar coordinates for point (x,y,z)
      real*8 x,y,z,r,theta,phi
      data PI /3.1415926/
      real*8 cosa
c     NOTE: Unlike most conventions, phi is the angle between R and the Z axis

      r=(x**2 + y**2 + z**2)**0.5
      if ((x .ne. 0.0) .and. (y .ne. 0.0)) then
         theta=atan(abs(y/x))
         if (y .gt. 0.0) then
            if (x .lt. 0.0) theta=PI-theta ! quadrant 2
         else
            if (x .lt. 0.0) theta=PI+theta ! quadrant 3
            if (x .gt. 0.0) theta=2.0*PI-theta ! quadrant 4
         end if
      else
         if (x .eq. 0.0) then
            if (y .ge. 0.0) then
               theta=PI/2.0
            else
               theta=3*PI/2.0
            end if
         end if
         if (y .eq. 0.0) then
            if (x .ge. 0) then
               theta=0.0
            else
               theta=PI
            end if
         end if
      end if
      if (r .eq. 0.0) then
         phi = 0.0
      else
         cosa=dble(z)/dble(r)
         phi=dacoscheck(cosa,ccc,1,6,'POLAR')
      end if
      end
