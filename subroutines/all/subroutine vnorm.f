      subroutine vnorm(x)
      real x(3)
c     /*- normalize the vector x -*/
      real rmag,amag
      rmag=amag(x)
      if (rmag .lt. 1e-12) then
c        write(*,*) 'vnorm: can''t normalize zero vector'
         continue
      else
         x(1)=x(1)/rmag
         x(2)=x(2)/rmag
         x(3)=x(3)/rmag
      end if
      end
