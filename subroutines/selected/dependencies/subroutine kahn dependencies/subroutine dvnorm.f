      subroutine dvnorm(x)
      real*8 x(3)
c     /*- normalize the vector x -*/
      real*8 rmag,dmag
      rmag=dmag(x)
      if (rmag .lt. 1e-12) then
         write(*,*) 'dvnorm: can''t normalize zero vector'
      else
         x(1)=x(1)/rmag
         x(2)=x(2)/rmag
         x(3)=x(3)/rmag
      end if
      end
