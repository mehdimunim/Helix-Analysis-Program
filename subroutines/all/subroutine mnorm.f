      subroutine mnorm(r)
C*****Normalizes the matrix r
      dimension r(3,3)
      do i=1,3
        rr=0.0
        do k=1,3
          rr=rr+r(k,i)**2
        end do
        rr=sqrt(rr)
        do k=1,3
          r(k,i)=r(k,i)/rr
        end do
      end do
      return
      end
