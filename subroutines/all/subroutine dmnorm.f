      subroutine dmnorm(r)
c#    MMC routine 052 lstmod: 11/13/90
c*****Normalizes the matrix r
      real*8 r,rr
      dimension r(3,3)
      do i=1,3
        rr=0.0
        do k=1,3
          rr=rr+r(k,i)**2
        end do
        rr=sqrt(sngl(rr))
        do k=1,3
          r(k,i)=r(k,i)/rr
        end do
      end do
      return
      end
