      subroutine unitmat(r)
      dimension r(3,3)
      do i=1,3
        do j=1,3
          r(i,j)=0.0
        end do
        r(i,i)=1.0
      end do
      return
      end
