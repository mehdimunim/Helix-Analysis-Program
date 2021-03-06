      subroutine matprod(r1,r2,r3)
      dimension r1(3,3),r2(3,3),r3(3,3),r(3,3)
      do i=1,3
        do j=1,3
          rr=0.0
          do k=1,3
            rr=rr+r1(i,k)*r2(k,j)
          end do
          r(i,j)=rr
        end do
      end do
      call trnsfr(r3,r,9)
      return
      end
