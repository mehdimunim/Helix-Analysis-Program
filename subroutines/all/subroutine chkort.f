      subroutine chkort(rot)
      dimension rot(3,3),rt(3,3)
c     Check rot for orthogonality
      do i=1,3
        do j=1,3
          sum=0.0
          do k=1,3
            sum=sum+rot(i,k)*rot(j,k)
          end do
          rt(i,j)=sum
        end do
      end do
      write (6,1000) rt
1000  format(' Orthogonality check:'/,(3f10.6))
      return
      end
