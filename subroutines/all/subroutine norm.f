      subroutine norm(d,fac)
      dimension d(3)
c     Normalize d
      sum=0.0
      do k=1,3
        sum=sum+d(k)**2
      end do
      if (sum .gt. 0.0) then
        do k=1,3
          d(k)=fac*d(k)/sqrt(sum)
        end do
      end if
      return
      end
