      subroutine zeroiti(ia,i0,n)
      dimension ia(n)
      do k=i0+1,n
        ia(k)=0
      end do
      return
      end
