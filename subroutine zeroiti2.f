      subroutine zeroiti2(ia,i0,n)
      integer*2 ia(n)
      do k=i0+1,n
        ia(k)=0
      end do
      return
      end
