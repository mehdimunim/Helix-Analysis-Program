      subroutine trnsfrd(a,b,n)
c*****Fast array transfer
      real*8 a(n),b(n)
      do i=1,n
        a(i)=b(i)
      end do
      return
      end
