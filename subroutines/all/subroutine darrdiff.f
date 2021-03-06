      subroutine darrdiff(a,b,c,n)
c*****c=a-b
      real*8 a(n),b(n),c(n)
      do i=1,n
        c(i)=a(i)-b(i)
      end do
      return
      end
