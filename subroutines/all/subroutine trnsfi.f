      subroutine trnsfi(ia,ib,n)
c*****Fast integer array transfer
      dimension ia(n),ib(n)
      do i=1,n
        ia(i)=ib(i)
      end do
      return
      end
