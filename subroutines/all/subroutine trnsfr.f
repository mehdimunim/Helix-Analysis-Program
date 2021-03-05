      subroutine trnsfr(a,b,n)
c*****Real array transfer
      dimension a(n),b(n)
      do i=1,n
        a(i)=b(i)
      end do
      return
      end
