      subroutine leftadjustn(in,out,n)
c*****Left-adjust a string of n characters
      character*(*) in
      character*(*) out
      nz=0
      i=1
      do while(i .lt. n .and. in(i:i) .eq. ' ')
        i=i+1
      end do
      if (i .eq. n .and. in(n:n) .eq. ' ') then
        call blankout(out,1,n)
        return
      end if
      nz=i-1
      do j=nz+1,n
        out(j-nz:j-nz)=in(j:j)
      end do
      do j=n-nz+1,n
        out(j:j)=' '
      end do
      return
      end
