      subroutine blankout(line,n1,n2)
      character*(*) line
      do i=n1,n2
        line(i:i)=' '
      end do
      return
      end
