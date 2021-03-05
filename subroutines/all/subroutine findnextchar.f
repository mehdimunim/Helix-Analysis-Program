      subroutine findnextchar(char,line,ic,len)
      character*1 char
      character*(*) line
c     Finds the next character char in line
c     print *,'FINDCHAR char=',char
      do while (line(ic:ic) .ne. char .and. ic .lt. len)
        ic=ic+1
      end do
      return
      end
