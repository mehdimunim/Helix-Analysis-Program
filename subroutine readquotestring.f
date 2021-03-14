      subroutine readquotestring(line,delim,i1,i2,ifail,len)
      character*(*) line
      character*1 delim
      i1=2
      ii=2
      ifail=0
      do while (line(ii:ii) .ne. delim .and. ii .lt. len)
          ii=ii+1
      end do
      if (line(ii:ii) .ne. delim) then
        ifail=1
        print *,'Missing closing quote'
      else
        i2=ii-1
      end if
      return
      end
