      subroutine condenselist(index,len,incr,iout)
      dimension index(len)
      character*80 line
c     print *,'CONDENSELIST len=',len
      i=1
      i1=index(1)
      i2=index(1)
      ic=1
      do while (i .lt. len)
        i=i+1
        if (index(i) .eq. index(i-1)+1 .and. i .lt. len) then
          i2=i2+1
        else
c         Range found
          if (i .eq. len) i2=index(i)
          call printrange(line,i1,i2,ic,incr,iout)
          i1=index(i)
          i2=index(i)
        end if
      end do
      if (ic .gt. 1) write (iout,1000) line(1:ic-2)
      return
1000  format(1x,a)
      end
