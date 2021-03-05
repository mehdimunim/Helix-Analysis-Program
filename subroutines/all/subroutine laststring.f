      subroutine laststring(line,icf,icl,lline,len)
      character*(*) line
      ic=1
      icl=1
      do while (icl .lt. lline)
        call nextstring(line,ic,icf,icl,len)
        ic=icl+1
      end do
      return
      end
