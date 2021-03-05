      subroutine nextchar(line,ifc,len)
      character*(*) line
c     Finds the next nonblank in line
      character*1 tab,ctrlM
      common /tab/ tab,ctrlM
      if (ifc .gt. len-1) return
      ifc1=ifc
      do i=ifc1,len-1
        ifc=i
c       if (line(i:i) .eq. ' ') print *,'nextc i=',i,' blank'
c       if (line(i:i) .eq. tab) print *,'nextc i=',i,' tab'
        if (line(i:i) .ne. ' ' .and. line(i:i) .ne. tab) then
c         if (line(i:i) .eq. '!') ifc=len
          return
        end if
      end do
      ifc=len
      return
      end
