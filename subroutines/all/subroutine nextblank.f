      subroutine nextblank(line,ifc,len)
      character*(*) line
c     Finds the next blank in line
      character*1 tab,ctrlM
      common /tab/ tab,ctrlM
      if (ifc .gt. len-1) return
      ifc1=ifc
      do i=ifc1,len-1
        ifc=i
        if (line(i:i) .eq. ' ' .or. line(i:i) .eq. tab) then
c         if (line(i:i) .eq. '!') ifc=len
          return
        end if
      end do
      ifc=len
      return
      end
