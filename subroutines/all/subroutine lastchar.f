      subroutine lastchar(line,ifc,len)
      character*(*) line
c     Finds the last nonblank in line
      do i=1,len
        ifc=len-i+1
        if (line(ifc:ifc) .ne. ' ')  return
      end do
      if (ifc .eq. 1 .and. line(1:1) .eq. ' ') ifc=0
      return
      end
