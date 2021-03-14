      subroutine nextstring(line,ifc,ic1,ic2,len)
      character*(*) line
      call nextchar(line,ifc,len)
      if (line(ifc:ifc) .eq. '"' .or. line(ifc:ifc) .eq. "'") then
        ifc0=ifc
c       Look for closing quote
        ifc=ifc+1
        ic1=ifc
        if (line(ifc0:ifc0) .eq. '"') then
          call findnextchar('"',line,ifc,len)
        else
          call findnextchar("'",line,ifc,len)
        end if
        if (ifc .eq. len) then
          print *,'No closing quote was found'
          stop
        end if
        ic2=ifc-1
        ifc=ifc+1
      else
        ic1=ifc
        call nextblank(line,ifc,len)
        ic2=ifc-1
      end if
      return
      end
