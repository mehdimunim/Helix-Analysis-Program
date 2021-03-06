      subroutine getseg4(segid4,line,isegcol1,nsegcol)
      character*4 segid4
      character*8 segid8
      character* 132 line
c     Get a max 4 character segment name
      if (nsegcol .le. 4) then
        segid4(1:nsegcol)=line(isegcol1:isegcol1+nsegcol-1)
      else
        segid8(1:nsegcol)=line(isegcol1:isegcol1+nsegcol-1)
        call leftadjustline(segid8,1,nsegcol)
        segid4(1:4)=segid8(1:4)
      end if
      return
      end
