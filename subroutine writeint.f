      subroutine writeint(line,icol,int,len)
      character*(*) line
c     Writes the integer int into line starting at icol
      ii=int
      if (ii .lt. 0) ii=-ii
      if (ii .eq. 0) then
        len=1
      else
        rii=ii+0.00001
        len=alog10(rii)+1
      end if
      if (int .lt. 0) len=len+1
      if (len .eq. 1) then
        write (line(icol:icol),101) int
      else if (len .eq. 2) then
        write (line(icol:icol+1),102) int
      else if (len .eq. 3) then
        write (line(icol:icol+2),103) int
      else if (len .eq. 4) then
        write (line(icol:icol+3),104) int
      else if (len .eq. 5) then
        write (line(icol:icol+4),105) int
      else if (len .eq. 6) then
        write (line(icol:icol+5),106) int
      else if (len .eq. 7) then
        write (line(icol:icol+6),107) int
      else if (len .eq. 8) then
        write (line(icol:icol+7),108) int
      else if (len .eq. 9) then
        write (line(icol:icol+8),109) int
      end if
      icol=icol+len
      return
101   format(i1)
102   format(i2)
103   format(i3)
104   format(i4)
105   format(i5)
106   format(i6)
107   format(i7)
108   format(i8)
109   format(i9)
      end
