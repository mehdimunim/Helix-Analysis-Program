      subroutine writeline(iout,line,icol1,icol2,nocr)
      character* 132 line
      lc=icol2
      if (icol2 .eq. 0) call lastchar(line,lc,132)
      if (nocr .eq. 0) then
        write (iout,2000) line(icol1:lc)
      else
        write (iout,1000) line(icol1:lc)
      end if
      return
1000  format(a,$)
2000  format(a)
      end
