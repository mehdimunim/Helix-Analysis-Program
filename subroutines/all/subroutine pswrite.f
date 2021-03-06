      subroutine pswrite(iout,i1,i2,lab,llab)
      character*(*) lab
c     Writes a PS command with the minimum number of blanks
      character*80 line
      icol=1
      call writeint(line,icol,i1,lenw)
      line(icol:icol)=' '
      icol=icol+1
      call writeint(line,icol,i2,lenw)
      line(icol:icol)=' '
      line(icol+1:icol+llab)=lab(1:llab)
      write (iout,1000) line(1:icol+llab)
      return
1000  format(a)
      end
