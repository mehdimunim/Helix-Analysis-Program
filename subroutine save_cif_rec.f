      subroutine save_cif_rec(label,llabel,clabel,lclabel,limcol,
     -  line,linein,ic1,ic2,len,noadjust,ierr)
      dimension limcol(2)
      character*(*) label,clabel
      character*132 line
      character*200 linein
c     print *,'SAVE_CIF_REC ',label(1:llabel),' ic1,i2=',ic1,ic2
c     print *,' REC:',linein(limcol(1):limcol(2)),' LIM:',limcol
c     if (ic1 .le. 111 .and. ic2 .ge. 111 .and.
c    -  label(1:llabel) .ne. 'Alt at id ') then
c       print *,label(1:llabel),' ic1,ic2=',ic1,ic2
c       stop
c     end if
      if (ic2-ic1+1 .lt. len) then
        write (6,2000) label(1:llabel),len,ic1,ic2,clabel(1:lclabel)
        len=ic2-ic1+1
        ierr=ierr+1
      end if
      line(ic1:ic1+len-1)=linein(limcol(1):limcol(2))
      if (noadjust .eq. 0) call rightadjustline(line,ic1,ic2)
      return
2000  format(' ERROR: record length for ',a,' (',i2,') exceeds range:',
     -  i3,' - ',i3,/,' Record will be truncated',/,
     -  ' Modify the corresponding ',a,' values in block data')
      end
