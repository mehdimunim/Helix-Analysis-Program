      subroutine leftadjustline(line,ic1,ic2)
      character*(*) line
      icol=ic1
      call nextchar(line,icol,132)
      nshift=icol-ic1
      if (nshift .gt. 0) then
        do i=icol,ic2
          line(i-nshift:i-nshift)=line(i:i)
        end do
        do i=1,nshift
          line(ic2-i+1:ic2-i+1)=' '
        end do
      end if
      return
      end
