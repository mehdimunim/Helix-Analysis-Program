      subroutine rightadjustline(line,ic1,ic2)
      character* 132 line
c     print *,'RIGHTADJUSTLINE ic1,ic2=',ic1,ic2,' line=',line(ic1:ic2)
      nshift=0
      do while (line(ic2-nshift:ic2-nshift) .eq. ' '
     -         .and. nshift .le. ic2-ic1)
        nshift=nshift+1
      end do
      do ic=ic1,ic2-nshift
        icc=ic2-ic+ic1
        line(icc:icc)=line(icc-nshift:icc-nshift)
      end do
      do ic=ic1,ic1+nshift-1
        line(ic:ic)=' '
      end do
c     print *,'nshift=',nshift,' line=',line(ic1:ic2)
      return
      end
