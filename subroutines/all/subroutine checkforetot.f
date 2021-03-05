      subroutine checkforetot(lmarker,line,iconf,etot,ietot,iverb)
      character* 132 line
c     print *,'DB ECHECK lmarker=',lmarker,' line=',line(1:40)
      icol=lmarker+1
      call nextchar(line,icol,80)
      if (icol .gt. 80) return
      icol1=icol
      call nextblank(line,icol,80)
      icol=0
      if (line(icol1:icol1+1) .eq. 'E=') then
        icol=icol1+2
      else if (line(icol1:icol1+6) .eq. 'energy=' .or.
     -    line(icol1:icol1+6) .eq. 'Energy=') then
        icol=icol1+7
      end if
      if (icol .gt. icol1) then
        call nextchar(line,icol,80)
        icol1=icol
        call nextblank(line,icol,80)
        read (line(icol1:icol-1),*,err=100) etot
        ietot=1
        if (iverb .gt. 0)
     -    print *,' Structure ',iconf,' Energy read=',etot
      end if
      return
100   print *,'ERROR: invalid energy value:',line(icol1:icol-1)
      return
      end
