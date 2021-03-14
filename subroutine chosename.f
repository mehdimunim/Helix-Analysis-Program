      subroutine chosename(namelist,nnames,lab,lablen,iname)
      character*(*) namelist(nnames),lab
100   write (6,1000) lab(1:lablen)
      do i=1,nnames
        if (i .lt. nnames) write (6,1001) namelist(i),i
        if (i .eq. nnames) write (6,1002) namelist(i),i
      end do
      call getint('Convention number',17,999999,1,nnames,iname,00)
      if (iname .lt. 1 .or. iname .gt. nnames) then
        print *,'Invalid choice'
        go to 100
      end if
      return
1000  format(' Possible ',a,' name conventions:')
1001  format(5x,a8,':',i3)
1002  format(5x,a8,':',i3,'  ',$)
      end
