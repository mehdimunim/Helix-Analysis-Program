      subroutine listformat(ityp)
      character*11 formatname
      common /formats/ iqdconv(20),formatname(19)
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      character*1 ans
      if (ityp .eq. 0) then
        call quiz(ans,ians,'d',' ',0,'file format',11,0,5,6,0)
        it=iqdconv(ians)
      else
        it=ityp
      end if
      write (6,1000) formatname(it)
      write (6,1002) incol(it)
      if (iinamcol(2,it) .ge. iinamcol(1,it))
     -   write (6,1001) 'atom name:           ',(iinamcol(k,it),k=1,2)
      if (iirescol(2,it) .ge. iirescol(1,it))
     -   write (6,1001) 'residue name:        ',(iirescol(k,it),k=1,2)
      if (iiseqncol(2,it) .ge. iiseqncol(1,it))
     -   write (6,1001) 'atom number:         ',(iiseqncol(k,it),k=1,2)
      if (iiresncol(2,it) .ge. iiresncol(1,it))
     -   write (6,1001) 'residue number:      ',(iiresncol(k,it),k=1,2)
      if (iisegcol(2,it) .ge. iisegcol(1,it))
     -   write (6,1001) 'segment (chain) ID:  ',(iisegcol(k,it),k=1,2)
      if (iiresidcol(2,it) .ge. iiresidcol(1,it))
     -   write (6,1001) 'residue ID (Charmm): ',(iiresidcol(k,it),k=1,2)
      if (iiccol(2,it) .ge. iiccol(1,it))
     -   write (6,1001) 'X,Y,Z coordinates:   ',(iiccol(k,it),k=1,2)
      if (iiqcol(2,it) .ge. iiqcol(1,it))
     -   write (6,1001) 'atomic charges:      ',(iiqcol(k,it),k=1,2)
      if (iipotcol(2,it) .ge. iipotcol(1,it))
     -   write (6,1001) 'potential identifier:',(iipotcol(k,it),k=1,2)
      if (iichemcol(2,it) .ge. iichemcol(1,it))
     -   write (6,1001) 'Chemical element:    ',(iichemcol(k,it),k=1,2)
      if (iialtcol(it) .gt. 0) write (6,1003) iialtcol(it)
      if (iiinscol(it) .gt. 0) write (6,1005) iiinscol(it)
      if (iidcol(it) .gt. 0) write (6,1004) iidcol(it)
      return
1000  format(' Record format specification for ',a,' input')
1001  format(' Column range for ',a,i3,' - ',i3)
1002  format(' Number of characters in a record:',11x,i3)
1003  format(' Alternate marker column:',20x,i3)
1004  format(' Deletion  marker column:',20x,i3)
1005  format(' Insertion marker column:',20x,i3)
      end
