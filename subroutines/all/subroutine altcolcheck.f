      subroutine altcolcheck(lineinp,idcol,iruntyp,ialtcol,ikeepfullalt,
     -  altnam,naltnam,naltrec,naltdel,ipredict,froccin,outfile,altfile,
     -  namleno,namlena,ncol,n,asterisk,maxrec)
      dimension froccin(maxrec)
      character*1 altnam(maxrec),asterisk
      character*132 lineinp
      character*200 outfile,altfile
c     print *,'ALTCOLCHECK n=',n,' ialtcol=',ialtcol
      if (lineinp(idcol:idcol) .ne. asterisk .and. iruntyp .ne. 21) then
        if (lineinp(ialtcol:ialtcol) .ne. ' ') then
c         Alternative location found
          if (froccin(n) .eq. 1.0) then
            if (ikeepfullalt .eq. -1) then
              print *,'Fully occupied alternative record found:'
              print *,lineinp
              call askyn(
     -        'Do you want to keep such alternative record mark',
     -        48,1,-1,ikeepfullalt,0,0)
            end if
            if (ikeepfullalt .eq. 0) lineinp(ialtcol:ialtcol)=' '
          end if
        end if
        if (lineinp(ialtcol:ialtcol) .ne. ' ') then
          if (naltnam .eq. 0) then
            naltnam=1
            naltrec=0
            altnam(naltnam)=lineinp(ialtcol:ialtcol)
            if (ipredict .eq. 0) then
              call changeext(outfile,altfile,namleno,namlena,
     -          'alt',3,1,0)
              call openfile(40,0,' ',1,'new',altfile,namlena,
     -          notfnd,0,1,1,0,0)
              write (6,2040) altfile(1:namlena)
              write (40,2044)
            end if
            iname=1
          else
c           See if name already occurred
            do iname=1,naltnam
              if (altnam(iname) .eq. lineinp(ialtcol:ialtcol))
     -          go to 22
            end do
            naltnam=naltnam+1
            iname=naltnam
            altnam(naltnam)=lineinp(ialtcol:ialtcol)
          end if
22        if (ipredict .eq. 0) then
            call writeline(40,lineinp,1,ncol,0)
            naltrec=naltrec+1
          else
c           Duplicates will be dropped
            if (iname .ne. 1) then
              naltdel=naltdel+1
              lineinp(idcol:idcol)=asterisk
            end if
          end if
        end if
      end if
      return
2040  format(' The input PDB file contains records with non-blank ',
     -  'character in column 17.',/,
     -  6x,'These records, marking alternate positions, will be ',
     -  'written on file',/,6x,a)
2044  format('REMARK RECORDS WITH ALTERNATE LOCATION MARKER')
      end
