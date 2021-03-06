      subroutine altdelcheck(n,naltnam,naltrec,altnam,ipredict,iruntyp,
     -  altcol,ialtcol,idcol,naltdel,ncol,line,index,idrop,linewr,
     -  altfile,namlena,asterisk,maxrec)
      dimension index(maxrec)
      character*1 altnam(50),altcol(maxrec),asterisk
      character*80 linewr
      character*132 line(maxrec)
      character*200 altfile
c     print *,'ALTDELCHECK n=',n,' ialtcol=',ialtcol,
c    -  ' ipredict=',ipredict
      if (naltnam .eq. 0 .or. iruntyp .eq. 21) return
      write (6,2030) naltrec,(altnam(in),in=1,naltnam)
      naltdel=0
      iblankalt=0
      if (ipredict .eq. 1) then
        write (6,2029) altnam(1)
        iblankalt=1
        altnam(1)=' '
      else
        linewr(1:39)='Do you want to drop records with mark  '
        ndrop=0
        do in=1,naltnam
          linewr(39:39)=altnam(in)
          call askyn(linewr,39,1,0,idrop,0,0)
          if (idrop .eq. 0) altnam(in)=' '
          ndrop=ndrop+idrop
        end do
        if (ndrop .eq. naltnam) then
          print *,'WARNING: ALL alternate records will be dropped'
          call askstop(0)
        else
          if (ndrop .eq. 0) then
            print *,'NOTE: ALL alternate records will be kept'
          else
            print *,'List of deleted records will be ADDED to the',
     -        ' file ',altfile(1:namlena)
          end if
        end if
        call askyn(
     -    'Do you want to blank out the alternative mark',45,1,+1,
     -    iblankalt,00,0)
        write (40,2045)
      end if
      do in=1,naltnam
        naltdel_c=0
        do i=1,n
c         Mark locations to be deleted
          if (line(index(i))(ialtcol:ialtcol) .ne. ' ') then
            if (line(index(i))(ialtcol:ialtcol) .eq. altnam(in))then
              naltdel_c=naltdel_c+1
              line(index(i))(idcol:idcol)=asterisk
              if (ipredict .eq. 0)
     -          call writeline(40,line(index(i)),1,ncol,0)
            end if
          end if
        end do
        if (altnam(in) .ne. ' ') then
          write (6,2043) naltdel_c,altnam(in)
          if (naltdel_c .eq. 0) print *,'PROGRAM ERROR: naltdel_c=0'
        end if
        naltdel=naltdel+naltdel_c
      end do
      nab=0
      if (iblankalt .eq. 1) then
        do i=1,n
          if (line(index(i))(ialtcol:ialtcol) .ne. ' ') then
            line(index(i))(ialtcol:ialtcol)=' '
            altcol(i)=' '
            nab=nab+1
          end if
        end do
      end if
      return
2029  format(' Alternative locations ',a,' will be kept and the others',
     -  ' dropped',/,' Alternative marks will be deleted')
2030  format(i4,' alternative location records have been found',/,
     -  ' The following characters were found in column 17:',50(1x,a1))
2043  format(i5,' records containing alternate location marked with ',
     -  a1,' will be deleted')
2045  format('REMARK RECORDS WITH ALTERNATE LOCATION MARKER TO BE ',
     -  'DELETED')
      end
