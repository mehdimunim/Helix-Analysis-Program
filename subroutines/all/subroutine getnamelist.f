      subroutine getnamelist(list,lench,nlist,label,lablen,maxlen)
      character*(*) list(maxlen),label
      character*80 line
      character*200 outfiletmp
      common /logging/ logfile,ipredict
      nlist=0
      line(1:3)='***'
      write (6,1004) lench
      do while (line(1:3) .ne. 'XYZ' .and. line(1:3) .ne. '   ')
        write (6,1000) label(1:lablen),nlist+1
        call blankout(line,1,80)
        read (5,1001) line
        call lastchar(line,iclast,80)
        if (logfile .gt. 0) write (logfile,1001) line(1:iclast)
        if (line(1:4) .eq. 'READ') then
c         Switch to reading from a file
          listinp=55
          namlenl=-1
          call openfile(listinp,1,'name list',9,'old',outfiletmp,
     -      namlenl,notfnd,0,1,1,0,0)
          do while (.true.)
            call blankout(line,1,80)
            read (listinp,1001,end=99) line
            if (line(1:3) .eq. 'XYZ') go to 99
            call lastchar(line,iclast,80)
            if (iclast .gt. lench) write (6,1002) lench
            nlist=nlist+1
            call blankout(list(nlist),1,lench)
            list(nlist)(1:iclast)=line(1:iclast)
          end do
99        close (listinp)
          line(1:3)='XYZ'
        else
          if (iclast .gt. lench) write (6,1002) lench
        end if
        if (line(1:3) .ne. 'XYZ') then
          if (nlist .eq. maxlen) then
            write (6,1003) maxlen
            call askstop(0)
            return
          else
            nlist=nlist+1
            call blankout(list(nlist),1,lench)
            list(nlist)(1:iclast)=line(1:iclast)
          end if
        end if
      end do
      write (6,5655) (label(1:lablen),i,list(i)(1:lench),i=1,nlist)
5655  format(1x,a,' #',i4,':',a)
      return
1000  format(1x,a,' #',i3,' (XYZ or enter to end list)=',$)
1001  format(a)
1002  format(' WARNING: Only the first ',i2,' characters will be used')
1003  format(' WARNING: maximum number of list elements (',i4,') has ',
     -  'been reached - input stops')
1004  format(' Type READ as the residue name to read a list from a ',
     -  'file',/,' Names are of maximum ',i2,' characters')
      end
