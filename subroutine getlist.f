      subroutine getlist(list,len,minval,maxval,nlist,maxlen)
      dimension list(nlist,maxlen)
      common /logging/ logfile,ipredict
      character*200 listfile
      character*80 listprompt,line
c     Get a list of integers (nlist=1) or a pair of integers (nlist=2)
c     from a file or from the terminal
      namlenl=-1
      listinp=55
      listprev=-999999
      iallpair=0
      irangeinp=0
      if (nlist .eq. 1) then
        write (6,1006) ' ',' in increasing order'
      else
        write (6,1006)
        call askyn(
     -    'Do you want to input just one number per line',
     -    45,1,-1,iallpair,99,0)
         if (iallpair .eq. 0) write (6,1006) ' one pair per line'
         if (iallpair .eq. 1) write (6,1006) ' one per line'
      end if
      write (6,1011)
      call getname(listfile,llistfile,'Name of the list file',21,80,
     -  '',0,1,0,0)
      if (llistfile .gt. 0) then
        call openfile(listinp,1,'list',4,'old',listfile,llistfile,
     -    notfnd,0,1,1,0,0)
      else
        listinp=5
        if (nlist .eq. 1) write (6,1002)
      end if
      maxvalchk=maxval
      if (maxval .eq. 999999) maxvalchk=0
      len=0
      if (nlist .eq. 1 .and. listinp .ne. 5) then
        call readfreelist(list,len,listinp,1,maxlen)
        go to 8
      else if (nlist .eq. 1) then
c       Input list interactively
100     do while (.true.)
          call blankout(line,1,80)
          call getname(line,lline,'0/r/#:',6,80,'',0,1,0,0)
          if (lline .eq. 0) go to 999
          if (line(1:1) .eq. 'r' .or. line(1:1) .eq. 'R') then
c           Input a range
            irangeinp=1
            call getrange(ifst,999999,ilst,999999,increment,1,'number',
     -        6,maxvalchk,0)
            if (ifst .le. listprev) then
              print *,'Ignoring this range - list has to be in ',
     -          'increasing order'
              ndrop=1
            else
              do ics=ifst,ilst,increment
                len=len+1
                list(1,len)=ics
              end do
              listprev=list(1,len)
              if (list(1,len) .ne. ilst) then
                write (6,1005) list(1,len),ilst
              end if
            end if
          else
            read (line,*,err=99,end=8) ics
            if (ics .eq. 0) then
              list(1,len+1)=0
              go to 8
            else
              len=len+1
              list(1,len)=ics
            end if
          end if
          ndrop=0
          if (ics .le. listprev) then
            print *,'Ignoring',ics,' - sorry, list has to be in ',
     -        'increasing order'
            print *,'LISTPREV=',listprev
            len=len-1
            ndrop=1
          end if
          if (ndrop .eq. 0) then
            if (len .eq. maxlen) then
              print *,'List truncated at the',maxlen,'-th element ',
     -         '- redimension Simulaid for longer list'
              go to 8
            end if
c           if (irangeinp .eq. 1) then
c             len=len+1
c             list(1,len)=ics
c             listprev=ics
c           end if
          end if
          listprev=list(1,len)
          irangeinp=0
        end do
99      print *,'Invalid input - ignored'
        go to 100
      else if (iallpair .eq. 0) then
c       Input pair list
200     do while (.true.)
          call blankout(line,1,80)
          if (listinp .eq. 5) then
            write (listprompt,1008) 'pair'
            call getname(line,lline,listprompt,28,80,'',0,1,0,0)
          else
            read (listinp,1004,end=998) line
          end if
          read (line,*,err=88) i1
          if (i1 .eq. 0) go to 8
          len=len+1
          read (line,*,err=88) (list(k,len),k=1,2)
          if (maxvalchk .gt. 0) then
            do k=1,2
              if (list(k,len) .gt. maxval .or.
     -            list(k,len) .lt. minval) then
                write (6,1010) list(k,len),minval,maxval
                len=len-1
                go to 200
              end if
            end do
          end if
          if (len .eq. maxlen) then
            print *,'List truncated at the',maxlen,'-th element ',
     -       '- redimension Simulaid for longer list'
            go to 8
          end if
        end do
88      print *,'Invalid input - ignored'
        go to 200
      else
c       Create pair list from input list
        lenpair=1
        ix1=1
300     do while (.true.)
          if (listinp .eq. 5) write (6,1008) 'element'
          call blankout(line,1,80)
          if (listinp .eq. 5) then
            write (listprompt,1008) 'element'
            call getname(line,lline,listprompt,31,80,'',0,1,0,0)
          else
            read (listinp,1004,end=999) line
          end if
          read (line,*,err=77) i1
          if (ix1 .eq. 1 .and. i1 .eq. 0) go to 777
          if (maxvalchk .gt. 0 .and.
     -        (i1 .gt. maxval .or. i1 .lt. minval)) then
            write (6,1010) i1,minval,maxval
            go to 300
          end if
          list(ix1,lenpair)=i1
          ix1=ix1+1
          if (ix1 .eq. 3) then
            ix1=1
            lenpair=lenpair+1
          end if
          if (lenpair .eq. maxlen) then
            print *,'List truncated at the',maxlen,'-th element ',
     -       '- redimension Simulaid for longer list'
            go to 777
          end if
        end do
77      print *,'Invalid input - ignored'
        go to 300
777     len=lenpair-1
        go to 8
      end if
999   write (6,1007) line
998   list(1,len+1)=0
8     if (listinp .ne. 5) close(listinp)
      if (len .lt. 100) write (6,1000) ((list(k,i),k=1,nlist),i=1,len)
      write (6,1009) len
      return
1000  format(' The list read:',/,(10i8))
1002  format(' Type 0 to finish the list',/,
     -  ' Type a number to be added to the list',/,
     -  ' Type R or r to be able to enter a range')
1004  format(a)
1005  format(' WARNING: Last number added to the list=',i6,
     -  ' (instead of ',i6,')')
1006  format(' List file is a free-formatted list of sequence numbers',
     -  a,/,a)
1007  format(' Premature end of the list file - last record read=',/,a)
1008  format(' Next ',a,' (or 0 to finish): ',$)
1009  format(' Finished reading',i6,' list items')
1010  format(' ERROR: List element ',i8,' is outside the alllowed ',
     -  'range [',i6,',',i8,']')
1011  format(' Hit return to input list interactively')
      end
