      subroutine getint(q,len,idef,noneg,maxval,in,ihelp)
      character*(*) q
      character*132 ansline,pline
      common /logging/ logfile,ipredict
c     print *,'GETINT idef,noneg,maxval=',idef,noneg,maxval
      lenq=len
      pline(1:1)=' '
      pline(2:lenq+1)=q(1:lenq)
c     Allow idef == in
      idefault=idef
      if (idefault .ne. 999999) then
        if (maxval .gt. 0 .and. idefault .gt. maxval) then
          write (6,2002)  idefault,maxval
          lenq=lenq+1
          pline(lenq:lenq)='='
        else
1003  format(i5)
1004  format(i4)
1005  format(i3)
1006  format(i2)
c         Put default on query
          pline(lenq+2:lenq+3)=' ['
          if (iabs(idefault) .ge. 10000) then
            write (pline(lenq+4:lenq+11),1002) iabs(idefault)
            lenq=lenq+11
          else if (iabs(idefault) .ge. 1000) then
            write (pline(lenq+4:lenq+8),1003) iabs(idefault)
            lenq=lenq+8
          else if (iabs(idefault) .ge. 100) then
            write (pline(lenq+4:lenq+7),1004) iabs(idefault)
            lenq=lenq+7
          else if (iabs(idefault) .ge. 10) then
            write (pline(lenq+4:lenq+6),1005) iabs(idefault)
            lenq=lenq+6
          else
            write (pline(lenq+4:lenq+5),1006) iabs(idefault)
            lenq=lenq+5
          end if
          if (ihelp .eq. 0) then
            pline(lenq+1:lenq+2)=']='
            lenq=lenq+2
          else
            pline(lenq+1:lenq+4)=',?]='
            lenq=lenq+4
          end if
        end if
      else
        if (ihelp .eq. 0) then
          pline(lenq+2:lenq+2)='='
          lenq=lenq+2
        else
          pline(lenq+2:lenq+5)='[?]='
          lenq=lenq+5
        end if
      end if
100   call writeline(6,pline,1,lenq,1)
      call blankout(ansline,1,132)
      read (5,1001,end=99) ansline
      ii=1
      call nextchar(ansline,ii,132)
      if (ansline(ii:ii) .eq. '?') then
        if (ihelp .eq. 0) then
          print *,'Sorry, no help for this'
        else
          call explanation(ihelp,0)
        end if
        go to 100
      end if
      i1=ii
      call nextblank(ansline,ii,132)
      i2=ii-1
      if (i1 .gt. i2) then
        if (idefault .eq. 999999) go to 100
        in=idefault
        if (logfile .gt. 0) write (logfile,3000)
      else
        ic1=i1
        if (ansline(i1:i1) .eq. '-' .or.  ansline(i1:i1) .eq. '+')
     -    ic1=i1+1
        do i=ic1,i2
          if (idigit(ansline(i:i),1) .ne. 1) go to 999
        end do
        read (ansline(i1:i2),*,err=998) in
        if (logfile .gt. 0) write (logfile,*) in
      end if
      if (noneg .eq. 1 .and. in .lt. 0) then
        write (6,2000) q(1:len)
        go to 100
      end if
      if (maxval .gt. 0 .and. in .gt. maxval) then
        write (6,2001) maxval
        go to 100
      end if
      return
99    print *,'ERROR: run out of data'
      stop
998   print *,'Invalid input for an integer:',ansline(i1:i2)
      print *,'I1,I2=',i1,i2
      go to 900
999   print *,'Invalid character for an integer:',ansline(i:i)
900   call lastchar(ansline,lc,132)
      if (ishexadecimal(ansline,lc) .eq. 1) then
        print *,'Looks like the input is hexadecimal'
        call askyn('Do you want to interpret it as hexadecimal',42,1,-1,
     -    ihx,0,0)
        if (ihx .eq. 1) then
          in=iconvhex(ansline,lc)
          print *,'Input string ',ansline(1:lc),' converted to ',in
          return
        end if
      end if
      go to 100

1001  format(a132)
1002  format(i8)
2000  format(' ERROR: negative number is invalid for ',a)
2001  format(' ERROR: number read exceeds limit (',i9,')')
2002  format(' PROGRAM ERROR: default value (',i8,') exceeds maximum ',
     -  'value (',i8,')')
3000  format(a)
      end
