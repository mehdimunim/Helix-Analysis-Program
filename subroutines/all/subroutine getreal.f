      subroutine getreal(q,len,default,r,noneg,ihelp)
      character*(*) q
      character*132 ansline,pline
      common /logging/ logfile,ipredict
      lenq=len
      pline(1:1)=' '
      pline(2:lenq+1)=q(1:lenq)
      if (default .ne. 999999.0) then
c       Put default on query
        pline(lenq+2:lenq+3)=' ['
        if (default .eq. 0.0) then
          write (pline(lenq+4:lenq+6),1007) default
          lenq=lenq+6
        else if (abs(default) .ge. 1000.0) then
          write (pline(lenq+4:lenq+11),1002) default
          lenq=lenq+11
        else if (abs(default) .ge. 100.0) then
          write (pline(lenq+4:lenq+10),1003) default
          lenq=lenq+10
        else if (abs(default) .ge. 10.0) then
          write (pline(lenq+4:lenq+9),1004) default
          lenq=lenq+9
        else if (abs(default) .ge. 1.0) then
          write (pline(lenq+4:lenq+8),1005) default
          lenq=lenq+8
        else
          write (pline(lenq+4:lenq+11),1006) default
          lenq=lenq+11
        end if
        if (ihelp .eq. 0) then
          pline(lenq+1:lenq+2)=']='
          lenq=lenq+2
        else
          pline(lenq+1:lenq+4)=',?]='
          lenq=lenq+4
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
        if (default .eq. 999999.0) go to 100
        r=default
        if (logfile .gt. 0) write (logfile,3000)
      else
        read (ansline(i1:i2),*,err=999) r
        if (noneg .eq. 1 .and. r .lt. 0.0) then
          print *,'ERROR: negative number is invalid'
          go to 100
        end if
        if (logfile .gt. 0) write (logfile,*) r
      end if
      return
99    print *,'ERROR: run out of data'
      stop
999   print *,'Invalid input for a real'
      go to 100
1001  format(a132)
1002  format(f8.1)
1003  format(f7.2)
1004  format(f6.2)
1005  format(f5.2)
1006  format(f8.5)
1007  format(f3.1)
3000  format(a)
      end
