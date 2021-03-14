      subroutine getintline(q,len,noneg,maxval,inputs,n,ihelp)
      dimension inputs(n)
      character*(*) q
      character*132 ansline,pline
      common /logging/ logfile,ipredict
      lenq=len
      pline(1:1)=' '
      pline(2:lenq+1)=q(1:lenq)
      if (ihelp .eq. 0) then
        pline(lenq+2:lenq+2)='='
        lenq=lenq+2
      else
        pline(lenq+2:lenq+5)='[?]='
        lenq=lenq+5
      end if
100   call writeline(6,pline,1,lenq,1)
      call blankout(ansline,1,132)
      read (5,1001) ansline
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
      call lastchar(ansline,ilc,132)
      if (logfile .gt. 0) write (logfile,3000) ansline(1:ilc)
      read (ansline,*,end=888,err=999) (inputs(k),k=1,n)
      do k=1,n
        if (inputs(k) .gt. maxval) then
          write (6,2001) inputs(k),maxval
          go to 100
        else if (noneg .eq. 1 .and. inputs(k) .lt. 0) then
          write (6,2000) inputs(k),q(1:len)
          go to 100
        end if
      end do
      return
888   print *,'Insufficient number of entries'
      go to 100
999   print *,'Invalid input for an integer:',ansline(1:ilc)
      go to 100
1001  format(a132)
2000  format(' ERROR: negative number read (',i9,') is invalid for ',a)
2001  format(' ERROR: number read (',i9,') exceeds limit (',i9,')')
3000  format(a)
      end
