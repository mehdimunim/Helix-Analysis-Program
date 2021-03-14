      subroutine getname(name,len,q,lqin,maxlen,defname,ld,nonameok,
     -  ihelp,itip)
      character*(*) q,defname,name
      common /logging/ logfile,ipredict
      character*80 qline,ansline,defn
100   lq=lqin
      defn(1:ld)=defname(1:ld)
      call blankout(name,1,maxlen)
      qline(1:lq)=q(1:lq)
      call blankout(ansline,1,80)
      nitems=0
      if (ld+ihelp+itip .gt. 0) then
        qline(lq+1:lq+2)=' ['
        lq=lq+2
        if (ld .gt. 0) then
          qline(lq+1:lq+ld)=defn(1:ld)
          lq=lq+ld
          nitems=nitems+1
        end if
        if (ihelp .gt. 0) then
          if (nitems .gt. 0) then
            qline(lq+1:lq+1)=','
            lq=lq+1
          end if
          qline(lq+1:lq+1)='?'
          lq=lq+1
          nitems=nitems+1
        end if
        if (itip .gt. 0) then
          if (nitems .gt. 0) then
            qline(lq+1:lq+1)=','
            lq=lq+1
          end if
          qline(lq+1:lq+1)='+'
          lq=lq+1
        end if
        qline(lq+1:lq+1)=']'
        lq=lq+1
      end if
      write (6,1000) qline(1:lq)
      read (5,1001) ansline
      call lastchar(ansline,lc,80)
      if (ansline(1:1) .eq. '?' .and. lc .eq. 1) then
        call explanation(ihelp,0)
        go to 100
      end if
      if (ansline(1:1) .eq. '+' .and. lc .eq. 1) then
        call explanation(0,itip)
        go to 100
      end if
      ifail=0
      if (ansline(1:1) .eq. '"') then
        call readquotestring(ansline,'"',i1,i2,ifail,80)
        if (ifail .gt. 0) go to 100
      else if (ansline(1:1) .eq. "'") then
        call readquotestring(ansline,"'",i1,i2,ifail,80)
        if (ifail .gt. 0) go to 100
      else
        ii=1
        call nextchar(ansline,ii,80)
        i1=ii
        call nextblank(ansline,ii,80)
        i2=ii-1
      end if
      if (i2-i1+1 .gt. maxlen) then
        write (6,1002) ansline(i1:i2),maxlen
        go to 100
      else if (i1 .gt. i2) then
        if (ld .gt. 0) then
c         Use default
          name(1:ld)=defn(1:ld)
          len=ld
        else if (nonameok .eq. 0) then
          print *,'Please, type in a name'
          go to 100
        else
          len=0
        end if
      else
        len=i2-i1+1
        name(1:len)=ansline(i1:i2)
      end if
      if (logfile .gt. 0) then
        nbl=0
        do ic=i1,i2
          if (ansline(ic:ic) .eq. ' ') nbl=nbl+1
        end do
        if (nbl .eq. 0) write (logfile,1001) ansline(i1:i2)
        if (nbl .gt. 0) write (logfile,1003) ansline(i1:i2)
      end if
      return
1000  format(1x,a,'=',$)
1001  format(a)
1002  format(' Name read (',a,') is longer than allowed (',i3,')')
1003  format('"',a,'"')
      end
