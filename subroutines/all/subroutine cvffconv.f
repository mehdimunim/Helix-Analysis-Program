      subroutine cvffconv(cvffname,charmmname,incvff,icvff,outnam,n,
     -  iused,ntypein,iambig,lstambig,nambig,noconv,ioutmsg,maxtyp)
      character*4 incvff,outnam,xxxx
      dimension incvff(n),icvff(n),outnam(n),iused(maxtyp),lstambig(n),
     -  iambig(maxtyp)
      character*2 cvffname(100)
      character*4 charmmname(100)
      data xxxx /'****'/
c     Convert cvff types into Charmm22 type (first approximation)
c      write (6,8888) cvffname
c8888  format(20(1x,a2))
      call zeroiti(iused,0,maxtyp)
      nambig=0
      noconv=0
      do i=1,n
        icvff(i)=0
        j=0
        do while (j .lt. ntypein .and. icvff(i) .eq. 0)
          j=j+1
          if (incvff(i)(1:2) .eq. cvffname(j)) then
            icvff(i)=j
            outnam(i)=charmmname(j)
            if (outnam(i)(1:1) .eq. '?') then
              write (6,101) incvff(i)(1:2),i
              if (ioutmsg .gt. 0) write (ioutmsg,101) incvff(i)(1:2),i
              noconv=noconv+1
            else
              iused(j)=iused(j)+1
            end if
            if (iambig(j) .gt. 0) then
              nambig=nambig+1
              lstambig(nambig)=i
            end if
          end if
        end do
        if (icvff(i) .eq. 0) then
          write (6,100) incvff(i)(1:2),i
          if (ioutmsg .gt. 0) write (ioutmsg,100) incvff(i)(1:2),i
          outnam(i)=xxxx
        end if
      end do
      return
100   format(' ERROR: atom type ',a2,' atom number',i4,' is not a',
     -  ' valid cvff type',/,' - check input file')
101   format(' PROBLEM: no conversion is given for atom type ',a2,
     -  ' atom number ',i4)
      end
