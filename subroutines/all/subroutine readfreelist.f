      subroutine readfreelist(list,len,inp,isort,maxlen)
      dimension list(maxlen)
      character*80 line
      len=0
      ic=0
      lc=0
      iprev=0
      nerr=0
      do while (.true.)
        if (ic .ge. lc) then
          call blankout(line,1,80)
          read (inp,1000,end=999) line
          call lastchar(line,lc,80)
          if (lc .eq. 0) return
          ic=1
        end if
        call nextstring(line,ic,ic1,ic2,80)
        ic=ic2+1
        iok=0
        read (line(ic1:ic2),*,err=777) i
        iok=1
        len=len+1
        if (len .gt. maxlen) then
          write (6,1001) maxlen
          return
        end if
        if (isort .eq. 1) then
          if (len .eq. 1 .or. i .gt. iprev) then
            list(len)=i
          else
            write (6,1002) i,iprev
            len=len-1
          end if
          iprev=i
        else
          list(len)=i
        end if
777     if (iok .eq. 0) then
          print *,'ERROR: Invalid list element:',line(ic1:ic2),
     -      ' ignored'
          nerr=nerr+1
          if (nerr .gt. 2*maxlen) then
            print *,'Too many invalid elements - stopping; iuinp=',inp
            stop
          end if
        end if
      end do
999   return
1000  format(a)
1001  format(' ERROR: list length exceeds',i6,' - redimension Simulaid')
1002  format(' ERROR: list element (',i6,') is not greater than the ',
     -  'previous value (',i6,')')
      end
