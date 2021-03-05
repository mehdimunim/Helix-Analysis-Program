      subroutine printrange(line,i1,i2,ic,incr,iout)
      character*(*) line
      if (i1 .eq. i2) then
c       Singleton
        if (ic .gt. 72-incr) then
          write (iout,1000) line(1:ic-1)
          ic=1
        end if
        call writeint(line,ic,i1,lenw)
      else
c       Range
        if (ic .gt. 68-incr) then
          write (iout,1000) line(1:ic-1)
          ic=1
        end if
        call writeint(line,ic,i1,lenw)
        line(ic:ic)='-'
        ic=ic+1
        call writeint(line,ic,i2,lenw)
      end if
      line(ic:ic)=','
      ic=ic+1
      return
1000  format(1x,a)
      end
