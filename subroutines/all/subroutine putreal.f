      subroutine putreal(out,len,val,ndec)
      character*(*) out
      character*25 out25
c     Puts a real number into out
      if (ndec .eq. 0) then
        write (out25,100) val
      else if (ndec .eq. 1) then
        write (out25,101) val
      else if (ndec .eq. 2) then
        write (out25,102) val
      else if (ndec .eq. 3) then
        write (out25,103) val
      else if (ndec .eq. 4) then
        write (out25,104) val
      else if (ndec .eq. 5) then
        write (out25,105) val
      else if (ndec .eq. 6) then
        write (out25,106) val
      else if (ndec .eq. 7) then
        write (out25,107) val
      else if (ndec .eq. 8) then
        write (out25,108) val
      else if (ndec .eq. 9) then
        write (out25,109) val
      else if (ndec .eq. 10) then
        write (out25,110) val
      else
        print *,'Illegal number of decimal digits in putmol:',ndec
        write (out25,110) val
      end if
c     Check for exceeding the length
      ic=1
      call nextchar(out25,ic,25)
      if (ic .le. 25-len) then
        out=out25(ic:ic+len-1)
      else
        out=out25(25-len+1:25)
      end if
c     print *,'PUTREAL val=',val,' out=',out,' len=',len
      return
100   format(f25.0)
101   format(f25.1)
102   format(f25.2)
103   format(f25.3)
104   format(f25.4)
105   format(f25.5)
106   format(f25.6)
107   format(f25.7)
108   format(f25.8)
109   format(f25.9)
110   format(f25.10)
      end
