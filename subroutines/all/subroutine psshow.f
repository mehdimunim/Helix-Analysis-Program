      subroutine psshow(iout,string,len)
      character*(*) string
      character*200 clean
      inc=0
      do ic=1,len
        if (string(ic:ic) .eq. '(') then
          clean(ic+inc:ic+inc)='\\'
          clean(ic+inc+1:ic+inc+1)='('
          inc=inc+1
        else if (string(ic:ic) .eq. ')') then
          clean(ic+inc:ic+inc)='\\'
          clean(ic+inc+1:ic+inc+1)=')'
          inc=inc+1
        else
          clean(ic+inc:ic+inc)=string(ic:ic)
        end if
      end do
      write (iout,1000) clean(1:len+inc)
      return
1000  format('(',a,' ) show')
      end
