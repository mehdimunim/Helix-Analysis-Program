      subroutine read_head_mae(lab,llab,n,line,len,inpt,iout)
      character*(*) lab,line
      n=0
      do while (n .eq. 0)
        call blankout(line,1,len)
        read (inpt,1000,end=991) line
        ic=1
        call nextchar(line,ic,len)
        if (ic .lt. len) then
          if (line(ic:ic+llab-1) .eq. lab(1:llab) .and.
     -        line(ic+llab:ic+llab) .eq. '[') then
c           Header found
            ic1=ic+llab+1
            ic=ic1
            call findnextchar(']',line,ic,len)
            if (ic .eq. len) go to 992
            read (line(ic1:ic-1),*,end=992,err=992) n
          end if
        end if
      end do
      return
991   write (iout,2001) lab(1:llab)
      return
992   write (iout,2000) 'Invalid number of items',line(1:50)
      return
1000  format(a)
2000  format(' ERROR: ',a,' in line',/,1x,a)
2001  format(' ERROR: run out of data while reading ',a,' header')
      end
