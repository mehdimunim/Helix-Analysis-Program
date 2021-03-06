      subroutine getnextcsv(line,int,real,iorr,ic,len)
      character*(*) line
      if (line(ic:ic) .eq. ',') ic=ic+1
      ic1=ic
      ic2=ic1
      real=0.0
      int=0
      call findnextchar(',',line,ic,len)
      if (line(ic:ic).eq. ',') then
        ic2=ic-1
c       print *,'GETNEXTR ic1,ic2=',ic1,ic2
        if (iorr .eq. 0) then
          read (line(ic1:ic2),*,end=999,err=999) int
        else
          read (line(ic1:ic2),*,end=999,err=999) real
        end if
      else
c       Last value, without comma
        ic=ic1
        call nextblank(line,ic,len)
        ic2=ic-1
        if (iorr .eq. 0) then
          read (line(ic1:ic2),*,end=999,err=999) int
        else
          read (line(ic1:ic2),*,end=999,err=999) real
        end if
      end if
      return
999   print *,'Invalid number:',line(ic1:ic2),'| ic1,ic2=',ic1,ic2
      stop
      end
