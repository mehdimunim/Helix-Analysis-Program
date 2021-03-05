      subroutine getnextint(line,int,ic,len)
      character*(*) line
      call nextchar(line,ic,len)
      ic1=ic
      call nextblank(line,ic,len)
      ic2=ic
      read (line(ic1:ic2-1),*,err=999) int
      return
999   print *,'Invalid integer:',line(ic1:ic2-1)
      stop
      end
