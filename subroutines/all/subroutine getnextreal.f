      subroutine getnextreal(line,r,ic,len)
      character*(*) line
      call nextchar(line,ic,len)
      ic1=ic
      call nextblank(line,ic,len)
      ic2=ic
c     print *,'GETNEXTR ic1,ic2=',ic1,ic2
      read (line(ic1:ic2-1),*,err=999) r
      return
999   print *,'Invalid real:',line(ic1:ic2-1)
      stop
      end
