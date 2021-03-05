      subroutine writepdbd(iout,c,ia,ir,atnam,resnam,segnam,frocc,bfac)
      real*8 c(3)
      character*4 atnam
      character*3 resnam
      character*1 segnam
      write (iout,1000) ia,atnam,resnam,segnam,ir,c,frocc,bfac
      return
1000  format('ATOM  ',i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3,f6.2,f6.1)
      end
