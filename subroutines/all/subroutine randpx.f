      subroutine randpx(n,rno)
c#    MMC routine 429 lstmod: 10/04/86
c*****Congruential random number generator, Forsythe's constants
      integer*4 ixo,iy
      common /rangen/ ixo
      dimension rno(n)
      do i=1,n
        iy=ixo*314159269+453806245
c       Eliminate bits over 31
        iy=ibclr(iy,31)
        rno(i)=float(iy)/2.1474836E+09
        ixo=iy
      end do
      return
      end
