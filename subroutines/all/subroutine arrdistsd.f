      subroutine arrdistsd(a1,a2,dist,dist2)
c#    MMC routine 277 lstmod: 03/29/00
c*****dist2 = sum(ai-bi)**2
      real*8 a1,a2,dist,dist2
      real*8 x,y,z
      dimension a1(3),a2(3)
      x=a1(1)-a2(1)
      y=a1(2)-a2(2)
      z=a1(3)-a2(3)
      dist2=x*x+y*y+z*z
      dist=dsqrt(dist2)
      return
      end
