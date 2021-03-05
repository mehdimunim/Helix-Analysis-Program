      subroutine rgbrainbowcolor(iout,frac)
c#    MMC routine 355 lstmod: 12/05/03
c*****Set the PS color to frac way on the rainbow scale (0<=frac<=1)
      if (frac .lt. 0.0) then
        red=0.0
        green=0.0
        blue=0.0
      else if (frac .le. 0.5) then
        red=2.0*(0.5-frac)
        green=2.0*frac
        blue=0.0
      else
        red=0.0
        green=2.0*(0.5-(frac-0.5))
        blue=2.0*(frac-0.5)
      end if
      write (iout,1000) red,green,blue
      return
1000  format(3f5.1,' setrgbcolor')
      end
