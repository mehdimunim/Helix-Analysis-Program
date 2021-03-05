      subroutine rgbcolor_sgi(iout,icol)
c     Set the PS color to the SGI colors (0-7); reverse 0 & 7
c     0:white; 1:red; 2:green; 3:yellow; 4:blue; 5:pink; 6; cyan; 7:orange
c     8:yellow-green 9:back
      if (icol .eq. 0) then
        write (iout,1000) 1.0,1.0,1.0
      else if (icol .eq. 1) then
        write (iout,1000) 1.0,0.0,0.0
      else if (icol .eq. 2) then
        write (iout,1000) 0.0,1.0,0.0
      else if (icol .eq. 3) then
        write (iout,1000) 1.0,1.0,0.0
      else if (icol .eq. 4) then
        write (iout,1000) 0.0,0.0,1.0
      else if (icol .eq. 5) then
        write (iout,1000) 1.0,0.0,1.0
      else if (icol .eq. 6) then
        write (iout,1000) 0.0,1.0,1.0
      else if (icol .eq. 7) then
        write (iout,1000) 1.0,0.5,0.0
      else if (icol .eq. 8) then
        write (iout,1000) 0.5,1.0,0.0
      else if (icol .eq. 9) then
        write (iout,1000) 0.0,0.0,0.0
      end if
      return
1000  format(3f5.1,' setrgbcolor')
      end
