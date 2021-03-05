      subroutine rgbcolor(iout,icolinp)
c     Set the PS color to the rainbow colors (0-7);
c     0:white; 1:red; 2: orange; 3: yellow; 4:green;
c     5:cyan; 6:blue; 7:violet; 8:indigo; 9:black 10: magenta
      common /colorinfo/ ncolcode,maxcolcode
      dimension icolconv(9,8)
      data icolconv /
     -  1,9,0,0,0,0,0,0,9,
     -  1,6,9,0,0,0,0,0,9,
     -  1,4,6,9,0,0,0,0,9,
     -  1,3,4,6,9,0,0,0,9,
     -  1,3,4,6,10,9,0,0,9,
     -  1,3,4,5,6,10,9,0,9,
     -  1,2,3,4,5,6,10,9,9,
     -  1,2,3,4,5,6,7,8,9/
      if (icolinp .eq. maxcolcode+1 .or. icolinp .eq. 0) then
        icol=icolinp
      else if (icolinp .lt. 0) then
        icol=-icolinp
      else
        icol=icolconv(icolinp,ncolcode)
      end if
      if (icol .eq. 0) then
        write (iout,1000) 1.0,1.0,1.0
      else if (icol .eq. 1) then
        write (iout,1000) 1.0,0.0,0.0
      else if (icol .eq. 2) then
        write (iout,1000) 1.0,0.5,0.0
      else if (icol .eq. 3) then
        write (iout,1000) 1.0,1.0,0.0
      else if (icol .eq. 4) then
        write (iout,1000) 0.0,1.0,0.0
      else if (icol .eq. 5) then
        write (iout,1000) 0.0,1.0,1.0
      else if (icol .eq. 6) then
        write (iout,1000) 0.0,0.0,1.0
      else if (icol .eq. 7) then
        write (iout,1000) 0.5,0.0,1.0
      else if (icol .eq. 8) then
        write (iout,1000) 1.0,0.0,0.5
      else if (icol .eq. 9) then
        write (iout,1000) 0.0,0.0,0.0
      else if (icol .eq. 10) then
        write (iout,1000) 1.0,0.0,1.0
      end if
      return
1000  format(3f5.1,' setrgbcolor')
      end
