      subroutine rrgbcolor(iout,icol,imax,iprt)
c     Set the color on a continous 0-1 scale
      if (iprt .eq. 1) write (iout,*) 'stroke'
      if (imax .eq. 1) then
        write (iout,1000) 0.0,0.0,0.0
      else
        rcol=float(icol-1)/float(imax-1)
        if (rcol .le. 0.25) then
          write (iout,1000) 1.0,4.0*rcol,0.0
        else if (rcol .le. 0.5) then
          write (iout,1000) (2.0-4.0*rcol),1.0,0.0
        else if (rcol .le. 0.75) then
          write (iout,1000) 0.0,1.0,4.0*(rcol-0.5)
        else
          write (iout,1000) 0.0,4.0-4.0*rcol,1.0
        end if
c       rcol=2.0*float(icol-1)/float(imax-1)
c       if (rcol .le. 1.0) then
c         write (iout,1000) (1.0-rcol),rcol,0.0
c       else
c         rcol=rcol-1.0
c         write (iout,1000) 0.0,(1.0-rcol),rcol
c       end if
      end if
      if (iprt .eq. 1) write (iout,*) 'newpath'
      return
1000  format(3f6.2,' setrgbcolor')
      end
