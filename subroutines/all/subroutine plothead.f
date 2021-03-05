      subroutine plothead(iout,xm,ym,title,ltitle,title1,ltitle1,
     -  filename,lfilename,filename2,lfilename2)
      character*(*) title,title1,filename,filename2
      yinc=0.0
      if (lfilename2 .gt. 0) then
        write (iout,1000) xm*0.04,ym*0.91+yinc
        write (iout,1003) 'Second file',filename2(1:lfilename2)
        yinc=yinc+ym*0.02
      end if
      if (lfilename .gt. 0) then
        write (iout,1000) xm*0.04,ym*0.91+yinc
        write (iout,1003) 'File',filename(1:lfilename)
        yinc=yinc+ym*0.02
      end if
      if (ltitle1 .gt. 0) then
        write (iout,1000) xm*0.04,ym*0.91+yinc
        call psshow(iout,title1,ltitle1)
        yinc=yinc+ym*0.02
      end if
      write (iout,1000) xm*0.04,ym*0.91+yinc
      itw=0
      if (ltitle .ge. 4) then
        if (title(1:4) .eq. '@#$%') then
          call psshow(iout,'Simulaid-generated plot',23)
          itw=1
        end if
      end if
      if (itw .eq. 0 .and. ltitle .ge. 76) then
        call lastchar(title,lct,ltitle)
        call psshow(iout,title,lct)
      end if
      return
1000  format(2f12.5,' m')
1003  format('(',a,' analyzed: ',a,' ) show')
      end
