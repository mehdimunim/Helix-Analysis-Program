      subroutine psheader(iout,title,ltitle,ix0,iy0,ixm,iym,npspages,
     -  ipspage)
      character*(*) title
      character*12 today
      common /today_date/ ltoday,today
c     print *,'PSHEADER npspages,ipspage=',npspages,ipspage
      if (npspages .gt. 0) then
        write (iout,1000) title(1:ltitle)
c       if (npspages .lt. 10) then
c         write (iout,1021) npspages
c       else if (npspages .lt. 10) then
c         write (iout,1022) npspages
c       else
c         write (iout,1023) npspages
c       end if
        write (iout,1030) ix0,iy0,ixm,iym
        if (ltoday .gt. 0) write (iout,1001) today(1:ltoday)
        write (iout,1040)
        write (iout,2000)
        ipspage=0
      end if
      return
1000  format('%!PS-Adobe-2.0 ',/,'%%Title: ',a,/,
     -  '%%Creator: SIMULAID')
1001  format('%%CreationDate: ',a)
c1021  format('%%Pages:',i1)
c1022  format('%%Pages:',i2)
c1023  format('%%Pages:',i3)
1030  format('%%BoundingBox:',4i8)
1040  format('%%EndComments')
2000  format('/m { moveto } def',/,'/l { lineto } def',/,
     -  '/np { newpath } def',/, '/sk { stroke } def',/,
     -  '/f { fill } def',/,'/lw { setlinewidth } def',/,
     -  '/r { rlineto } def')
      end
