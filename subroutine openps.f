      subroutine openps(iout,xm,ym,title,ltitle,title1,ltitle1,
     -  filename,lfilename,filename2,lfilename2,npspages,ipspage)
      character*(*) title,title1,filename,filename2
c     print *,'OPENPS'
      ixm=xm+25.0
      iym=ym+25.0
      call psheader(iout,title,ltitle,-10,-10,ixm,iym,npspages,ipspage)
      write (iout,2000)
      call plothead(iout,xm,ym,title,ltitle,title1,ltitle1,
     -  filename,lfilename,filename2,lfilename2)
      return
2000  format('/m { moveto } def',/,'/l { lineto } def',/,
     -  '/np { newpath } def',/, '/sk { stroke } def',/,
     -  '/f { fill } def',/,'/lw { setlinewidth } def',/,
     -  '/Helvetica findfont',/,'12 scalefont',/,'setfont',/,
     -  '25 25 translate')
      end
