      subroutine checkdir(dirname,ldirname,iunit,iopen)
      character*(*) dirname
      character*200 filename
      filename(1:ldirname)=dirname(1:ldirname)
      filename(ldirname+1:ldirname+7)='/xyzXYZ'
      open(unit=iunit,status='new',file=filename(1:ldirname+7),
     -  iostat=iopen,form='formatted')
      if (iopen .eq. 0) close(iunit,status='delete')
      return
      end
