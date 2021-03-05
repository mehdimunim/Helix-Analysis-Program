      subroutine isambertraj(filename,namlen,itrajtyp)
c     itrajtyp: 0: not aan Amber trajectory file
c               1: likely to be an Amber tarjectory file
      character*200 filename
      character* 132 line
      call openfile(97,0,' ',1,'old',filename,namlen,notfound,0,1,1,1,0)
      rewind 97
      itrajtyp=0
      read (97,1000,end=99,err=99) line
      call blankout(line,1,132)
      read (97,1000,end=99,err=99) line
      call lastchar(line,ifc,132)
      if (ifc .le. 72 .or. ifc .gt. 80) go to 99
      ndec=0
      do ic=1,ifc
       if (idigit(line(ic:ic),2) .eq. 0) go to 99
       if (line(ic:ic) .eq. '.') ndec=ndec+1
      end do
      if (ndec .ne. 10) go to 99
      itrajtyp=1
99    close (97)
      return
1000  format(a)
      end
