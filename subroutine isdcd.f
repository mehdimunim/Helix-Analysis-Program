      subroutine isdcd(filename,namlen,idcdtyp)
c     idcdtyp: 0: not a DCD file; 1: coordinate DCD file; 2: velocity DCD file
      character*200 filename
      character*4 header
      call openfile(97,0,' ',1,'old',filename,namlen,notfound,0,2,1,1,0)
      rewind 97
      idcdtyp=0
      read (97,end=99,err=99) header
      if (header .eq. 'CORD') idcdtyp=1
      if (header .eq. 'VELD') idcdtyp=2
99    close (97)
      return
      end
