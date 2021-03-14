      subroutine print_auc(lab,llab,ntracks,nframe,bondname,lbondname)
      character*(*) lab,bondname(ntracks)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      parameter (MAXFRAMES2=MAXFRAMES/2)
      parameter (MAXFRAMES21=MAXFRAMES2+1)
      common /aucw/ auc_all(MAXFRAMES21,MAXCOPY)
      call openfile(91,0,' ',1,'new',lab,llab,notfnd,0,1,1,0,0)
      ntrackw=min0(ntracks,MAXCOPY)
      write (91,1003) (it,it=1,ntrackw)
      write (91,1004) (it,bondname(it)(1:lbondname),it=1,ntrackw)
      write (91,1001) 0,(1.0,it=1,ntrackw)
      do if=1,nframe/2
        write (91,1001) if,(auc_all(if,it),it=1,ntrackw)
      end do
      write (91,1002) (auc_all(nframe/2+1,it),it=1,ntrackw)
      close (91)
      return
1001  format(i7,400f7.3)
1002  format('#AUCSUM',400f7.1)
1003  format('#      ',400i7)
1004  format('#      ',400(i4,':',a))
      end
