      subroutine readxcluster(natoms,c,inpt,nconfig)
      dimension c(3,natoms)
      character*132 inpline
c     X-cluster new coordinates
c     print *,'natoms,inpt,nconfig=',natoms,inpt,nconfig
      read (inpt,1001) inpline
c     print *,'inpline=',inpline
      read (inpline(1:6),1008,err=999) nn
      if (-nn .ne. natoms) then
        print *,'ERROR: Number of atoms in the first ',
     -    'structure=',natoms
        print *,'Number of atoms in the X-cluster structure=',-nn
        print *,'Nconfig=',nconfig
        stop
      end if
      do i=1,natoms
        read (inpt,1001) inpline
        read (inpline(6:41),1103) (c(k,i),k=1,3)
      end do
      return
999   print *,'Invalid number of atoms reading in xcluster'
      stop
1001  format(a132)
1008  format(i6)
1103  format(3f12.5)
      end
