      subroutine printrmsd2(iout,n1,n2,lab,llab)
      character*(*) lab
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-MAX2D*MAX2D)
      common /nnwork/ rmsd2d(MAX2D,MAX2D),fill(IFILL2)
      write (iout,*) lab(1:llab)
      do i=1,n1
        write (iout,1000) i,(rmsd2d(i,j),j=1,n2)
      end do
      return
1000  format(i4,' rmsd2:',10f8.2)
      end
