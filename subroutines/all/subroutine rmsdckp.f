      subroutine rmsdckp(iunit,ireadwrite)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-MAX2D*MAX2D)
      common /nnwork/ rmsd2d(MAX2D,MAX2D),fill(IFILL2)
      if (ireadwrite .eq. 1) then
        write (iunit) rmsd2d
      else
        read (iunit,end=999) rmsd2d
      end if
      return
999   print *,'ERROR: 2D-RMMSD checkpoint file was too short'
      stop
      end
