      subroutine read_write_ccc(ndials,iout,irw)
      parameter (MAXPHI=400,MAX2D=5000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-(MAX2D*MAX2D+MAX2D))
      common /nnwork/ ccc(MAX2D,MAX2D),itemp(MAX2D),fill(IFILL2)
      if (irw .eq. 1) then
        write (iout) ndials
        do i=1,ndials
          write (iout) (ccc(i,j),j=1,ndials)
        end do
      else
        rewind iout
        read (iout) ndials
        do i=1,ndials
          read (iout) (ccc(i,j),j=1,ndials)
        end do
      end if
      return
      end
