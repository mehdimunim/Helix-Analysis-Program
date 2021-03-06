      subroutine lookforneig(atnam4,label,ing,ilabel,nfound,ica,ierr)
      character*4 atnam4,label
      if (atnam4 .eq. label) then
        if (ilabel .gt. 0) then
          write (6,1000) ica,label,ilabel,ing
          ierr=ierr+1
        else
          ilabel=ing
          nfound=nfound+1
        end if
      end if
      return
1000  format(' ERROR: CA atom # ',i6,' has extra ',a,' neighbour:',2i7)
      end
