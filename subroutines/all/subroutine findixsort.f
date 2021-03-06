      subroutine findixsort(list,ifrst,nlist,intg,ix,itry)
      dimension list(nlist)
c     Find intg in the sorted list
      imin=ifrst
      imax=nlist
      ix=0
      if (imin .gt. imax) then
        itry=imin
        return
      end if
      do while (imax-imin .gt. 1)
        itry=(imax+imin)/2
        if (list(itry) .eq. intg) then
          ix=itry
          return
        else if (list(itry) .gt. intg) then
          imax=itry
        else
          imin=itry
        end if
      end do
      itry=imax
      if (imax .eq. imin) then
        if (list(itry) .eq. intg) ix=itry
      else
        if (list(imin).eq. intg) then
          ix=imin
        else if (list(imax).eq. intg) then
          ix=imax
        end if
      end if
      return
      end
