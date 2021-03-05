      subroutine findlim(index,ifst,ilst,max)
      dimension index(max)
      i=max
      do while (i .gt. 0 .and. index(i) .eq. 0)
        i=i-1
      end do
      ilst=i
      if (ilst .eq. 0) return
      i=1
      do while (i .lt. ilst .and. index(i) .eq. 0)
        i=i+1
      end do
      ifst=i
      return
      end
