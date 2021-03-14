      subroutine indexit(index,ifirst,ilast,incr)
      dimension index(ilast)
      do i=ifirst,ilast
        index(i)=i+incr
      end do
      return
      end
