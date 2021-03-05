      subroutine reverseindex(indexrev,indexorig,ifres,ifirst,ilast,max)
      dimension indexrev(max),indexorig(max),ifres(max)
      call zeroiti(indexrev,0,max)
      do i=ifirst,ilast
        indexrev(indexorig(ifres(i)))=i
      end do
      return
      end
