      subroutine newdist(r,c,list,c0,n,maxrec)
      dimension r(maxrec),c(3,maxrec),list(maxrec),c0(3)
      do i=1,n
        li=list(i)
        r(i)=dist2(c(1,li),c0)
      end do
      return
      end
