      subroutine lpshiftc(c,list,nlist,del,c0,np,maxrec)
      dimension c(3,maxrec),list(maxrec),del(3),c0(3)
c     Generate linear or planar shift
      do k=1,3
        del(k)=0.0
        do ip=1,np
          del(k)=del(k)+c(k,list(nlist+ip))
        end do
        del(k)=del(k)/np-c0(k)
      end do
      call norm(del,1.0)
      return
      end
