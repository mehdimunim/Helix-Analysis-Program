      subroutine swapl(list,r,nlist,lmax,nlincr,maxrec)
      dimension list(maxrec),r(maxrec)
c     Shift the atom lmax to the end of the list,
c     change nlist by nlincr (-1 or 0)
c     print *,'nlist,lmax=',nlist,lmax
c     print *,'SWAPL nlist,lmax=',nlist,lmax
      call swapi4(list(nlist),list(lmax))
      rr=r(nlist)
      r(nlist)=r(lmax)
      r(lmax)=rr
      nlist=nlist+nlincr
      return
      end
