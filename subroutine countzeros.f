      subroutine countzeros(z,n,nzero)
      dimension z(n)
      nzero=0
      do i=1,n
        if (z(i) .eq. 0.0) nzero=nzero+1
      end do
      return
      end
