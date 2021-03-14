      subroutine extract(c,index,ndim,n,nfinal)
      dimension c(ndim,n),index(n)
      nfinal=0
      do ia=1,n
        if (index(ia) .eq. 0) then
          nfinal=nfinal+1
          if (nfinal .lt. ia) call trnsfr(c(1,nfinal),c(1,ia),ndim)
        end if
      end do
      return
      end
