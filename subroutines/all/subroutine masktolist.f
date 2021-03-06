      subroutine masktolist(index,mask,n,nfinal,isave)
      dimension index(n),mask(n)
c     Create a list of indices in index where the value of mask is isave
      nfinal=0
      do i=1,n
        if (mask(i) .eq. isave) then
          nfinal=nfinal+1
          index(nfinal)=i
        end if
      end do
      return
      end
