      subroutine clean_ng(nfrst,n,nng,ing,maxn)
      dimension nng(maxn),ing(maxn,maxn)
c     Remove the eliminated nodes from the nn list of the rest
      do i=nfrst,n
        ndel=0
        do j=1,nng(i)
          if (ing(j,i) .le. 0) then
            ndel=ndel+1
          else
            ing(j-ndel,i)=ing(j,i)
          end if
        end do
        nng(i)=nng(i)-ndel
      end do
      return
      end
