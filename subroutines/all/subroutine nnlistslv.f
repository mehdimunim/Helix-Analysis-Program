      subroutine nnlistslv(nfirst,n,iatnum,c,nneig,nhneig,ineig,maxng)
      dimension nneig(n),nhneig(n),ineig(maxng,n),iatnum(n),c(3,n)
      call zeroiti(nneig,nfirst-1,n)
      do i=nfirst,n
        do j=i+1,n
          r2=dist2(c(1,i),c(1,j))
          call decidebondcut(iatnum(i),iatnum(j),rlim)
          if (r2 .le. rlim .and. iatnum(i)+iatnum(j) .gt. 2) then
            nneig(i)=nneig(i)+1
            nneig(j)=nneig(j)+1
            ineig(nneig(i),i)=j
            ineig(nneig(j),j)=i
            if (iatnum(i) .eq. 1) nhneig(j)=nhneig(j)+1
            if (iatnum(j) .eq. 1) nhneig(i)=nhneig(i)+1
          end if
        end do
      end do
      return
      end
