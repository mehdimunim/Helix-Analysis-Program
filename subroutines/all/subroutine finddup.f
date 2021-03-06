      subroutine finddup(p,ixdup,n,nunique)
      dimension p(3,n),ixdup(n),ixshift(48)
      do i=2,n
        do j=1,i-1
          if (ixdup(j) .eq. j) then
            if (dist2(p(1,i),p(1,j)) .le. 0.001) then
              ixdup(i)=j
              go to 100
            end if
          end if
        end do
100     continue
      end do
      ndel=0
      do i=1,n
        if (ixdup(i) .ne. i) then
c         Duplicate - delete
          ndel=ndel+1
        else
          call trnsfr(p(1,i-ndel),p(1,i),3)
          ixshift(i)=i-ndel
        end if
      end do
      do i=1,n
        ixdup(i)=ixshift(ixdup(i))
      end do
      nunique=n-ndel
c      write (6,7777) n,nunique,(ixdup(i),i=1,n)
c7777  format(' FINDDUP n,nunique=',2i6,(' ixdup=',10i3,3x,10i3))
      return
      end
