      subroutine pbcreset(c,n,crep,cell,ncell,ixyzexcld,ixyzincld,img)
      dimension c(3,n),crep(3),cell(3,ncell)
c     Finds the PBC cell of the representative atom and shifts the whole molec
c     back to the central cell
c     print *,'PBCRESET ncell,ixyzexcld,ixyzincld=',
c    -  ncell,ixyzexcld,ixyzincld
      if (n .eq. 0) return
      call genimdist123dim(crep,cell,1,ncell,ixyzexcld,ixyzincld,
     -  img,rmin2)
      do ia=1,n
        call arrdiff(c(1,ia),cell(1,img),c(1,ia),3)
      end do
      return
      end
