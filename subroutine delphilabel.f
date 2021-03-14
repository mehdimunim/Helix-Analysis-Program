      subroutine delphilabel(c,n,nslt,xstart,ystart,zstart,gx,gy,gz,cv)
      dimension c(3,n),cv(n)
      parameter (MAXPHI=400)
      common /nnwork/ phimap(MAXPHI,MAXPHI,MAXPHI)
      nl=0
      do ia=1,nslt
        call interpolate(c(1,ia),c(2,ia),c(3,ia),gx,gy,gz,
     -    xstart,ystart,zstart,phi)
        cv(ia)=phi
      end do
      return
      end
