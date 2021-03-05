      subroutine rot90(ix,r)
      dimension r(3,3)
c     Generate rotation matrix for 90 degree rotation around the ix axis
      call unitmat(r)
      iy=mod(ix,3)+1
      iz=mod(ix+1,3)+1
      r(iy,iy)=0.0
      r(iz,iz)=0.0
      r(iy,iz)=1.0
      r(iz,iy)=-1.0
      return
      end
