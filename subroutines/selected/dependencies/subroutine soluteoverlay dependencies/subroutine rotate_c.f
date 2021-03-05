      subroutine rotate_c(c,n,rot,cnew,lab,llab)
      dimension c(3,n),rot(3,3),cnew(3,n),cn(3)
      character*(*) lab
      iverb=0
      call check_rotmat(rot,lab,llab,ifail,iverb)
c     Orient the molecule
      do im=1,n
        cn(1)=rot(1,1)*c(1,im)+rot(1,2)*c(2,im)+rot(1,3)*c(3,im)
        cn(2)=rot(2,1)*c(1,im)+rot(2,2)*c(2,im)+rot(2,3)*c(3,im)
        cn(3)=rot(3,1)*c(1,im)+rot(3,2)*c(2,im)+rot(3,3)*c(3,im)
        call trnsfr(cnew(1,im),cn,3)
c       cnew(1,im)=cn1
c       cnew(2,im)=cn2
c       cnew(3,im)=cn3
      end do
      return
      end
