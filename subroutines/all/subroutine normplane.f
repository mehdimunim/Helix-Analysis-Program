      subroutine normplane(a,b,c,rn)
      dimension a(3),b(3),c(3),rn(3)
c     rn is the normal to the plane of a,b, and c
      dimension d1(3),d2(3)
      do k=1,3
        d1(k)=a(k)-b(k)
        d2(k)=c(k)-b(k)
      end do
      call vprd(d1,d2,rn)
      call norm(rn,1.0)
      return
      end
