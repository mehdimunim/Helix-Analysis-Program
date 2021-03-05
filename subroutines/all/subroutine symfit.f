      subroutine symfit(c,co,c1,n,ioppbc,bestrot)
      dimension c(3,n),co(3,n),c1(3,n),bestrot(3,3)
c     Find the closest symmetrical image of the optimized ort to the original
      dimension nrot(3,8),kk(3),rot0(3,3),rot(3,3,3),rotopt(3,3)
      data nrot /4,4,4, 2,2,2, 4,4,4, 6,2,2, 6,2,2, 4,4,4, 4,4,4, 4,4,4/
      pi2=2.0*acos(-1.0)
      d0=dist2n(c,co,n)
      ncha=0
      do k=1,3
        call unitmat(rot(1,1,k))
        k1=mod(k,3)+1
        k2=mod(k+1,3)+1
        angle=pi2/float(nrot(k,ioppbc))
        rot(k1,k1,k)=cos(angle)
        rot(k2,k2,k)=rot(k1,k1,k)
        rot(k1,k2,k)=sin(angle)
        rot(k2,k1,k)=-rot(k1,k2,k)
      end do
      do k1=1,nrot(1,ioppbc)
        kk(1)=k1
        do k2=1,nrot(2,ioppbc)
          kk(2)=k2
          do k3=1,nrot(3,ioppbc)
            kk(3)=k3
            call unitmat(rot0)
            do k=1,3
              do i=1,kk(k)
                call matprod(rot(1,1,k),rot0,rot0)
              end do
            end do
            call rotate_c(c,n,rot0,c1,'SYMFIT1',7)
c           call checkchir(c1,n,2,1,3,4,6,isg)
c           print *,'k1,k2,k3=',k1,k2,k3,' isg=',isg
            d1=dist2n(co,c1,n)
            if (d1 .lt. d0-1.0) then
              ncha=ncha+1
              call trnsfr(rotopt,rot0,9)
              d0=d1
c              write (6,1781) ncha,k,d0,d1,rot0
c1781          format(' nch,k=',2i3,' d0,d1=',2e12.5,' rot0=',(/,3f8.4))
            end if
          end do
        end do
      end do
      if (ncha .gt. 0) then
        call rotate_c(c,n,rotopt,c,'SYMFIT',6)
        call matprod(rotopt,bestrot,bestrot)
      end if
      return
      end
