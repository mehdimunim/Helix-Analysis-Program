      subroutine rescale_bl(c1,c2,cn,ro,rn)
      dimension c1(3),c2(3),cn(3)
      dimension e(3)
      do k=1,3
        e(k)=(c2(k)-c1(k))/ro
      end do
      do k=1,3
        cn(k)=c1(k)+e(k)*rn
      end do
      return
      end
