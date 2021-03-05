      subroutine combinevol(cmin1,cmax1,cmin2,cmax2,cmin12,cmax12,vol)
      dimension cmin1(3),cmax1(3),cmin2(3),cmax2(3),cmin12(3),cmax12(3)
      vol=1.0
      do k=1,3
        cmin12(k)=amin1(cmin1(k),cmin2(k))
        cmax12(k)=amax1(cmax1(k),cmax2(k))
        vol=vol*(cmax12(k)-cmin12(k))
      end do
      return
      end
