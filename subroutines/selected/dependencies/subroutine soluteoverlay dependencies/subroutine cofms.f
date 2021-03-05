      subroutine cofms(csa,rmass,natm,aw)
c#    MMC routine 021 lstmod: 06/18/93
c*****compute center of mass from global atomic coordinates
      dimension csa(3,natm),rmass(3),aw(natm)
c     print *,'COFMS natm=',natm
c     print *,'COFMS aw=',aw
c     write (6,1003) csa
c1003  format(' csa=',/,(3f10.5))
      rmx=0.
      rmy=0.
      rmz=0.
      wmol=0.
      do i=1,natm
        wmol=wmol+aw(i)
        rmx=rmx+aw(i)*csa(1,i)
        rmy=rmy+aw(i)*csa(2,i)
        rmz=rmz+aw(i)*csa(3,i)
      end do
      if (wmol .eq. 0.0) wmol=1.0
      rmass(1)=rmx/wmol
      rmass(2)=rmy/wmol
      rmass(3)=rmz/wmol
      return
      end
