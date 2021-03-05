      subroutine grelcd(cs,csr,natms,aw)
c#    MMC routine 017 lstmod: before 07/15/85
c*****compute the c.o.m. centered coordinates from cs into csr
      dimension cs(3,natms),csr(3,natms),aw(natms),cmcrd(3)
      call cofms(cs,cmcrd,natms,aw)
      do i=1,natms
        do k=1,3
          csr(k,i)=cs(k,i)-cmcrd(k)
        end do
      end do
      return
      end
