      subroutine shiftmol(c,n,c0,cnew,fac)
      dimension c(3,n),cnew(3,n),c0(3)
C     Center the molecule at c0
c     print *,'SHIFTMOL n=',n,' fac=',fac,' c0=',c0
      do i=1,n
        do k=1,3
          cnew(k,i)=c(k,i)+fac*c0(k)
        end do
      end do
      return
      end
