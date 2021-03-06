      subroutine dvprd(a,b,c)
c#    MMC routine 048  lstmod: 11/13/90
c*****Computes the vector product a x b and saves it into c
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
      return
      end
