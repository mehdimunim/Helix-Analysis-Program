      subroutine vprd(a,b,c)
c#    MMC routine  61 lstmod: 02/06/86
c*****Computes the vector product a x b and saves it into c
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
      return
      end
