      subroutine vprod(r,i,j,k)
c#    MMC routine 025 lstmod: 02/06/86
c*****Computes the vector product of the columns i and j into k
      dimension r(3,3)
      r(1,k)=r(2,i)*r(3,j)-r(2,j)*r(3,i)
      r(2,k)=r(3,i)*r(1,j)-r(3,j)*r(1,i)
      r(3,k)=r(1,i)*r(2,j)-r(1,j)*r(2,i)
      return
      end
