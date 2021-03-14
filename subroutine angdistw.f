      subroutine angdistw(c1,c2,c3,roh1,roh2,rhh,ahoh)
c*****Calculate the angle c2-c1-c3 and distance c1-c2; c1-c3
      dimension c1(3),c2(3),c3(3)
      real*8 cosa,rroh1,rroh2,droh1,droh2
      rroh1=(c1(1)-c2(1))**2+(c1(2)-c2(2))**2+(c1(3)-c2(3))**2
      rroh2=(c1(1)-c3(1))**2+(c1(2)-c3(2))**2+(c1(3)-c3(3))**2
      rrhh=(c2(1)-c3(1))**2+(c2(2)-c3(2))**2+(c2(3)-c3(3))**2
      droh1=dsqrt(rroh1)
      droh2=dsqrt(rroh2)
      rhh=sqrt(rrhh)
      cosa=(rroh1+rroh2-rrhh)/(2.d0*droh1*droh2)
      ahoh=dacoscheck(cosa,ccc,1,6,'ANGDISTW')
      roh1=droh1
      roh2=droh2
      return
      end
