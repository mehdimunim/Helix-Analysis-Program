      subroutine angcomp(v0,d0,v,ang)
      real*8 v0(3),d0(3),v(3),v01(3),ddot
c     Calculate the angle between v1 and v2.
c     Obtain the sign by requiring that d0 and d1 correspond to the z axis
      ang=dacoscheck(ddot(v0,v)/dsqrt(ddot(v0,v0)*ddot(v,v)),ccc,1,
     -  6,'ANGCOMP')
      call dcross(v0,v,v01)
      if (ddot(v01,d0) .lt. 0.d0) ang=-ang
      return
      end
