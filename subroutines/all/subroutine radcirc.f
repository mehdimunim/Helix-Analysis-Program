      subroutine radcirc(a,b,c,r)
      dimension a(3),b(3),c(3)
c     r is the radius of the circle going through a,b, and c
      call angdistw(b,a,c,rHB,rb,rac,angabc)
      sabc=sin(angabc)
      if (sabc .lt. 0.00001) then
        r=999999.9
      else
        r=rac/(sin(angabc)*2.0)
      end if
      return
      end
