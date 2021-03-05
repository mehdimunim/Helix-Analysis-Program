      subroutine makewat(c,na,roh,ahohd,radtodeg)
      dimension c(3,na)
c     Create a standard water
      do i=1,na
        c(i,1)=0.0
        c(3,i)=0.0
      end do
      ahoh=ahohd/(radtodeg*2.0)
      xw=roh*cos(ahoh)
      yw=roh*sin(ahoh)
      c(1,2)=xw
      c(1,3)=xw
      c(2,2)=yw
      c(2,3)=-yw
      return
      end
