      subroutine angles(e1,e2,e3,ca1,ca2,ca3)
      ee1=sqrt(e1)
      ee2=sqrt(e2)
      ee3=sqrt(e3)
      ca1=(e2+e3-e1)/(2.0*ee2*ee3)
      ca2=(e1+e3-e2)/(2.0*ee1*ee3)
      ca3=(e1+e2-e3)/(2.0*ee1*ee2)
      return
      end
