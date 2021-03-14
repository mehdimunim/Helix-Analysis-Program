      subroutine copyintgreal(intgv,realv,itor)
c*****Copy the bitpattern between an integer and a real
      if (itor .eq. 0) then
c       Put the real into the integer
        call copyirr(realv,0)
        call copyiri(intgv,1)
      else
c       Put the integer into the real
        call copyiri(intgv,0)
        call copyirr(realv,1)
      end if
      return
      end
