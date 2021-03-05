      subroutine copyiri(intgv,iread)
      common /janus/ intgvar
      if (iread .eq. 1) intgv=intgvar
      if (iread .eq. 0) intgvar=intgv
      return
      end
