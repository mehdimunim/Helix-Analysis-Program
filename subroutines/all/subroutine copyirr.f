      subroutine copyirr(realv,iread)
      common /janus/ realvar
      if (iread .eq. 1) realv=realvar
      if (iread .eq. 0) realvar=realv
      return
      end
