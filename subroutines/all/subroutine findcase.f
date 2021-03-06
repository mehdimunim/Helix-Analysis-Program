      subroutine findcase(charin,icase)
      character*1 charin
      character*1 abc,idig,digits,hexdigits
      common /charactersets/ ihex(25),abc(26,2),idig(10),digits(14),
     -  hexdigits(25)
c     Returns icase=2: lower; 1: upper; 0: neither
      icase=0
      do ic=1,26
        if (charin .eq. abc(ic,1)) icase=1
        if (charin .eq. abc(ic,2)) icase=2
      end do
      return
      end
