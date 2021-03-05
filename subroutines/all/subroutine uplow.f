      subroutine uplow(charin,charout,iuptolow,noabc)
      character*1 charin,charout
      character*1 abc,idig,digits,hexdigits
      common /charactersets/ ihex(25),abc(26,2),idig(10),digits(14),
     -  hexdigits(25)
c     If iuptolow=1: up -> low; iuptolow=2: low -> up
      do ic=1,26
        if (charin .eq. abc(ic,iuptolow)) then
          charout=abc(ic,3-iuptolow)
          noabc=0
          return
        end if
      end do
      noabc=1
      return
      end
