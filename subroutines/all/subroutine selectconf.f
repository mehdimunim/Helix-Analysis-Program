      subroutine selectconf(numsel,ninconf,ifirst,increment,
     -  iconfsel,nextconfsel,nframeref,icsel,ifail,maxconfsel)
      dimension iconfsel(maxconfsel)
      if (numsel .eq. 0) then
        if (ninconf .ge. ifirst .and.
     -    mod(ninconf-ifirst,increment) .eq. 0) icsel=1
      else
c       See if this configuration is on the list
        if (ninconf .eq. iconfsel(nextconfsel)) then
          icsel=1
          nextconfsel=nextconfsel+1
        else if (ninconf .gt. iconfsel(numsel)) then
          if (nframeref .le. 1)
     -      write (6,2060) ninconf,iconfsel(numsel)
          ifail=1
        end if
      end if
      return
2060  format(' Structure read (',i9,') is beyond the last structure ',
     -  'requested:',i9,/,' Trajetory scan stopped')
      end
