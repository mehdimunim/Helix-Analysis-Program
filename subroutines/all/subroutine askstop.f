      subroutine askstop(idef)
      common /logging/ logfile,ipredict
      if (ipredict .eq. 0) then
        call askyn('Do you want to continue',23,1,idef,icont,20,0)
        if (icont .eq. 0) stop
      else
        print *,'Run continues as predictable input was requested'
      end if
      return
      end
