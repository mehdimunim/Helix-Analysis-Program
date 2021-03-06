      function dacoscheck(dargin,argin,idbp,iout,lab)
c*****Make sure cos is in the [-1,1] range, return dacos(arg)
      real*8 dargin,arg
      character*(*) lab
      if (idbp .eq. 1) then
        arg=dargin
      else
        arg=argin
      end if
      if (arg .lt. -1.0d0) then
        if (arg .ge. -1.01d0) then
          arg=-1.0d0
          dacoscheck=dacos(-1.0d0)
          return
        else
          write (iout,1000) arg,lab
          dacoscheck=0.d0
        end if
      else if (arg .gt. 1.0d0) then
        if (arg .le. 1.01d0) then
          arg=1.0d0
          dacoscheck=dacos(+1.0d0)
          return
        else
          write (iout,1000) arg,lab
          dacoscheck=dacos(+1.0d0)
        end if
      else
        dacoscheck=dacos(arg)
        return
      end if
      return
1000  format(' PROGRAM ERROR: invalid sin or cos value:',f10.6,
     -  ' in subroutine ',a)
      end
