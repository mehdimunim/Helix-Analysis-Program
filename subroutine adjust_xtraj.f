      subroutine adjust_xtraj(xtraj,ifirst,ilast,increment,fact,iadjust)
      dimension xtraj(ilast)
c     print *,'ADJUST_XTRAJ ifirst,ilast,increment=',
c    -  ifirst,ilast,increment
      if (iadjust .eq. 0) return
      iadjust=0
      i=0
      if (ifirst .gt. 1 .or. increment .gt. 1) then
        do ii=ifirst,ilast,increment
          i=i+1
          xtraj(i)=i*increment*fact
          MAXI=i
        end do
      end if
      return
      end
