      subroutine blockfromcum(bl,cum,xcum,n)
c#    MMC routine 257 lstmod: 03/20/01
c*****Extract block averages from cumulative sum
      real*8 cum,xcum,denom
      dimension bl(n),cum(n),xcum(n)
      if (n .lt. 1) return
      if (xcum(1) .eq. 0.d0) bl(1)=0.d0
      if (xcum(1) .ne. 0.d0) bl(1)=cum(1)/xcum(1)
      do i=2,n
        denom=xcum(i)-xcum(i-1)
        if (denom .eq. 0.d0) bl(i)=bl(i-1)
        if (denom .ne. 0.d0) bl(i)=(cum(i)-cum(i-1))/denom
      end do
      return
      end
