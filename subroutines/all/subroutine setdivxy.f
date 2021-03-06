      subroutine setdivxy(armin,armax,nxydiv,xydiv,xymin)
      if (armin .ne. 0.0) then
c      print *,'armin(inp)=',armin
       arminimp=armin
c       First, round off the minimum
        call roundlim(armin,xymindiv,nxymindiv)
        armin=xymindiv*(nxymindiv-1)
        do while (armin .gt. arminimp)
          armin=armin-xymindiv
        end do
      end if
      range=armax-armin
      call roundlim(range,xydiv,nxydiv)
      nmin=armin/xydiv
      xymin=nmin*xydiv
      return
      end
