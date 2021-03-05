      subroutine rainbowscale(iplot,ix0,ixwid,iy0,n,rn,rmin,rmax,
     -  label,llabel)
      character*(*) label
      write (iplot,1019) 15
      ixprev=ix0
      do i=1,100
        call rrgbcolor(iplot,i,100,1)
        ix=ix0+ixwid*float(i)/float(100)
        write (iplot,1002) ixprev,iy0
        write (iplot,1009) ix,iy0
        ixprev=ix
      end do
      write (iplot,1011)
      call rrgbcolor(iplot,i,1,0)
      write (iplot,1002) ix0-30,iy0-25
      if (n .eq. 0 .and. rn .eq. 0.0) then
c       Use rcmin,rcmax
        write (iplot,1014) rmin
        ixlabwid=8.0*alog10(amax1(10.1,rmax))
        write (iplot,1002) ix0+ixwid-30-ixlabwid,iy0-25
        write (iplot,1014) rmax
      else
c       Use 0 - n/rn
        write (iplot,1013) 0
        if (rn .gt. float(n)) then
          ixlabwid=8.0*alog10(amax1(10.1,rn))
          write (iplot,1002) ix0+ixwid-30-ixlabwid,iy0-25
          write (iplot,1014) rn
        else
          ixlabwid=8.0*(alog10(float(n))+2.0)
          write (iplot,1002) ix0+ixwid-ixlabwid,iy0-25
          write (iplot,1013) n
        end if
      end if
      write (iplot,1002) ix0+ixwid/2-50,iy0-25
      call psshow(iplot,label,llabel)
      return
1002  format(i4,i5,' moveto')
1009  format(i4,i5,' lineto')
1011  format('stroke')
1013  format('(',i10,') show')
1014  format('(',f10.4,') show')
1019  format(i5,' setlinewidth')
      end
