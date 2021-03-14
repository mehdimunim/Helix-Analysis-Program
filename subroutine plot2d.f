      subroutine plot2d(iplot,xy,n,nfravgt,xmn,xdv,nxdv,ymn,ydv,nydv,
     -  title,ntit,plotdesc,lplotdesc,remark,lremark,xlab,lxlab,
     -  ylab,lylab,timelab,ltimelab,iperiod,iprt,iout,npspages,ipspage,
     -  noclose,maxplot)
      dimension xy(2,maxplot)
      character*(*) title,plotdesc,remark,xlab,ylab,timelab
c     print *,'PLOT2D xmn,xdv,nxdv=',xmn,xdv,nxdv
c     print *,'PLOT2D ymn,ydv,nydv=',ymn,ydv,nydv
c     print *,'START PLOT2D iout,npspages,ipspage=',
c    - iout,npspages,ipspage
      write (iplot,*) '% START plot2d: ',ylab(1:lylab)
      if (nfravgt .eq. 0)
     -  call getint(
     -    'Number of snapshots to average in the 2D trace plots',43,
     -    1,1,nframe,nfravgt,142)
      incr=max0(1,nfravgt)
      ix0=110
      iy0=150
      ixwid=400
      iyhgt=400
      call psheader(iplot,title,ntit,0,0,612,792,npspages,ipspage)
      ipspage=ipspage+1
      write (iplot,1020) ipspage
      write (iplot,1001) ix0+ixwid/2-50,iy0-40
      call psshow(iplot,xlab,lxlab)
      write (iplot,1001) 45,iy0+iyhgt/2-50
      write (iplot,1018) 90
      call psshow(iplot,ylab,lylab)
      write (iplot,1018) -90
      ixcent=80
      if (ntit .gt. 50) ixcent=0
      write (iplot,1002) ix0+ixcent,iy0+iyhgt+10
      write (iplot,1017)
      call psshow(iplot,title,ntit)
c     write (iplot,1002) ix0+ixwid,iy0+iyhgt+10
      write (iplot,1010) ix0,iy0,ixwid,iyhgt,-ixwid
      if (lplotdesc .gt. 0) then
c       Print remark (if any)
        write (iplot,1002) ix0,iy0+iyhgt+25
        call psshow(iplot,plotdesc,lplotdesc)
      end if
      if (lremark .gt. 0) then
c       Print remark (if any)
        write (iplot,1002) ix0,iy0+iyhgt-20
        call psshow(iplot,remark,lremark)
      end if
      write (iplot,1011)
c     Find X,Y range
      call arminmax2(xy,1,n,2,xmin,xmax,ymin,ymax,0,2)
      if (nxdv .gt. 0) then
        xmin=xmn
        xdiv=xdv
        nx=nxdv
        xmax=xmin+nx*xdiv
      else
        nx=10
        xdiv=(xmax-xmin)/10.0
      end if
      if (nydv .gt. 0) then
        ymin=ymn
        ydiv=ydv
        ny=nydv
        ymax=ymin+ydiv*ny
      else
        ny=10
        ydiv=(ymax-ymin)/ny
      end if
      lx=alog10(amax1(abs(xmin),abs(xmax)))
c     Draw ticks, write axis values
      if (xmin .ne. 0.0) then
c       Plot initial value too
        ixmin=0
      else
        ixmin=1
      end if
c     Draw ticks, write axis values
      write (iplot,1007) 'Draw ticks, axis values'
      if (ymin .ne. 0.0) then
c       Plot initial value too
        iymin=0
      else
        iymin=1
      end if
      ly=alog10(amax1(abs(ymin),abs(ymax)))
      do i=iymin,ny
c       Left tick
        write (iplot,1012)
        write (iplot,1002) ix0,iy0+i*iyhgt/ny
        write (iplot,1003) 5,0
c       write (iplot,1011)
c       Right tick
c       write (iplot,1012)
        write (iplot,1002) ix0+ixwid,iy0+i*iyhgt/ny
        write (iplot,1003) -5,0
        write (iplot,1011)
        write (iplot,1002) ix0-60,iy0-5+i*iyhgt/ny
        if (ny .lt. 8 .or. mod(i,2) .eq. 0) then
          if (ly .lt. 7) write (iplot,1005) ymin+i*ydiv
          if (ly .ge. 7) write (iplot,1006) ymin+i*ydiv
        end if
      end do
      write (iplot,1011)
      if (xmin .ne. 0.0) then
c       Plot initial value too
        ixmin=0
      else
        ixmin=1
      end if
      do i=ixmin,nx
c       Lower tick
        write (iplot,1012)
        write (iplot,1002) ix0+i*ixwid/nx,iy0
        write (iplot,1003) 0,5
c       write (iplot,1011)
c       Upper tick
c       write (iplot,1012)
        write (iplot,1002) ix0+i*ixwid/nx,iyhgt+iy0
        write (iplot,1003) 0,-5
        write (iplot,1011)
        write (iplot,1002) ix0-30+i*ixwid/nx,iy0-15
        if (nx .lt. 8 .or. mod(i,2) .eq. 0) then
          if (lx .lt. 7) write (iplot,1005) xmin+i*xdiv
          if (lx .ge. 7) write (iplot,1006) xmin+i*xdiv
        end if
      end do
c     Plot traces
      write (iplot,1007) 'Plot traces'
      if (xmax .gt. xmin) then
        xfac=ixwid/(xmax-xmin)
      else
        xfac=0.0
      end if
      if (ymax .gt. ymin) then
        yfac=iyhgt/(ymax-ymin)
      else
        yfac=0.0
      end if
      ixrange=xdiv*nx*xfac*0.8
      iyrange=ydiv*ny*yfac*0.8
      call rrgbcolor(iplot,1,n,1)
      newpath=1
      iprev=0
      intcol=n/100+1
      nn=max0(1,n/incr)
      xav=0.0
      yav=0.0
      iprev=0
      ixprev=ix0+xfac*(xy(1,1)-xmin)
      iyprev=iy0+yfac*(xy(2,1)-ymin)
      do i=1,n
        xav=xav+xy(1,i)
        yav=yav+xy(2,i)
        if (mod(i,incr) .eq. 0) then
          xav=xav/float(incr)
          yav=yav/float(incr)
          ix=ix0+xfac*(xav-xmin)
          iy=iy0+yfac*(yav-ymin)
          xav=0.0
          yav=0.0
          if (i .gt. 1) then
            idraw=1
            if (iperiod .eq. 1) then
c             Don't draw lines for moves under periodicity
              if (iabs(ix-ixprev) .gt. ixrange .or.
     -            iabs(iy-iyprev) .gt. iyrange ) idraw=0
            end if
            if (idraw .eq. 1) then
              if (newpath .eq. 1) write (iplot,1002) ixprev,iyprev
              if (float(i-iprev)/float(n) .ge. 0.1) then
c               Draw 2-color line
                ixx=(ix+ixprev)/2
                iyy=(iy+iyprev)/2
                write (iplot,1009) ixx,iyy
                call rrgbcolor(iplot,i,n,1)
                write (iplot,1002) ixx,iyy
                write (iplot,1009) ix,iy
              else
                write (iplot,1009) ix,iy
              end if
            end if
            if (nn .lt. 200 .or. mod(i,intcol) .eq. 0 .or.
     -          idraw .eq. 0) then
              call rrgbcolor(iplot,i,n,1)
              newpath=1
            else
              newpath=0
            end if
          end if
          iprev=i
          ixprev=ix
          iyprev=iy
        end if
      end do
      write (iplot,1011)
      call rainbowscale(iplot,ix0,ixwid,iy0-100,n,0.0,0.0,0.0,
     -  timelab,ltimelab)
      write (iplot,1098)
      if (noclose .eq. 0) close (iplot)
      if (iprt .gt. 0) write (iout,2001) (i,(xy(k,i),k=1,2),i=1,n)
      return
1001  format('/Helvetica findfont',/,'11 scalefont',/,'setfont',/,
     -  i4,i5,' moveto')
1002  format(i5,i6,' moveto')
1003  format(i5,i6,' rlineto')
1005  format('(',f10.2,') show')
1006  format('(',e10.4,') show')
1007  format('% ',a)
1009  format(i5,i6,' lineto')
1010  format('% Drawing of graph boundaries',/,'newpath',/,
     -  i3,1x,i3,' moveto',/,i4,' 000 rlineto',/,'000 ',i4,' rlineto',/,
     -  i5,' 000 rlineto',/,'closepath',/,'stroke')
1011  format('stroke')
1012  format('newpath')
1017  format('(System: ) show')
1018  format(i5,' rotate')
1020  format('%%Page: 1 ',i4)
1098  format('-90 rotate',/,'-612 0 translate',/,'showpage')
c1099  format('%%Trailer')
2001  format(i11,2e14.6)
      end
