      subroutine plot2fun(iplot,nplot,x,y12,sd12,n,xmn,xdv,nxdv,y1mn,
     -  y1dv,ny1dv,y2mn,y2dv,ny2dv,title,ntit,remark,lremark,xlab,lxlab,
     -  y1lab,ly1lab,y2lab,ly2lab,itrajfile,iprt,iout,nplotdim,
     -  nosd,nsdfreq_in,iconn,iyinc,npspages,ipspage,noclose,landscape,
     -  modhgt)
      dimension x(n),y12(nplotdim,n),sd12(nplotdim,n)
      character*(*) title,remark,xlab,y1lab,y2lab
      data ny2 /0/,ly2 /0/,iy2min /0/,y2div /0.0/,
     -  y1fac /0.0/,y2fac /0.0/,icol1 /7/,icol2 /7/
c     iyinc: increment the first index of y12 with iyinc (only for nplot=1)
c     print *,'PLOT2FUN n=',n,' nplot=',nplot,' yl1=',y1lab(1:ly1lab)
c     print *,'y2mn,y2dv,ny2dv=',y2mn,y2dv,ny2dv
c     print *,'XLAB  =',xlab  (1:lxlab  )
c     print *,'PLOT2FUN NOSD=',nosd
      if (n .lt. 1) return
      write (iplot,8791) '% START plot2fun: ',y1lab(1:ly1lab)
8791  format(a,a)
      if (iyinc .gt. 0 .and. nplot .gt. 1) then
        print *,'PROGRAM ERROR: iyinc > 0 and nplot > 1 in plot2fun'
        print *,'Called for ',title(1:ntit)
        return
      end if
c     do i=1,n
c        write (40,9789) x(i),y12(1,i),sd12(1,i),y12(2,i),sd12(2,i)
c9789    format(f6.0,' y1=',f9.4,' sd1=',f9.4,' y2=',f9.4,' sd2=',f9.4)
c     end do
      nsdfreq=nsdfreq_in
      ix0=150
      iy0=100
      ixwid=450
      if (modhgt .eq. 1 .and. ny1dv .gt. 0) then
        iyhgt=min0(420,20*ny1dv)
      else
        iyhgt=420
      end if
      if (nplot .eq. 2) then
        icol1=6
        icol2=1
      else
        icol1=9
      endif
      call psheader(iplot,title,ntit,0,0,612,792,npspages,ipspage)
      ipspage=ipspage+1
      write (iplot,1019) ipspage
      if (landscape .eq. 1) write (iplot,1001)
      write (iplot,1002) ix0+ixwid/2-50,iy0-40
      call psshow(iplot,xlab,lxlab)
      ixcent=80
      if (ntit .gt. 50) ixcent=0
      write (iplot,1002) ix0+ixcent,iy0+iyhgt+10
      call psshow(iplot,title,ntit)
      if (itrajfile .gt. 0) then
        write (iplot,1002) ix0+ixcent,iy0+iyhgt+25
        call write_traj_lim(iplot,' ',1,itrajfile,incr_tr,1)
      end if
c     write (iplot,1002) ix0+ixwid,iy0+iyhgt+10
      write (iplot,1011)
      write (iplot,1010) ix0,iy0,ixwid,iyhgt,-ixwid
c     Find X range
      if (nxdv .gt. 0) then
        xmin=xmn
        nx=nxdv
        if (xdv .eq. 0.0) then
          xmax=x(n)
          xdiv=(xmax-xmin)/float(nx)
        else
          xdiv=xdv
          xmax=xmin+nx*xdiv
        end if
      else
        xmin=0.0
        xmax=x(n)
        call roundlim(xmax,xdiv,nx)
      end if
c     Find Y1 range
      call arminmax2(y12,1,n,nplot,y1min,y1max,y2min,y2max,iyinc,
     -  nplotdim)
c     print *,'y1min,y1max=',y1min,y1max,' ny1dv=',ny1dv
c     print *,'y2min,y2max=',y2min,y2max,' ny2dv=',ny2dv
c     print *,'IYINC=',iyinc
      if (iyinc .eq. 1) then
c       Use 2nd column of y12 only
        if (ny2dv .gt. 0) then
          y1min=y2mn
          y1div=y2dv
          ny1=ny2dv
          y1max=y1min+y1div*ny1
        else
          ny1=10
          y1div=(y1max-y1min)/ny1
        end if
      else
        if (ny1dv .gt. 0) then
          y1min=y1mn
          y1div=y1dv
          ny1=ny1dv
          y1max=y1min+y1div*ny1
        else
          ny1=10
          y1div=(y1max-y1min)/ny1
        end if
      end if
      ly1=alog10(amax1(abs(y1min),abs(y1max)))
c     Draw ticks, write axis values
      iy1min=0
      if (nplot .eq. 2) then
        if (ny2dv .gt. 0) then
          y2min=y2mn
          y2div=y2dv
          ny2=ny2dv
          y2max=y2min+y2div*ny2
        else
          ny2=10
          y2div=(y2max-y2min)/ny2
        end if
        ly2=alog10(amax1(abs(y2min),abs(y2max)))
c       Draw ticks, write axis values
        iy1min=0
      end if
c     Draw y1 tics
      do i=iy1min,ny1
        write (iplot,1012)
        write (iplot,1002) ix0,iy0+i*iyhgt/ny1
        write (iplot,1003) 5,0
        write (iplot,1011)
      end do
      call rgbcolor(iplot,-icol1)
      write (iplot,1012)
c     Print y1 axis tick values
      do i=iy1min,ny1
        write (iplot,1002) ix0-50,iy0-5+i*iyhgt/ny1
        if (ly1 .lt. 0) then
          write (iplot,1023) y1min+i*y1div
        else if (ly1 .lt. 2) then
          write (iplot,1022) y1min+i*y1div
        else if (ly1 .lt. 4) then
          write (iplot,1007) y1min+i*y1div
        else if (ly1 .lt. 7) then
          iy1=y1min+i*y1div
          write (iplot,1005) iy1
        else
          write (iplot,1006) y1min+i*y1div
        end if
      end do
      write (iplot,1011)
      call rgbcolor(iplot,9)
      if (nplot .eq. 1) then
        iy2min=iy1min
        ny2=ny1
      end if
      write (iplot,8922) iy2min,ny2,ix0,iy0,iyhgt
8922  format('% iy2min,ny2,ix0,iy0,iyhgt=',5i6)
      do i=iy2min,ny2
        write (iplot,1012)
        write (iplot,1002) ix0+ixwid,iy0+i*iyhgt/ny2
        write (iplot,1003) -5,0
        write (iplot,1011)
      end do
      if (nplot .eq. 2) then
        call rgbcolor(iplot,-icol2)
        write (iplot,1012)
        do i=iy2min,ny2
          write (iplot,1002) ix0+ixwid+3,iy0-5+i*iyhgt/ny2
          if (ly2 .lt. 0) then
            write (iplot,1023) y2min+i*y2div
          else if (ly2 .lt. 2) then
            write (iplot,1022) y2min+i*y2div
          else if (ly2 .lt. 4) then
            write (iplot,1007) y2min+i*y2div
          else if (ly2 .lt. 7) then
            iy2=y2min+i*y2div
            write (iplot,1005) iy2
          else
            write (iplot,1006) y2min+i*y2div
          end if
        end do
        write (iplot,1011)
        call rgbcolor(iplot,9)
      end if
c     if (xmin .ne. 0.0) then
c       Plot initial value too
c       ixmin=0
c     else
c       ixmin=1
c     end if
      xrange=nx*xdiv
      do i=0,nx
        write (iplot,1012)
        write (iplot,1002) ix0+i*ixwid/nx,iy0
        write (iplot,1003) 0,5
        write (iplot,1002) ix0+i*ixwid/nx,iy0+iyhgt
        write (iplot,1003) 0,-5
        write (iplot,1011)
        write (iplot,1002) ix0-30+i*ixwid/nx,iy0-15
        if (xrange .gt. 10.01) then
          ix=xmin+i*xdiv+0.01
          write (iplot,1013) ix
        else
          xi=xmin+i*xdiv
          write (iplot,1008) xi
        end if
      end do
      write (iplot,1011)
      call rgbcolor(iplot,-icol1)
      write (iplot,1012)
      if (iconn .eq. 1) then
        write (iplot,1002) ix0+10,iy0+iyhgt-15
        write (iplot,1003) 40,0
        write (iplot,1016) 0,-3
      else
        write (iplot,1020) ix0+10,iy0+iyhgt-12,' 2 0 360 arc'
        write (iplot,1016) 0,-3
      end if
      call psshow(iplot,':',1)
      if (iyinc .eq. 0) call psshow(iplot,y1lab,ly1lab)
      if (iyinc .eq. 1) call psshow(iplot,y2lab,ly2lab)
      write (iplot,1011)
      if (lremark .gt. 0) then
c       Print remark (if any)
        write (iplot,1012)
        call rgbcolor(iplot,9)
        write (iplot,1002) ix0+10,iy0+iyhgt-30
        call psshow(iplot,remark,lremark)
        write (iplot,1011)
        call rgbcolor(iplot,-icol1)
      end if
      write (iplot,1012)
c     Plot graphs
      if (y1max .gt. y1min) y1fac=iyhgt/(y1max-y1min)
      if (y1max .eq. y1min) y1fac=0.0
      iprev=0
      if (nsdfreq .eq. 0) nsdfreq=max0(1,n/30)
      mody2=1
      if (nsdfreq .eq. 1) mody2=0
      xmaxgraph=amax1(xmax,nx*xdiv)
      write (iplot,1004) 'first function'
      do i=1,n
        if (x(i) .ge. xmin .and. x(i) .le. xmax .and.
     -      y12(iyinc+1,i) .ge. y1min .and.
     -      y12(iyinc+1,i) .le. y1max) then
          ix=ix0+ixwid*(x(i)-xmin)/(xmaxgraph-xmin)
          iy=iy0+y1fac*(y12(iyinc+1,i)-y1min)
          if (iconn .eq. 1) then
            if (iprev .eq. 0) then
              write (iplot,1002) ix,iy
              iprev=1
            else
              write (iplot,1009) ix,iy
            end if
          else
c           Just draw a circle
            write (iplot,1012)
            write (iplot,1020) ix,iy,' 2 0 360 arc'
            write (iplot,1011)
          end if
          if (nosd .eq. 0) then
c           Plot error bars
            if (mod(i,nsdfreq) .eq. 0) then
              iysd=y1fac*sd12(iyinc+1,i)
              if (iysd .gt. 0) then
                write (iplot,1003) 0,iysd
                write (iplot,1003) 3,0
                write (iplot,1003) -6,0
                write (iplot,1003) 3,0
                write (iplot,1003) 0,-2*iysd
                write (iplot,1003) 3,0
                write (iplot,1003) -6,0
                write (iplot,1003) 3,0
                write (iplot,1003) 0,iysd
              end if
            end if
          end if
        else
          iprev=0
        end if
      end do
      write (iplot,1004) 'second function'
      write (iplot,1011)
      if (nplot .eq. 2) then
        call rgbcolor(iplot,-icol2)
        if (iconn .eq. 1) write (iplot,1015)
        write (iplot,1012)
        if (iconn .eq. 1) then
          write (iplot,1002) ix0+ixwid/2,iy0+iyhgt-15
          write (iplot,1003) 40,0
          write (iplot,1016) 0,-3
        else
          write (iplot,1020) ix0+ixwid/2,iy0+iyhgt-12,' 2 0 360 arc'
          write (iplot,1021)
          write (iplot,1011)
          write (iplot,1012)
          write (iplot,1002) ix0+ixwid/2,iy0+iyhgt-13
        end if
        call psshow(iplot,':',1)
        call psshow(iplot,y2lab,ly2lab)
        write (iplot,1011)
        write (iplot,1012)
        if (y2max .gt. y2min) y2fac=iyhgt/(y2max-y2min)
        if (y2max .eq. y2min) y2fac=0.0
        iprev=0
        do i=1,n
c         write (6,6734) i,x(i),xmin,xmax,y12(2,i),y2min,y2max
c6734     format(i4,' X=',f8.3,' XMIN/MAX=',2f8.3,' Y12=',f8.3,
c    -      ' Y2MIN/MAX=',2f8.3)
          if (x(i) .ge. xmin .and. x(i) .le. xmax .and.
     -        y12(2,i) .ge. y2min .and. y12(2,i) .le. y2max) then
            ix=ix0+ixwid*(x(i)-xmin)/(xmaxgraph-xmin)
            iy=iy0+y2fac*(y12(2,i)-y2min)
            if (iconn .eq. 1) then
              if (iprev .eq. 0) then
                write (iplot,1002) ix,iy
                iprev=1
              else
                write (iplot,1009) ix,iy
              end if
            else
c             Just draw a circle
              write (iplot,1012)
              write (iplot,1020) ix,iy,' 2 0 360 arc'
              write (iplot,1021)
              write (iplot,1011)
            end if
            if (nosd .eq. 0) then
c             Plot error bars
              if (mod(i,nsdfreq) .eq. mody2) then
                iysd=y2fac*sd12(2,i)
                if (iysd .gt. 0) then
                  write (iplot,1003) 0,iysd
                  write (iplot,1003) 3,0
                  write (iplot,1003) -6,0
                  write (iplot,1003) 3,0
                  write (iplot,1003) 0,-2*iysd
                  write (iplot,1003) 3,0
                  write (iplot,1003) -6,0
                  write (iplot,1003) 3,0
                  write (iplot,1003) 0,iysd
                end if
              end if
            end if
          else
            iprev=0
          end if
        end do
        write (iplot,1011)
        call rgbcolor(iplot,9)
      end if
      if (landscape .eq. 1) write (iplot,1097)
      write (iplot,1098)
      if (noclose .eq. 0) close (iplot)
      if (iprt .gt. 0) then
        do i=1,n
          write (iout,2001) i,x(i),(y12(iyinc+k,i),k=1,nplot)
        end do
      end if
      return
1001  format('/Helvetica findfont',/,'11 scalefont',/,'setfont',/,
     -  '612 0 translate',/,'90 rotate')
1002  format(i5,i6,' m')
1003  format(i5,i6,' r')
1004  format('% Plot ',a)
1005  format('(',i8,') show')
1006  format('(',e10.4,') show')
1007  format('(',f8.2,') show')
1008  format('(',f10.2,') show')
1009  format(i5,i6,' l')
1010  format('% Drawing of graph boundaries',/,'newpath',/,
     -  i3,1x,i3,' moveto',/,i4,' 000 rlineto',/,'000 ',i4,' rlineto',/,
     -  i5,' 000 rlineto',/,'closepath',/,'stroke')
1011  format('sk')
1012  format('np')
1013  format('(',i10,') show')
1015  format('[2] 0 setdash')
1016  format(i5,i6,' rmoveto')

1019  format('%%Page: 1 ',i4)
1020  format(2i5,1x,a)
1021  format('f')
1022  format('(',f8.4,') show')
1023  format('(',f8.6,') show')
1097  format('-90 rotate',/,'-612 0 translate')
1098  format('showpage')
c1099  format('%%Trailer')
2001  format(i11,3e14.6)
      end
