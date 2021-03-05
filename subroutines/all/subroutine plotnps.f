      subroutine plotnps(x,y,nxmax,nymax,nf,imf,iml,ifg,r0,cx,
     -  x00,xd,nx,y00,yd,ny,iprt,ntit,tit,ntit2,tit2,xlab,lenx,
     -  fclab,lfclab,imarkx,imarks,marks,iplot,ipspage,npspages,
     -  inperr,iout)
c#    MMC routine 354 lstmod: 04/05/08
c*****Postscript plot of n functions of the same variable
      dimension x(nxmax),y(nymax),imf(nf),iml(nf),ifg(nf),lfclab(nf),
     -  imarks(nxmax)
      character*(*) xlab,fclab(nf),tit,tit2
      character*1 marks(9)
c     character*80 ident
c     common /title/ ident(2)
c     nf: number of functions to be plotted
c     x: the x coordinates of the functions to be plotted
c     ifg(if): First value in x for the if-th function
c     imf(if), iml(if): Y {imf(if) - iml(if)} the value of the if-th function
c     r0,cx: x coordonate labels are transformed as r0+x(i)*cx
c     y00,yd: y scale minimum and unit, yd=0 => program finds them
c     iprt: if .ne. 0, print the function values;
c     tit: string containing the title; ntit: number of chars in tit
c     print *,'ny,iprt,ntit,ntit2=',ny,iprt,ntit,ntit2
c     print *,'nxmax,nymax,nf=',nxmax,nymax,nf
c     print *,'tit=',tit(1:ntit)
c     print *,'tit2=',tit2(1:ntit2)
c     print *,'fclabs=',(fclab(i)(1:lfclab(i)),i=1,nf)
c     print *,'PLOTNPS ipspage,npspages,iout=',ipspage,npspages,iout
      if (iplot .eq. 0) return
      if (iprt .gt. 0) then
        write (iout,5002)
        do if=1,nf
          write (iout,5005) if,x(ifg(if)),x(ifg(if)+iml(if)-imf(if)),
     -      (y(ig),ig=imf(if),iml(if))
          if (iprt .gt. 1) write (iout,5006) (x(ig),ig=imf(if),iml(if))
        end do
      end if
      ix0=100
      iy0=90
      ixwid=600
      iyhgt=480
      call psheader(iplot,tit,ntit,0,0,612,792,npspages,ipspage)
      write (iplot,1001)
      if (nf .gt. 1) call rgbrainbowcolor(iplot,-1.0)
      write (iplot,1002) ix0+300-ntit*3,iy0+iyhgt+10
      write (iplot,1004) tit(1:ntit)
      write (iplot,1002) ix0+300,iy0-40
      iytit=iy0+iyhgt-15
      if (ntit2 .gt. 0) then
        write (iplot,1002) ix0+10,iytit
        write (iplot,1004) tit2(1:ntit2)
        iytit=iytit-15
      end if
c     write (iplot,1002) ix0+10,iytit
c     write (iplot,1004) ident(1)
c     iytit=iytit-15
c     write (iplot,1002) ix0+10,iytit
c     write (iplot,1004) ident(2)
c     write (iplot,1010) ix0,iy0,ixwid,iyhgt,-ixwid
      if (lenx .gt. 0) then
        write (iplot,1002) ix0+300,iy0-30-lenx/5
        write (iplot,1004) xlab(1:lenx)
      end if
      iyhgt=iyhgt-40
      if (xd .eq. 0.0) then
        write (iout,5001)
        inperr=inperr+1
        return
      else
        xmin=x00
        xdiv=xd
        xmax=xmin+nx*xdiv
      end if
      if (yd .eq. 0.0) then
        ymin=y(1)
        ymax=ymin
        do i=1,iml(nf)
          if (y(i) .lt. ymin) ymin=y(i)
          if (y(i) .gt. ymax) ymax=y(i)
        end do
        ny=10
        ydiv=(ymax-y00)/10.0
      else
        ymin=y00
        ydiv=yd
        ymax=ymin+ny*ydiv
      end if
c     Drow graph boundary box
      write (iplot,*) 'np'
      write (iplot,1002) ix0,iy0
      write (iplot,1007) ix0,iy0+iyhgt
      write (iplot,1007) ix0+ixwid,iy0+iyhgt
      write (iplot,1007) ix0+ixwid,iy0
      write (iplot,1007) ix0,iy0
      write (iplot,*) 'sk'
c     Draw ticks, write axis values
      ly=alog10(abs(ydiv))
      do i=1,ny
        write (iplot,1008)
        write (iplot,1002) ix0,iy0+i*iyhgt/ny
        write (iplot,1003) 5,0
        if (ly .lt. 7) then
          write (iplot,1002) ix0-50,iy0-5+i*iyhgt/ny
          write (iplot,1005) ymin+i*ydiv
        else
          write (iplot,1002) ix0-60,iy0-5+i*iyhgt/ny
          write (iplot,1006) ymin+i*ydiv
        end if
        write (iplot,1002) ix0+ixwid,iy0+i*iyhgt/ny
        write (iplot,1003) -5,0
        write (iplot,1011)
        write (iplot,1002) ix0+ixwid+10,iy0-5+i*iyhgt/ny
        if (ly .lt. 7) write (iplot,1005) ymin+i*ydiv
        if (ly .ge. 7) write (iplot,1006) ymin+i*ydiv
      end do
      if (imarkx .gt. 0) then
c       Mark residues
        maxx=0
        do if=1,nf
          nxif=ifg(if)+iml(if)-imf(if)
          if (nxif .gt. maxx) maxx=nxif
        end do
c       write (6,9292) (imarks(i),i=1,120 )
c9292    format(50i1)
        do i=1,maxx
          if (imarks(i) .gt. 0) then
            ix=(x(i)-xmin)*ixwid/(xmax-xmin)
            write (iplot,1002) ix0+ix,iy0+iyhgt
            write (iplot,1004) marks(imarks(i))
          end if
        end do
      end if
      write (iplot,1011)
      do i=1,nx
        write (iplot,1008)
        write (iplot,1002) ix0+i*ixwid/nx,iy0
        write (iplot,1003) 0,5
        write (iplot,1011)
        write (iplot,1008)
        write (iplot,1002) ix0+i*ixwid/nx,iy0+iyhgt
        write (iplot,1003) 0,-5
        write (iplot,1011)
        write (iplot,1002) ix0-30+i*ixwid/nx,iy0-15
        write (iplot,1013) r0+i*cx*xdiv
      end do
      write (iplot,1011)
c     Plot graphs
      if (ymax .gt. ymin) then
        yfac=iyhgt/(ymax-ymin)
      else
        yfac=0.0
      end if
c     iyhgt=480
      iylab=iy0+iyhgt-15
      ixinc=max0(1,nx/15)
      maxfclab=0
      do if=1,nf
        if (maxfclab .lt. lfclab(if)) maxfclab=lfclab(if)
      end do
      do if=1,nf
        istarted=0
        if (nf .gt. 1)
     -    call rgbrainbowcolor(iplot,float(if-1)/float(nf-1))
        write (iplot,1016) 1+2*if
        write (iplot,1008)
        write (iplot,1002) ix0+3*ixwid/4,iylab
        write (iplot,1012) if
        write (iplot,1002) ix0+3*ixwid/4+15,iylab+4
        write (iplot,1003) max0(ixwid/12,ixwid/4-6*maxfclab-30),0
        write (iplot,1002) ix0+3*ixwid/4+15+max0(ixwid/12,
     -    ixwid/4-6*maxfclab-30),iylab
        write (iplot,1015) fclab(if)(1:lfclab(if))
        write (iplot,1011)
        iylab=iylab-15
        write (iplot,1008)
c       print *,'if,imf(if),iml(if)=',if,imf(if),iml(if)
        do ig=imf(if),iml(if)
c         write (77,*) 'if,ig,y(ig)=',if,ig,y(ig)
          if (y(ig) .ne. 0.0) then
            ix=ix0+(x(ifg(if)+ig-imf(if))-xmin)*ixwid/(xmax-xmin)
            iy=iy0+yfac*(y(ig)-ymin)
            if (istarted .eq. 0) then
              write (iplot,1002) ix,iy
              istarted=1
            else
              write (iplot,1007) ix,iy
            end if
          end if
        end do
        write (iplot,1011)
      end do
      if (nf .gt. 1) call rgbrainbowcolor(iplot,-1.0)
      write (iplot,1014)
      write (iplot,1099)
      return
1001  format('/Helvetica findfont',/,'11 scalefont',/,'setfont',/,
     -  '612 0 translate',/,'90 rotate')
1002  format(i4,i4,' moveto')
1003  format(i5,i5,' rlineto')
1004  format('(',a,') show')
1005  format('(',f10.3,') show')
1006  format('(',e10.4,') show')
1007  format(i5,i5,' lineto')
1008  format('newpath')
c1010  format('% Drawing of graph boundaries',/,'newpath',/,
c     -  i3,1x,i3,' moveto',/,i4,' 000 rlineto',/,'000 ',i4,' rlineto',/,
c     -  i5,' 000 rlineto',/,'closepath',/,'stroke')
1011  format('stroke')
1012  format('(',i2,') show')
1013  format('(',f10.2,') show')
1014  format('-90 rotate',/,'-612 0 translate')
1015  format('( : ',a,') show')
1016  format('[',i3,' 3] 0 setdash')
1099  format('showpage')
5001  format(' ***** PROGRAM ERROR: plotnps needs nonzero xdiv')
5002  format(1x,//,(1x,a))
5005  format(' Plot ',i5,' xfirst=',e12.5,' xlast=',e12.5,' y=',/,
     -  (10e13.6))
5006  format(//)
      end
