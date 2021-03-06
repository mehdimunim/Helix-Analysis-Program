      subroutine plotdssp(ips,ifss,ilss,itypss,nss,icx,maxconf,framefac,
     -  nxd,nyd,itit,ntit,xlab,nxlab,ylab,nylab,nframe,ifrdssp,ilrdssp,
     -  ixresno)
c*****Postscript plot of the development of DSSP assignments
      dimension ifss(nss),ilss(nss),itypss(nss),ixresno(ilrdssp)
      character*(*) itit,xlab,ylab
      dimension shiftl(9)
      character*1 typc
      character*21 ssname
      common /dsspnames/ lssname(9),ssname(9),typc(9)
      data shiftl /0.0,0.19,0.41,0.64,0.0,0.3,0.45,0.6,0.75/
c     nxd, nyd: increment between tics in the x, y axes, resp
c     itit: array containing the title; ntit: number of characters in itit
c     xlab, ylab: labels of the x, y axes
c     nxlab, ylab: number of characters in the labels of the x,y axes
c     print *,'PLOT nss,icx,maxconf=',nss,icx,maxconf
c     print *,'PlOT nxd,nyd,nframe=',nxd,nyd,nframe
      xm=675
      ym=495
      xmm=0.9*xm
      ymm=0.9*ym
c     xm0=0.075*xm
c     ym0=0.09*ym
      xm0=0.076*xm
      ym0=0.085*ym
      ixd=maxconf/nxd
      if (mod(maxconf,nxd) .ne. 0) ixd=ixd+1
      maxrsd=ilrdssp-ifrdssp+1
      iyd=maxrsd/nyd
      if (mod(maxrsd,nyd) .ne. 0) iyd=iyd+1
      xmax=nxd*ixd
      ymax=nyd*iyd
c     print *,'ixd,iyd,xmin,xmax=',ixd,iyd,xmin,xmax
      if (nframe .eq. 1) then
        call psheader(ips,itit,ntit,-30,-130,830,830,1,ipspage)
        write (ips,1006)
        write (ips,3000)
        write (ips,*) '-90 rotate'
        write (ips,1001) -xm*1.1,0.1*ym,' translate'
        write (ips,*) 'np'
        write (ips,1001) xm0,ym0,' m'
        write (ips,1001) xm0,ym0+ymm,' l'
        write (ips,1001) xm0+xmm,ym0+ymm,' l'
        write (ips,1001) xm0+xmm,ym0,'   l'
        write (ips,1001) xm0,ym0,' l'
        write (ips,*) 'sk'
        write (ips,*) 'np'
c       write (ips,1001) xm0+0.3*xmm,ym0+1.03*ymm,' m'
        write (ips,1001) xm0,ym0+1.03*ymm,' m'
        call psshow(ips,itit,ntit)
        write (ips,1001) xm0+0.50*xmm,ym0+1.02*ymm,' m'
        write (ips,*) 'sk'
        write (ips,*) 'np'
c       write (ips,1001) xm0+0.45*xmm,ym0-0.10*ymm,' m'
        write (ips,1001) xm0+0.45*xmm,ym0-0.07*ymm,' m'
        call psshow(ips,xlab,nxlab)
        write (ips,*) 'sk'
        write (ips,1001) xm0-0.06*xmm,ym0+0.4*ymm,' translate'
        write (ips,*) '90 rotate'
        write (ips,*) 'np'
        write (ips,1001) 0.0,0.0,' m'
        call psshow(ips,ylab,nylab)
        write (ips,*) 'sk'
        write (ips,*) '-90 rotate'
        write (ips,1001) -(xm0-0.06*xmm),-(ym0+0.4*ymm),' translate'
        write (ips,*) 3,' setlinewidth'
        yinc=0.10*ymm
        do is=1,9
          if (is .eq. 5) yinc=yinc+0.05*ymm
          call rgbcolor(ips,is)
          write (ips,*) 'np'
c         write (ips,1001) xm0+shiftl(is)*xmm,ym0-0.15*ymm,' m'
          write (ips,1001) xm0+shiftl(is)*xmm,ym0-yinc,' m'
          call psshow(ips,ssname(is),lssname(is))
          call psshow(ips,':',1)
          write (ips,1005) 0.01*xmm,0.006*ymm,' rmoveto'
          write (ips,1005) 0.040*xmm,0.0,' rlineto'
          write (ips,1005) 0.01*xmm,-0.006*ymm,' rmoveto'
          write (ips,*) 'sk'
        end do
        write (ips,*) 'np'
        write (ips,*) 1,' setlinewidth'
        call rgbcolor(ips,9)
        ifloatx=0
        if (xlab(1:1) .ne. ' ' .and. framefac*float(ixd) .lt. 100.0)
     -    ifloatx=1
        do ix=1,nxd
          write (ips,*) 'np'
          write (ips,1001) xm0+xmm*float(ix)/float(nxd),ym0,' m'
          write (ips,1001) 0.0,+0.01*ymm,' rlineto'
          write (ips,*) 'sk'
          write (ips,*) 'np'
          write (ips,1001) xm0+0.02*xmm+xmm*float(2*ix-1)/float(2*nxd),
     -      ym0-0.04*ymm,' m'
          if (ifloatx .eq. 0) then
            ixl=framefac*float(ix*ixd)
            write (ips,1002) ixl
          else
            write (ips,1007) framefac*float(ix*ixd)
          end if
        end do
        do iy=1,nyd
          write (ips,*) 'np'
          write (ips,1001) xm0,ym0+float(iy)/float(nyd)*ymm,' m'
          write (ips,1001) +0.01*xmm,0.0,' rlineto'
          write (ips,*) 'sk'
          write (ips,*) 'np'
          write (ips,1001) xm0-0.07*xmm,ym0-0.01*ymm+
     -      float(iy)/float(nyd)*ymm,' m'
          write (ips,1002) ixresno(ifrdssp-1+iy*iyd)
        end do
      end if
      rx=xm0+icx*xmm/xmax
      write (ips,*) 1,' setlinewidth'
      write (ips,*) 'np'
      itprev=0
      ym0=ym0-(ifrdssp-1)*ymm/ymax
      do iss=1,nss
        if (ilss(iss) .gt. ifrdssp .and. ifss(iss) .lt. ilrdssp) then
          if (itypss(iss) .ne. itprev) then
            if (itprev .ne. 0) write (ips,1003)
            call rgbcolor(ips,itypss(iss))
            itprev=itypss(iss)
          end if
          write (ips,1001) rx,ym0+max0(ifrdssp,ifss(iss))*ymm/ymax,'m'
          write (ips,1001) rx,ym0+min0(ilrdssp,ilss(iss))*ymm/ymax,'l'
        end if
      end do
      write (ips,*) 'sk'
      return
1001  format(2f8.1,1x,a)
1002  format('(',i8,') show',/,'sk')
1003  format('sk',/,'np')
1005  format(2f8.1,1x,a)
1006  format('%%Page: 1 1')
1007  format('(',f8.2,') show',/,'sk')
3000  format('/m { moveto } def',/,'/l { lineto } def',/,
     -  '/np { newpath } def',/, '/sk { stroke } def',/,
     -  '/f { fill } def',/,'/lw { setlinewidth } def',/,
     -  '/Helvetica findfont',/,'12 scalefont',/,'setfont')
      end
