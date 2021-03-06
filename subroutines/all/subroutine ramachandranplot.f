      subroutine ramachandranplot(nres,ips,xm,iallrama)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
c     print *,'RAMACHANDRANPLOT nres,ips,maxpres=',nres,ips,maxpres
      xmm=0.9*xm
      ymm=0.9*xm
      xm0=0.09*xm
      ym0=0.20*xm
      if (nframe .le. 1) then
        write (ips,*) 'np'
        write (ips,1001) xm0,ym0,' m'
        write (ips,1001) xm0,ym0+ymm,' l'
        write (ips,1001) xm0+xmm,ym0+ymm,' l'
        write (ips,1001) xm0+xmm,ym0,'   l'
        write (ips,1001) xm0,ym0,' l'
        write (ips,*) 'sk'
        write (ips,1001) xm0,ym0-0.03*ymm,'m'
        write (ips,*) '(-180) show'
        write (ips,1001) xm0+xmm-0.05*xmm,ym0-0.03*ymm,'m'
        write (ips,*) '(+180) show'
        write (ips,1001) xm0+xmm/2.0-0.06*xmm,ym0-0.03*ymm,'m'
        write (ips,*) '(Phi) show'
        write (ips,1001) xm0-0.07*xmm,ym0+0.01*ymm,'m'
        write (ips,*) '(-180) show'
        write (ips,1001) xm0-0.07*xmm,ym0+ymm-0.02*ymm,'m'
        write (ips,*) '(+180) show'
        write (ips,1001) xm0-0.07*xmm,ym0+ymm/2.0-0.02*ymm,'m'
        write (ips,*) '(Psi) show'
        write (ips,*) 'np'
        write (ips,1001) xm0,ym0+1.01*ymm,' m'
      end if
      if (nframe .gt. 0) call rrgbcolor(ips,nframe,nframetot,0)
      nresplot=nxselres
      if (iallrama .eq. 1) nresplot=nres
      i=1
      do while (i .le. nresplot)
        ii=i
        if (iallrama .eq. 0) ii=ixselres(i)
        if (res(1,i,maxpres) .lt. 999.0) then
          if (nframe .eq. 0) call rrgbcolor(ips,i,nres,0)
c         write (77,*) i,' phi,psi=',(res(k,i,maxpres),k=1,2)
          xx=xm0+xmm*(res(1,i,maxpres)+180.0)/360.0
          yy=ym0+ymm*(res(2,i,maxpres)+180.0)/360.0
          write (ips,*) 'np'
          write (ips,1001) xx,yy,' 1 0 360 arc'
          write (ips,*) 'sk'
        end if
        i=i+1
      end do
      if (nframe .eq. 0) then
        ixmm=xmm
        ixm0=xm0
        call rainbowscale(ips,ixm0,ixmm,25,nres,0.0,0.0,0.0,'N(res)',6)
        write (ips,*) 'showpage'
        close (ips)
      end if
      return
1001  format(2f8.1,1x,a)
      end
