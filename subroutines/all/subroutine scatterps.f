      subroutine scatterps(ips,xmin,xmax,xm,ymin,ymax,ym,x,y,ix1,iy1,n,
     -  title,lentit,noclose)
      dimension x(2,n),y(2,n)
      character*(*) title
      write (ips,1005)
      x0=0.1*xm
      y0=0.1*ym
      write (ips,1001) x0,y0
      write (ips,1002) x0,ym
      write (ips,1002) xm,ym
      write (ips,1002) xm,y0
      write (ips,1002) x0,y0
      write (ips,1000)
      write (ips,1005)
      write (ips,1001) x0,0.01*ym
      write (ips,1004) xmin
      write (ips,1001) 0.9*xm,0.01*ym
      write (ips,1004) xmax
      write (ips,1001) 0.0,0.0
      write (ips,1004) ymin
      write (ips,1001) 0.0,ym
      write (ips,1004) ymax
      write (ips,1001) x0,1.01*ym
      call psshow(ips,title,lentit)
      write (ips,1000)
      do i=1,n
        xx=x0+(xm-x0)*(x(ix1,i)-xmin)/(xmax-xmin)
        yy=y0+(ym-y0)*(y(iy1,i)-ymin)/(ymax-ymin)
        write (ips,1005)
        write (ips,1003) xx,yy,3,0,360
        write (ips,1000)
      end do
      write (ips,*) 'showpage'
      if (noclose .eq. 0) close (ips)
      return
1000  format('stroke')
1001  format(2f10.1,' moveto')
1002  format(2f10.1,' lineto')
1003  format(2f10.1,3i5,' arc')
1004  format('(',f10.2,') show')
1005  format('newpath')
      end
