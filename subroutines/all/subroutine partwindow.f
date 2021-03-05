      subroutine partwindow(xm,ym,ymfac,nd,edge,x0,y0,ndincr,ndprow,
     -  nrows,maxdials)
      dimension x0(maxdials),y0(maxdials)
c     The window/page is xm by ym, (0,0) is the bottom left corner
c     The subroutine finds the upper left corners (x0(i),y0(i)) of nd cubes
c     that best fill the window/page
      y00=ym*ymfac
      nleft=nd
c     print *,'PARTWIN nd,ndprow=',nd,ndprow
      ndone=0
      nrows=0
      do while (nleft .gt. 0)
        ndo=min0(nleft,ndprow)
c       x00=(ndprow-ndo)*(xm-ndo*edge)/2.0
        x00=(xm-ndo*edge)/2.0
        do ir=1,ndo
          ndone=ndone+1
          x0(ndincr+ndone)=x00+(ir-1)*edge
          y0(ndincr+ndone)=y00
        end do
        nrows=nrows+1
        nleft=nleft-ndo
        y00=y00-edge*1.2
      end do
c     write (6,1000) xm,ym,edge,(i,x0(i),y0(i),i=1,nd)
      return
c1000  format(' xm=',f10.4,' ym=',f10.4,' edge=',f10.5,/,
c     -  i4,' x0=',f10.4,' y0=',f10.4)
      end
