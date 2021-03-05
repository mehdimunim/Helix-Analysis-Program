      subroutine angledials(c,n,nangsel,ixang123,index,line,pi,
     -  ir1,ir2,ic1,ic2,iw0,mxrec)
      dimension c(3,n),index(n),ixang123(4,nangsel)
      character* 132 line(mxrec)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      do it=1,nangsel
        phirad=angleijk(c,n,ixang123(1,it),ixang123(2,it),
     -    ixang123(3,it),iw0)
        phi=phirad*(180.0/pi)
        res(1,nframe,it)=cos(phirad)
        res(2,nframe,it)=sin(phirad)
        write (iw0,1000) it,phi,(line(index(ixang123(k,it)))(ic1:ic2),
     -    ixang123(k,it),k=1,3),
     -    line(index(ixang123(2,it)))(ir1:ir2)
      end do
      return
1000  format(' A',i2,'=',f8.3,' (',a4,i6,')',2(' - (',a4,i6,')'),
     -  ' R2:',a)
      end
