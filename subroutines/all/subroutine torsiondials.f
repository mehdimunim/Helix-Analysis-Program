      subroutine torsiondials(c,n,ntorsel,ixtor1234,index,line,pi,
     -  ir1,ir2,ic1,ic2,iw0,mxrec)
      dimension c(3,n),index(n),ixtor1234(4,ntorsel)
      character* 132 line(mxrec)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
c     print *,'TORSIONDIALS n,mxrec=',n,mxrec
      call trajlimtest(nframe,MAXFRAMES)
      do it=1,ntorsel
        phirad=dihangl(c,ixtor1234(1,it),ixtor1234(2,it),
     -    ixtor1234(3,it),ixtor1234(4,it),0,mxrec)
        phi=phirad*(180.0/pi)
        res(1,nframe,it)=cos(phirad)
        res(2,nframe,it)=sin(phirad)
        write (iw0,1000) it,phi,(line(index(ixtor1234(k,it)))(ic1:ic2),
     -    ixtor1234(k,it),k=1,4),
     -    line(index(ixtor1234(3,it)))(ir1:ir2)
      end do
      return
1000  format(' T',i3,'=',f8.3,' (',a4,i6,')',3(' - (',a4,i6,')'),
     -  ' R3:',a)
      end
