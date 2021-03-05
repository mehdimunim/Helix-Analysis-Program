      subroutine printhelix(iw0,axisini,axisend,cent,rms,helixlen,
     -  axisdir,angles,decidebend,nup,ndown,nrun,nnear,rcirc,turnperres,
     -  anglesn,ihx,radtodeg)
      real*8 axisdir(3),axisini(3),axisend(3),cent(3),rms
      dimension angles(3),anglesn(3)
      character*9 decidebend
      write (iw0,1000) ihx,axisini,axisend,helixlen,axisdir,
     -  (radtodeg*angles(k),k=1,3),rms,decidebend,nup,ndown,nrun-1,
     -  nnear,rcirc,radtodeg*turnperres,cent,(radtodeg*anglesn(k),k=1,3)
      return
1000  format(' HX#',i3,' S=',3f9.4,' E=',3f9.4,' Len=',f5.2,/,
     -  ' D=',3f10.6,' D-X,D-Y,D-Z angles=',3f7.2,/,
     -  ' RMS=',f5.2,' Shape:',a9,' Nup/dn=',2i3,' Ncross=',i2,
     -  ' Nax=',i2,' Rc=',f6.1,' TPR=',f6.2,/,
     -  ' C=',3f10.5,' N-X,N-Y,N-Z angles=',3f7.2)
      end
