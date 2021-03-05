      subroutine savekinkdat(nmem,bend,wobble,faceshift,psr,radtodeg)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      parameter (MAXCOPY6=MAXCOPY-6)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,bdxy(2,MAXFRAMES),wbxy(2,MAXFRAMES),
     -  fsxy(2,MAXFRAMES),wfsxy(2,MAXFRAMES),psxy(2,MAXFRAMES),
     -  turnpr(2,MAXFRAMES),scalardat(2,MAXFRAMES,MAXCOPY6),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      bdxy(2,nframe)=sin(bend/radtodeg)
      bdxy(1,nframe)=cos(bend/radtodeg)
            wbxy(2,nframe)=sin(wobble/radtodeg)
      wbxy(1,nframe)=cos(wobble/radtodeg)
      fsxy(2,nframe)=sin(faceshift/radtodeg)
      fsxy(1,nframe)=cos(faceshift/radtodeg)
      wfs=wobble-faceshift
      wfsxy(2,nframe)=sin(wfs/radtodeg)
      wfsxy(1,nframe)=cos(wfs/radtodeg)
      if (nmem .gt. 0) then
        psxy(2,nframe)=sin(psr/radtodeg)
        psxy(1,nframe)=cos(psr/radtodeg)
      end if
      return
      end
