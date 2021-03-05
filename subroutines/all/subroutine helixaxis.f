      subroutine helixaxis(c,nslt,iw0,calph,axisdir,axisini,axisend,
     -  helixcent,perpvec,camod,anglechangeref,circ,rn,axfact,axtol,rot,
     -  rms,helixlen,angles,decidebend,nup,ndown,nrun,nnear,rcirc,
     -  turnperres,anglesn,irot,incrot,nrep,irefang,nreshx,icaahx,ihx,
     -  nhxres,idebughx,radtodeg,pi,maxhlx)
      dimension c(3,nslt),rot(3,3),anglechangeref(maxhlx),icaahx(maxhlx)
      real*8 camod(3,maxhlx),axfact(maxhlx),perpvec(3,maxhlx),
     -  calph(3,maxhlx),axisdir(3),axisini(3),axisend(3),helixcent(3),
     -  circ(3),rn(3)
      real*8 rms
      dimension angles(3),anglesn(3)
      character*9 decidebend
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      character*60 message
      real*8 ddistsq
      character*9 decide(6)
      data decide/'Too short','Bent     ','Random   ',
     -  'Alternate','Too long ','         '/
c     print *,'HELIXAXIS maxhlx,nreshx,nrep=',maxhlx,nreshx,nrep
      iprintkahn=0
      it=0
      do ir=1,nreshx
        do k=1,3
          calph(k,ir)=c(k,icaahx(ir))
        end do
      end do
      if (irot .eq. 1) then
        do ir=1,nreshx
          call dsmatvec(rot,calph(1,ir),calph(1,ir))
        end do
      end if
      message=
     - 'Helix                                                       '
      call kahn(calph,nreshx,.true.,axisdir,axisini,axisend,rms,
     -  iprintkahn, message,maxhlx)
      rmsa=rms
      if (idebughx .gt. 0) then
        write (77,7011) 'axis direction',axisdir
        write (77,7011) 'initial position',axisini
        write (77,7011) 'end position',axisend
        write (77,*) 'Helix quality rms=',rms
        write (77,7012) nreshx,(icaahx(ir),ir=1,nreshx)
      end if
      if (nrep .le. 1) then
        do ir=1,nreshx
          call calcperp(axisini,axisdir,calph(1,ir),camod(1,ir),
     -      perpvec(1,ir),it)
        end do
        call checkbend(calph,axisdir,camod,axfact,perpvec,helixcent,
     -    nreshx,nup,ndown,nrun,nnear,axtol,rcirc,circ,rn,ibtyp,
     -    idebughx)
      end if
      if (nrep .le. 1) then
        do k=1,3
          anglesn(k)=dacoscheck(rn(k),ccc,1,6,'HELIXNORMAL')
        end do
        helixlen=dsqrt(ddistsq(axisini,axisend))
        do k=1,3
          angles(k)=dacoscheck(axisdir(k),ccc,1,6,'HELIXAXIS')
        end do
        call calcturnperres(turnperres,nreshx,incrot,perpvec,axisdir,
     -    anglechangeref,irefang,pi,maxhlx)
        nframes=max0(1,nframe)
        incr=(ihx-1)*nhxres
        call trajlimtest(nframe,MAXFRAMES)
        res(1,nframes,incr+7)=helixlen
        res(2,nframes,incr+7)=rcirc
        do k=1,3
          res(1,nframes,incr+k)=cos(angles(k))
          res(2,nframes,incr+k)=sin(angles(k))
        end do
        res(1,nframes,incr+6)=cos(turnperres)
        res(2,nframes,incr+6)=sin(turnperres)
        res(1,nframes,incr+13)=axisdir(1)
        res(2,nframes,incr+13)=axisdir(2)
        res(1,nframes,incr+14)=axisdir(3)
        res(2,nframes,incr+14)=helixcent(1)
        res(1,nframes,incr+15)=helixcent(2)
        res(2,nframes,incr+15)=helixcent(3)
        decidebend=decide(ibtyp)
        if (iw0 .gt. 0) call printhelix(iw0,axisini,axisend,helixcent,
     -    rms,helixlen,axisdir,angles,decidebend,nup,ndown,nrun,nnear,
     -    rcirc,turnperres,anglesn,ihx,radtodeg)
      end if
      return
7011  format(' Helix ',a,'=',3f10.4)
7012  format(' # of HX res=',i2,(' icaa:',20i5))
      end
