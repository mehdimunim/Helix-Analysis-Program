      subroutine helixcomp(nslt,nres,calph0,dir0,perpvec0,start0,end0,
     -  cent0,rn0,camod,axfact,calph,dir,perpvec,start,end,cent,
     -  anglechange,anglechangeref,axtol,iw0,torsion,
     -  rotation,nrep,indexax,incrot,ireorienthx,idebughx,radtodeg,
     -  pi,cc2,nreshx,icaahx,ihx,nhxres,maxhx,maxat)
      dimension anglechange(maxhx),anglechangeref(maxhx),indexax(3),
     -  cc2(3,maxat),icaahx(maxhx)
      real*8 calph0(3,maxhx),dir0(3),start0(3),end0(3),cent0(3),rn0(3),
     -  perpvec0(3,maxhx),camod(3,maxhx),axfact(maxhx),calph(3,maxhx),
     -  dir(3),start(3),end(3),cent(3),perpvec(3,maxhx)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      real*8 org(3),dev(3),startg(3),endg(3),circ(3),rn(3),xx(3),yy(3),
     -  startw(3),endw(3),dirw(3),dirdotdir0,ddot,rms,devs(3),deve(3)
      dimension rot(3,3),corig(3,3),ccurr(3,3),shift(3),angles(3),
     -  anglesn(3)
      character*9 decidebend
c     Compute the rotation angle of a helix form the average angle between
c     the corresponding perpendiculars to the helix axis from the Calphas
      incrhx=(ihx-1)*nhxres
      iverbort=1
      if (nframe .gt. 10 .or. nrep .gt. 1) iverbort=0
c     print *,'HELIXC nframe,nrep,iverbort=',nframe,nrep,iverbort
      if (nrep .eq. 1 .and. nframe .eq. 1 .and. idebughx .gt. 0)
     -  write (iw0,1001) 'initial',cent0,
     -  ((calph0(k,ir),k=1,3),ir=1,nres)
      it=0
c     cc2 will be the possibly translated/overlaid frame
      call helixaxis(cc2,nslt,0,calph,dirw,startw,endw,cent,perpvec,
     -  camod,anglechangeref,circ,rn,axfact,axtol,rot,
     -  rms,helixlen,angles,decidebend,nup,ndown,nrun,nnear,rcirc,
     -  turnperres,anglesn,0,incrot,nrep,1,nreshx,icaahx,ihx,nhxres,
     -  idebughx,radtodeg,pi,MAXHX)
      call trnsfrd(start,startw,3)
      call trnsfrd(end,endw,3)
      call trnsfrd(dir,dirw,3)
      if (nrep .eq. 1 .and. idebughx .gt. 0) write (iw0,1001)
     -  'current',cent,((calph(k,ir),k=1,3),ir=1,nres)
c     Shift the helix so that the start is at the origin
      do ir=1,nres
        call dvdif(calph(1,ir),start,calph(1,ir))
      end do
      call dvdif(end,start,end)
      do k=1,3
        shift(k)=start(k)
      end do
      call zeroitd(start,3)
      call trnsfrd(startg,start,3)
      call trnsfrd(endg,end,3)
      if (nrep .eq. 1) then
        dirdotdir0=ddot(dir,dir0)
        rotation=dacoscheck(dirdotdir0,ccc,1,6,'HELIXCOMP')
        res(1,nframe,incrhx+5)=cos(rotation)
        res(2,nframe,incrhx+5)=sin(rotation)
        call dvdif(cent,cent0,dev)
        res(1,nframe,incrhx+8)=dsqrt(ddot(dev,dev))
        res(2,nframe,incrhx+8)=dev(indexax(1))
        res(1,nframe,incrhx+9)=dev(indexax(2))
        res(2,nframe,incrhx+9)=dev(indexax(3))
        call dvdif(startw,start0,devs)
        res(1,nframe,incrhx+10)=devs(indexax(2))
        res(2,nframe,incrhx+10)=devs(indexax(3))
        call dvdif(endw,end0,deve)
        res(1,nframe,incrhx+11)=deve(indexax(2))
        res(2,nframe,incrhx+11)=deve(indexax(3))
        if (idebughx .gt. 0 .and. dirdotdir0 .gt. 0.001)
     -   write (iw0,1001) 'current',cent,((calph(k,ir),k=1,3),ir=1,nres)
        if (ireorienthx .eq. 1) then
c         For better result, rotate dir onto dir0. The rotation axis is the
c         normal to the dir0-dir plane
          call dcross(dir0,dir,org)
          do k=1,3
            corig(k,1)=0.0
            ccurr(k,1)=0.0
            corig(k,2)=dir0(k)
            ccurr(k,2)=dir(k)
            corig(k,3)=org(k)
            ccurr(k,3)=org(k)
          end do
          call ormat(rot,ccurr,corig,3,iverbort,iw0)
          do ir=1,nres
            call dsmatvec(rot,calph(1,ir),calph(1,ir))
          end do
          call trnsfrd(dir,dir0,3)
          if (idebughx .gt. 0) write (iw0,1001) 'current reoriented',
     -      cent,((calph(k,ir),k=1,3),ir=1,nres)
        end if
        do ir=1,nres
          call calcperp(start,dir,calph(1,ir),org,perpvec(1,ir),it)
        end do
        nflip=0
100     rotav=0
        rotav2=0
        changemin=100.0
        changemax=0.0
        nflipprev=nflip
c       write (79,*) 'NRES-',nres
        do ir=1,nres
          call angcomp(perpvec0(1,ir),dir0,perpvec(1,ir),angchange)
c          write (79,9671) ir,angchange,(dir0(k),k=1,3),
c     -      (perpvec(k,ir),k=1,3),(perpvec0(k,ir),k=1,3)
c9671      format(' ir=',i3,' dang=',f10.4,' dir0=',3f10.5,/,
c     -      ' perpvec=',3f10.5,' perpvec0=',3f10.5)
          rotav=rotav+angchange
          rotav2=rotav2+angchange**2
          if (angchange .lt. changemin) changemin=angchange
          if (angchange .gt. changemax) changemax=angchange
          anglechange(ir)=angchange
        end do
        torsion=rotav/nres
        sd=sqrt(abs(rotav2/nres-torsion**2))
cd77        write (77,6533) nframe,torsion*180.0/pi,changemin*180.0/pi,
cd77     -    changemax*180.0/pi,(anglechange(ir)*180.0/pi,ir=1,nres)
cd776533    format(' Nframe=',i6,' avg,min,max=',3f10.3,/,(10f8.2))
        if (changemax-changemin .gt. pi/2.0) then
c         Likely to have some sign flips
          do ir=1,nres
            do k=1,3
              xx(k)=-perpvec(k,ir)
            end do
            call angcomp(perpvec0(1,ir),dir0,xx,angchange)
            if (abs(angchange-torsion) .lt.
     -        abs(anglechange(ir)-torsion)) then
              call trnsfrd(perpvec(1,ir),xx,3)
              nflip=nflip+1
            end if
          end do
          if (nflip .le. nres .and. nflip .gt. nflipprev) go to 100
        end if
        if (nflip .gt. 0) then
c         Recalculate TPR
          call calcturnperres(turnperres,nres,incrot,perpvec,dir,
     -      anglechangeref,1,pi,MAXHX)
          res(1,nframe,incrhx+6)=cos(turnperres)
          res(2,nframe,incrhx+6)=sin(turnperres)
        end if
        res(1,nframe,incrhx+4)=cos(torsion)
        res(2,nframe,incrhx+4)=sin(torsion)
        call dcross(dir0,rn0,xx)
        call dcross(rn0,xx,yy)
        do k=1,3
          rot(1,k)=yy(k)
          rot(2,k)=xx(k)
          rot(3,k)=rn0(k)
        end do
c       Keep the normal from oscillating 180 degrees
        if (ddot(rn,rn0) .lt. 0.d0) call dvmul(rn,-1.0d0,rn)
        call dsmatvec(rot,rn,xx)
        res(1,nframe,incrhx+12)=xx(1)
        res(2,nframe,incrhx+12)=xx(2)
        call printhelix(iw0,startw,endw,cent,rms,helixlen,dirw,angles,
     -    decidebend,nup,ndown,nrun,nnear,rcirc,turnperres,anglesn,ihx,
     -    radtodeg)
        rn_rn0=radtodeg*dacoscheck(ddot(rn,rn0),ccc,1,6,'HELIXNORMALS')
        write (iw0,1000) radtodeg*torsion,radtodeg*sd,radtodeg*rotation,
     -    rn_rn0,dev,dsqrt(ddot(dev,dev))
      end if
      if (nrep .lt. 0) return
      return
1000  format(' Rotation=',f9.2,' SD=',f7.2,
     -  ' Local tilt=',f9.2,' N/Nr angle:',f8.3,/,
     -  ' X,Y,Z displacements=',3f8.2,' Absolute displacement=',f8.2)
1001  format(' ',a,' helix cent:',3f12.3,/,' Alpha C:',(3f12.3))
      end
