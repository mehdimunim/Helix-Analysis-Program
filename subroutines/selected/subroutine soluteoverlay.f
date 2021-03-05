      subroutine soluteoverlay(isubcrm,ioverlay,nslt,nsegslt,c,cres,cc1,
     -  cc2,crmslt0,crmslt,atw,atw1,overlaysds,overlaysh,molsltlim,
     -  itemp,idebughx,iw0,maxrsd,maxat)
      dimension c(3,maxat),cres(3,maxat),cc1(3,maxat),cc2(3,maxat),
     -  crmslt0(3),crmslt(3),atw(nslt),atw1(nslt),molsltlim(3,maxrsd),
     -  itemp(maxat),overlaysds(maxat),overlaysh(maxat)
      dimension crmshift(3),com1(3),com2(3),rot(3,3)
c     print *,'SOLUTEOVERLAY nsegslt,nslt,maxat=',nsegslt,nslt,maxat
      call zeroit(crmshift,3)
      if (isubcrm+ioverlay .eq. 0) then
c       No translation or rotation
        call trnsfr(cc2,c,3*nslt)
      else
        call cofms(c,crmslt,nslt,atw)
        call arrdiff(crmslt,crmslt0,crmshift,3)
        if (ioverlay .eq. 0) then
c         Just translation
          do ia=1,nslt
            call arrdiff(c(1,ia),crmshift,cc2(1,ia),3)
          end do
        else if (ioverlay .eq. 1) then
c         Overlay of the whole solute
          call indexit(itemp,1,nslt,0)
          call bestoverlay(nslt,itemp,itemp,cres,c,atw,0.d0,
     -      cc1,cc2,atw1,rot,com1,com2,idebughx,0.001,iw0,maxat)
          call shiftmol(c,nslt,com2,cc2,-1.0)
          call rotate_c(cc2,nslt,rot,cc2,'HELIX',5)
          call shiftmol(cc2,nslt,com1,cc2,+1.0)
          overlaysd=sdsumix(nslt,cres,cc2,atw,0,itemp,devmax,maxat)
        else
c         Overlay separately each solute molecule
          do is=1,nsegslt
            nats=molsltlim(2,is)-molsltlim(1,is)+1
            ifat=molsltlim(1,is)
            call indexit(itemp,1,nats,0)
            call bestoverlay(nats,itemp,itemp,cres(1,ifat),c(1,ifat),
     -        atw(ifat),0.d0,cc1(1,ifat),cc2(1,ifat),atw1,rot,com1,com2,
     -        idebughx,0.001,iw0,maxat)
            call shiftmol(c(1,ifat),nats,com2,cc2(1,ifat),-1.0)
            call rotate_c(cc2(1,ifat),nats,rot,cc2(1,ifat),'HELIXs',6)
            call shiftmol(cc2(1,ifat),nats,com1,cc2(1,ifat),+1.0)
            overlaysds(is)=sdsumix(nats,cres(1,ifat),cc2(1,ifat),
     -        atw(ifat),0,itemp,devmax,maxat)
            overlaysh(is)=sqrt(dist2(com1,com2))
          end do
        end if
      end if
      if (ioverlay .eq. 0) then
        write (iw0,1004) crmshift
      else if (ioverlay .eq. 1) then
        write (iw0,1004) crmshift,' RMSD=',overlaysd
      else if (ioverlay .eq. 2) then
        write (iw0,1003)
     -   (overlaysh(is),sqrt(overlaysds(is)),is=1,nsegslt)
      else if (isubcrm .gt. 0) then
        write (iw0,1004) crmshift
      end if
      return
1003  format(' (Molecular shift, RMSD):',4(' (',f8.2,',',f8.2,')'))
1004  format(' Solute COM shift=',3f8.2,a,f8.2)
      end
