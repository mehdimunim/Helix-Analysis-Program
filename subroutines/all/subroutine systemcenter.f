      subroutine systemcenter(n,nmolslt,nmolsltnoion,molsltlim,c,ct,it1,
     -  atw,cell,ncell,cellalt,icellalt,ixyzexcld,ixyzincld,nslt,naslv,
     -  iacent,imcenter_r,imcenter,nconf,mxat,maxrsd)
      dimension c(3,mxat),ct(3,mxat),molsltlim(3,maxrsd),cell(3,ncell),
     -  cellalt(3,ncell),it1(mxat),atw(mxat)
      dimension c00(3),c0(3),cmin(3),cmax(3)
cDB      common /DEBUG/ cx0(3,10000),
cDB     -  cx1(3,10000),cx2(3,10000),cx3(3,10000),cx4(3,10000)
c     print *,'SYSTEMCENT n,nmolslt,nslt,naslv=',n,nmolslt,nslt,naslv
c     print *,'SYSTEMCENT nmolsltnoion,iacent,maxrsd,mxat=',
c    -  nmolsltnoion,iacent,maxrsd,mxat
c     print *,'SYSTEMCENT ncell,ixyzexcld,ixyzincld=',
c    -  ncell,ixyzexcld,ixyzincld
c     write (78,*) 'SYSTEMCENT'
c     write (78,8822) (im,(molsltlim(k,im),k=1,3),im=1,nmolslt)
c8822 format(' im=',i4,' molsltlim=',3i6)
c     Bring together the solute molecules
c     Find the center of the first solute molecule and shift it to <0,0,0>
cDB      call trnsfr(cx0,c,3*nslt)
      imcenter=imcenter_r
c     print *,'IACENT,IMCENTER=',iacent,imcenter
      if (iacent .eq. 0) then
        if (nconf .le. 1) then
          ndrop=0
          do ia=1,nslt
            if (atw(ia) .eq. 0.0) ndrop=ndrop+1
          end do
          print *,'Solute centering will use ',nslt-ndrop,' atoms'
        end if
        if (imcenter .eq. 0) then
          if (nmolslt .eq. 1) then
            imcenter=1
          else
c           Default center: largest solute molecule's center
            maxmolats=0
            imc_def=0
            do im=1,nmolslt
              nmem=molsltlim(2,im)-molsltlim(1,im)+1
              if (nmem .gt. maxmolats) then
                maxmolats=nmem
                imc_def=im
              end if
            end do
            call askyn(
     -      'Do you want to specify the solute molecule at the center',
     -        56,1,-1,imcenterset,0,0)
            if (imcenterset .eq. 1) then
              call getint('Solute molecule number',22,1,1,imc_def,
     -          imcenter,0)
            else
              imcenter=imc_def
              write (6,1003) imcenter
            end if
          end if
        end if
        if (molsltlim(3,imcenter) .eq. 0) then
          call extension(c,it1,0,molsltlim(1,imcenter),
     -      molsltlim(2,imcenter),cmin,cmax,c00,0,1,v)
          ncentspec=0
        else if (molsltlim(3,imcenter) .eq. -1) then
cxx       Calculate COM
          call zeroit(c00,3)
          atwsum=0.0
          do ia=molsltlim(1,imcenter),molsltlim(2,imcenter)
            do k=1,3
              c00(k)=c00(k)+atw(ia)*c(k,ia)
            end do
            atwsum=atwsum+atw(ia)
          end do
          if (atwsum .gt. 0.0) then
            do k=1,3
              c00(k)=c00(k)/atwsum
            end do
          end if
        else
          call trnsfr(c00,c(1,molsltlim(3,imcenter)),3)
          ncentspec=1
        end if
        if (nconf .le. 1 .and. nmolslt .gt. 1)
     -    write (6,1001) 'molecule',imcenter
      else
        if (iacent .gt. molsltlim(2,nmolslt)) then
          write (6,1002) iacent,molsltlim(2,nmolslt)
          stop
        end if
        call trnsfr(c00,c(1,iacent),3)
        imcenter=1
        do while (molsltlim(2,imcenter) .lt. iacent)
          imcenter=imcenter+1
        end do
        if (nconf .le. 1) write (6,1001) 'atom',iacent
      end if
      if (ixyzexcld .gt. 0) c00(ixyzexcld)=0.0
      if (ixyzincld .gt. 0) then
        do k=1,3
          if (k .ne. ixyzincld) c00(k)=0.0
        end do
      end if
      call shiftmol(c,n,c00,c,-1.0)
cDB      call trnsfr(cx1,c,3*nslt)
c     call savepdb(88,'MOLEC_IMCENTER_CENTERED.pdb',27,c,n,0)
      noshift=0
      nmolfst=1
      call molreset(nmolfst,nmolslt,nmolshift,c,ct,molsltlim,it1,
     -  cell,ncell,cellalt,icellalt,imcenter,maxrsd,mxat)
c     call savepdb(88,'AFTER_MOLRESET.pdb',18,c,n,0)
cDB      call trnsfr(cx2,c,3*nslt)
c     print *,'After  molreset'
c     do is=nmolfst,nmolslt
c       write (6,*) 'molsltlim=',molsltlim(3,is)
c       call extension(c,it1,0,molsltlim(1,is),molsltlim(2,is),
c    -      cmin,cmax,c0,1,0,v)
c     end do
c     Repeat, to bring in 2nd neighbor cell members
c     call molreset(nmolfst,nmolslt,nmolshift,c,ct,molsltlim,it1,
c    -  cell,ncell,cellalt,icellalt,imcenter,maxrsd,mxat)
      if (noshift .eq. 1 .and. nmolshift .gt. 0 .and. nconf .le. 10)
     -  print *,'WARNING: ',nmolshift,' centered molecules were ',
     -    'shifted again after recentering'
      nsltnoion=molsltlim(2,nmolsltnoion)
c     call extension(c,it1,0,1,nsltnoion,cmin,cmax,c0,0,0,v)
      if (iacent .eq. 0 .and. noshift .eq. 0 .and. nmolslt .gt. 1) then
        if (nsltnoion .gt. 0) then
          call cofms(c,c0,nsltnoion,atw)
          call shiftmol(c,n,c0,c,-1.0)
          call arrsum(c0,c00,c0,3)
          if (nconf .le. 1) write (6,1000) c0
          nmolfst=1
          noshift=1
          call ionreset(nmolfst,nmolslt,nmolshift,c,molsltlim,it1,atw,
     -      cell,ncell,ixyzexcld,ixyzincld,maxrsd,mxat)
        end if
      end if
cDB      call trnsfr(cx3,c,3*nslt)
c     Check for ions outside the cell
      nreset=0
      do is=nmolsltnoion+1,nmolslt
        call pbcreset(c(1,molsltlim(1,is)),
     -    molsltlim(2,is)-molsltlim(1,is)+1,c(1,molsltlim(1,is)),
     -    cell,ncell,ixyzexcld,ixyzincld,img)
        if (img .gt. 1) then
          nreset=nreset+1
          call pbcreset(c(1,molsltlim(1,is)),
     -      molsltlim(2,is)-molsltlim(1,is)+1,c(1,molsltlim(1,is)),
     -      cell,ncell,ixyzexcld,ixyzincld,img)
        end if
      end do
cDB      call trnsfr(cx4,c,3*nslt)
c     Reset solvents
      nsw=(n-nslt)/naslv
      do iw=1,nsw
        call extension(c,it1,0,nslt+(iw-1)*naslv+1,nslt+iw*naslv,
     -    cmin,cmax,c0,0,0,v)
cx      call pbcreset(c(1,nslt+(iw-1)*naslv+1),naslv,c0,
cx   -    cell,ncell,ixyzexcld,ixyzincld,img)
        call genimdist(c0,cell,1,ncell,img,d2)
        if (img .gt. 1)
     -    call shiftmol(c(1,nslt+(iw-1)*naslv+1),naslv,cell(1,img),
     -      c(1,nslt+(iw-1)*naslv+1),-1.0)
      end do
      return
1000  format(' Aggregated solute molecules are shifted by ',3f8.3,
     -  ' A')
1001  format(' Solute is centered around solute ',a,i6)
1002  format(' ERROR: solute center requested (',i8,') is outside the ',
     -  'solute atom range (1 - ',i8,')')
1003  format(' Largest solute molecule (',i6,') selected for center')
      end
