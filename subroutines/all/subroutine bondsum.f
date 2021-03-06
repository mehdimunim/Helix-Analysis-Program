      subroutine bondsum(nhbtot,nhbfilt,indexbond,ixres,atnames,
     -  resnames,isegno,iresno,ifres,ihb,iout,mxbonds,maxrsd,
     -  maxrec)
c     Get the atom and residue bond sums
      dimension indexbond(mxbonds),ixres(maxrec),ihb(mxbonds),
     -  isegno(maxrec),iresno(maxrec),ifres(maxrec)
      character*8 atnames(maxrec),resnames(maxrsd)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (7*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),ihb_a(MAXBONDS),
     -  ihb_r(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
c     print *,'BONDSUM nhbtot,nhbfilt=',nhbtot,nhbfilt
      call askyn(
     -  'Do you want to limit the bond sums to the filtered bonds',56,
     -  1,-1,iusefilt,0,0)
      nhbuse=nhbtot
      if (iusefilt .eq. 1) nhbuse=nhbfilt
      maxat=0
      maxres=0
      do ibb=1,nhbuse
        ib=ibb
        if (iusefilt .eq. 1) ib=indexbond(ibb)
        ib1=ihbpair(1,ib)
        ib2=ihbpair(2,ib)
        if (maxat .lt. ib1) maxat=ib1
        if (maxat .lt. ib2) maxat=ib2
        ir1=ixres(ib1)
        ir2=ixres(ib2)
        if (maxres .lt. ir1) maxres=ir1
        if (maxres .lt. ir2) maxres=ir2
      end do
      if (maxat .gt. MAXBONDS) then
        write (6,2000) maxat,MAXBONDS
        return
      end if
      call zeroiti(nhb_atot,0,maxat)
      call zeroiti(nhb_rtot,0,maxres)
      lnam=0
      lres=0
      do ifr=1,nframetot
        call readbitc(ires(1,ifr),ihb,nhbtot,30,MAXITEMS)
c       Get the atom and residue bond sums
        call zeroiti(ihb_a,0,maxat)
        call zeroiti(ihb_r,0,maxres)
        do ibb=1,nhbuse
          ib=ibb
          if (iusefilt .eq. 1) ib=indexbond(ibb)
          if (ihb(ib) .eq. 1) then
            ib1=ihbpair(1,ib)
            ib2=ihbpair(2,ib)
            ihb_a(ib1)=1
            ihb_a(ib2)=1
            call lastchar(atnames(ib1),lc,8)
            if (lnam .lt. lc) lnam=lc
            call lastchar(atnames(ib2),lc,8)
            if (lnam .lt. lc) lnam=lc
            ir1=ixres(ib1)
            ir2=ixres(ib2)
            ihb_r(ir1)=1
            ihb_r(ir2)=1
            call lastchar(resnames(ir1),lc,8)
            if (lres .lt. lc) lres=lc
            call lastchar(resnames(ir2),lc,8)
            if (lres .lt. lc) lres=lc
          end if
        end do
        do ia=1,maxat
          nhb_atot(ia)=nhb_atot(ia)+ihb_a(ia)
        end do
        do ir=1,maxres
          nhb_rtot(ir)=nhb_rtot(ir)+ihb_r(ir)
        end do
      end do
      write (iout,2003) 'atom'
      do ia=1,maxat
        if (nhb_atot(ia) .gt. 0) write (iout,2001) ia,
     -    atnames(ia)(1:lnam),resnames(ixres(ia))(1:lres),
     -    iresno(ia),isegno(ia),nhb_atot(ia),
     -    100.0*float(nhb_atot(ia))/float(nframetot)
      end do
      write (iout,2003) 'residue'
      do ir=1,maxres
        if (nhb_rtot(ir) .gt. 0) write (iout,2002) ir,
     -    resnames(ir)(1:lres),isegno(ifres(ir)),iresno(ifres(ir)),
     -    nhb_rtot(ir),100.0*float(nhb_rtot(ir))/float(nframetot)
      end do
      return
2000  format(' ERROR: Number of atoms (',i7,') exceeds the value of ',
     -  'the parameter MAXBONDS ',/,8x,'(',i7,') - atom and residue ',
     -  'bond sum calculation will be skipped')
2001  format(' Atom',i6,1x,a,1x,a,' resno=',i5,' C/S:',i2,
     -  ' # of frames bonded=',i5,' % bonded=',f6.2)
2002  format(' Residue',i6,1x,a,' S/C:',i2,' resno=',i5,
     -  ' # of frames bonded=',i5,' % bonded=',
     -  f6.2)
2003  format(/,' === List of the number and percent of frames where ',
     -  'each ',a,' formed a bond')
      end
