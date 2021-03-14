      subroutine selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
     -  iselfanc,iqfsel2,nhbneig,ineig,nbfound,nbresfound,nmc,nhbdist,
     -  rhbdist,bondname,lbondname,ianc_anc,nosameseg,isegno,c,itemp,
     -  itempres,nbonds,ifail,iout,maxng,maxrec,mxbonds)
      dimension ixres(maxrec),ianchor(nanchor),indexa(maxrec),
     -  nhbneig(maxrec),ineig(maxng,maxrec),ianc_anc(mxbonds),
     -  nhbdist(mxbonds),rhbdist(mxbonds),isegno(maxrec),c(3,maxrec),
     -  itemp(mxbonds),itempres(mxbonds)
      character*(*) bondname
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
c     dimension ixx(3000)
c     ianchor: list of anchor atoms
c     indexa(ia): one if ia is an anchor atom;
c                 -1 if it is non-anchor atom but allowed by the charge filter
c     iqfsel2: If > 0, apply the charge filter to both atoms in the H bond
c     ianchor2: If > 0, both atoms forming the H bond have to be anchors
c     iselfanc: If = 0, exclude anchor-anchor hydrogen bonds
      nbitmax=30*MAXITEMS
c     write (iout,*) 'SELECTBOND nframe,nbresfound,ianchor2,iselfanc=',
c    -  nframe,nbresfound,ianchor2,iselfanc
      call trajlimtest(nframe,MAXFRAMES)
      call zeroiti(itemp,0,mxbonds)
      call zeroiti(itempres,0,mxbonds)
      ndistsum=0
      do i=1,nbfound
        ndistsum=ndistsum+nhbdist(i)
      end do
      nbonds=0
      do iia=1,nanchor
        iat=ianchor(iia)
        if (indexa(iat) .eq. 0) print *,'PROGRAM ERROR: ',iia,
     -    'th anchorlist=',iat,' is not marked as anchor'
        do in=1,nhbneig(iat)
          ian=ineig(maxng-in+1,iat)
          idoit=1
c         See if end is also an anchor
          if (ianchor2 .eq. 1) then
            if (indexa(ian) .lt. 1) idoit=0
          else if (iselfanc .eq. 0) then
            if (indexa(ian) .gt. 0) idoit=0
          else
            if (iqfsel2 .eq. 1 .and. indexa(ian) .eq. 0) idoit=0
          end if
          if (nosameseg .eq. 1 .and. isegno(ian) .eq. isegno(iat))
     -      idoit=0
          if (idoit .eq. 1 .and.
     -        (indexa(ian) .lt. 1 .or. ian .gt. iat)) then
            if (ian .gt. iat) then
              ib1=iat
              ib2=ian
            else
              ib2=iat
              ib1=ian
            end if
            nbonds=nbonds+1
c           write (iout,7943) ian,iat,indexa(ian),ib1,ib2
c7943       format(' ian,iat=',2i5,' indexa(ian)=',i4,' ib1,2=',2i5)
            ihb=1
            do while (ihb .le. nbfound .and.
     -        (ib1 .ne. ihbpair(1,ihb) .or. ib2 .ne. ihbpair(2,ihb)))
              ihb=ihb+1
            end do
            if (ihb .gt. nbfound) then
              nbfound=nbfound+1
              if (ihb .lt. MAXBONDS) then
                ihbpair(1,nbfound)=ib1
                ihbpair(2,nbfound)=ib2
                nhbdist(ihb)=1
                if (indexa(ian) .gt. 0) ianc_anc(ihb)=1
              else
                write (6,1005) bondname(1:lbondname),mxbonds
                write (iout,1005) bondname(1:lbondname),mxbonds
                if (nmc .gt. 0) then
                  write (6,1008) nmc
                  write (iout,1008) nmc
                else
                  percdone=100.0*float(iia)/float(nanchor)
                  write (6,1009) percdone
                  write (iout,1009) percdone
                end if
                call askyn('Do you want to continue without tracking',
     -            40,0,1,istopscan,112,0)
                ifail=2*istopscan-1
                return
              end if
            else
              nhbdist(ihb)=nhbdist(ihb)+1
            end if
            d12=sqrt(dist2(c(1,ib1),c(1,ib2)))
            rhbdist(ihb)=rhbdist(ihb)+d12
            itemp(ihb)=1
            ir1=ixres(ib1)
            ir2=ixres(ib2)
            ihb=1
            do while (ihb .le. nbresfound .and.
     -        (ir1 .ne. ihb_pair_res(1,ihb) .or.
     -         ir2 .ne. ihb_pair_res(2,ihb)))
              ihb=ihb+1
            end do
            if (ihb .gt. nbresfound) then
              nbresfound=nbresfound+1
              ihb_pair_res(1,nbresfound)=ir1
              ihb_pair_res(2,nbresfound)=ir2
              ihb_pair_res(3,nbresfound)=1
              ihb_pair_res(3,nbresfound+1)=0
              itempres(nbresfound)=1
            else if (itempres(ihb) .eq. 0) then
              ihb_pair_res(3,ihb)=ihb_pair_res(3,ihb)+1
              itempres(ihb)=1
            end if
          end if
        end do
      end do
      ndistsum=0
      itsum=0
      do i=1,nbfound
        ndistsum=ndistsum+nhbdist(i)
        itsum=itsum+itemp(i)
      end do
      if (nbfound .le. nbitmax)
     -  call savebitc(ires(1,nframe),itemp,nbfound,30,MAXITEMS)
      if (nbfound .eq. nbitmax+1) then
        write (6,1005) bondname(1:lbondname),nbitmax
        write (iout,1005) bondname(1:lbondname),nbitmax
      end if
      return
1005  format(' ERROR: maximum number of ',a,' bonds to store (',i5,
     -  ') is exceeded',/,8x,'Reduce the number of anchors or',/,
     -  8x,'increase the size of the array res in the common block ',
     -  '/analres/')
1008  format(8x,'Last frame processed=',i9)
1009  format(8x,'Processed ',f5.1,' % of the anchors')
      end
