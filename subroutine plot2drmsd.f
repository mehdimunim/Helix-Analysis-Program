      subroutine plot2drmsd(nrep,iw0,iw1,xtraj,maxn,title,title2,
     -  ltitle2,ltrajfile2,xlab,lxlab,ncolcode,maxcolcode,iedit,noopt2d,
     -  limresrange,normsavg,rmsdmin,rmsdmax,absdevmin,absdevmax,
     -  rmsdmn,rmsdmx,indexa,indexr,ixshuffle,ixshuffleref,ilastclx,
     -  nclx,ilastcly,ncly,ym_2d,itemp1,itemp2,temp,matplotonly,isymm,
     -  noplotdist,iskip2dplot,isortmat,noclose,ipspage,ifindbest)
      character*80 title
      character*(*) title2
      character*(*) xlab
      dimension xtraj(maxn),indexa(maxn),indexr(maxn),ixshuffle(maxn),
     -  ixshuffleref(maxn),ilastcly(maxn),ilastclx(maxn),
     -  itemp1(maxn),itemp2(maxn),temp(maxn)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      dimension kc(1,1)
      real*8 dc(1,1)
      dimension nrmsd(100),xrmsd(100),refrmsd(100)
      character*4 yclab(1)
      character*80 title2add
      data nyclab /1/,lyclab /1/
c     print *,'PLOT2DRMSD isymm,nframe,nframeref=',
c    -                   isymm,nframe,nframeref
c     print *,'PLOT2DRMSD noclose,matplotonly=',noclose,matplotonly
c     print *,'PLOT2DRMSD ioclose,matplotonly=',noclose,matplotonly
      cdiammin=100000.0
      cdiammax=0.0
      if (ifindbest .gt. 0) then
        call indexit(indexa,1,nframe,0)
        write (6,1015)
        call findbestrep(iw0,0,0,nframe,indexa,irdavmn,irdmxmn,
     -    0,0.0,0.0,0.0,cdiam,cdiammin,cdiammax,0,'RMSD',4,6,MAX2D)
      end if
      rcmax=0.0
      rcmin=100000.0
      nz_rmsd=0
      do ix=1,nframe
        if (isymm .eq. 1) then
          do iy=ix+1,nframe
            if (rmsd2d(ix,iy) .gt. rcmax) rcmax=rmsd2d(ix,iy)
            if (rmsd2d(ix,iy) .lt. rcmin) rcmin=rmsd2d(ix,iy)
            if (rmsd2d(ix,iy) .eq. 0.0) nz_rmsd=nz_rmsd+1
          end do
          rmsd2d(ix,ix)=0.0
        else
          do iy=1,nframeref
            if (rmsd2d(ix,iy) .gt. rcmax) rcmax=rmsd2d(ix,iy)
            if (rmsd2d(ix,iy) .lt. rcmin) rcmin=rmsd2d(ix,iy)
            if (rmsd2d(ix,iy) .eq. 0.0) nz_rmsd=nz_rmsd+1
          end do
        end if
      end do
      write (6,1011) 'RMSD',rcmin,rcmax
      write (iw0,1011) 'RMSD',rcmin,rcmax
      if (rcmin .eq. 0.0) then
        write (6,1017) 'WARNING','minimum','some',' ',nz_rmsd
        write (iw0,1017) 'WARNING','minimum','some',' ',nz_rmsd
      end if
      if (rcmax .eq. 0.0) then
        write (6,1017) 'ERROR','maximum','all'
        write (iw0,1017) 'ERROR','maximum','all'
        print *,'Exiting plotting'
        return
      end if
      if (absdevmin .lt. 100000.0) then
        write (6,1011) 'absolute deviation',absdevmin,absdevmax
        write (iw0,1011) 'absolute deviation',absdevmin,absdevmax
      end if
      if (isymm .eq. 1) then
        rmsdmn=0.0
      else if (rmsdmin .gt. 0.0) then
        write (6,1010) 'minimum',rmsdmin
        write (iw0,1010) 'minimum',rmsdmin
        rmsdmn=rmsdmin
      else
        rmsdmn=rcmin
      end if
      if (rmsdmax .ne. 0.0) then
        write (6,1010) 'maximum',rmsdmax
        write (iw0,1010) 'maximum',rmsdmax
        rmsdmx=rmsdmax
      else
        rmsdmx=rcmax
      end if
      nframetop=max0(nframe,nframeref)
      if (iskip2dplot .eq. 0) then
        inc=max0(1,500/nframetop)
        scalefac=amin1(1.0,500.0/float(nframetop))
        iyrange=scalefac*(nframeref*inc)+15
        iytop=ym_2d*0.87
        iydel=max0(150,iytop-iyrange)
c       print *,'IYDEL,IYTOP,IYRANGE,NFRAMEREF=',
c    -    iydel,iytop,iyrange,nframeref
        iouttemp=91
        call contractmat(rmsd2d,nframe,nframeref,nframeplot,
     -    nframerefplot,navg,ixshuffle,ixshuffleref,itemp1,itemp2,
     -    temp,isortmat,iouttemp,MAX2D)
c       iymax=iydel+iyrange
c       iymax=iymax+15
        iydel=iydel-15
        if (limresrange .gt. 0) then
          write (iw1,1005) 25,iytop
          call psshow(iw1,
     -      'RMSD was calculated on a limited residue set',44)
          iydel=iydel-15
          iytop=iytop-15
        end if
        write (iw1,1005) 25,iytop
        write (title2add,1019) navg
        iydel=iydel-15
        iytop=iytop-15
        ladd=30
        if (ltitle2 .gt. 0) then
          title2add(1:ladd+ltitle2)=title2(1:ltitle2)//' '//
     -      title2add(1:ladd)
          ladd=ladd+ltitle2
        end if
        call psshow(iw1,title2add(1:ladd),ladd)
        iydel=iydel-15
        iytop=iytop-15
        itrajname=1
        if (ltrajfile2 .gt. 0) then
          nframe2=nframe
          iydel=iydel-15
          iytop=iytop-15
          itrajname=3
        end if
        if (isortmat .eq. 1 .and. navg .gt. 1) then
          call indexit(ixshuffle,1,nframe,0)
          call indexit(ixshuffleref,1,nframeref,0)
        end if
        call plotmat(iw1,kc,rmsd2d,dc,nframeplot,nframerefplot,0,0,0,0,
     -    navg,nrep,25,iydel,00,iytop,rmsdmn,rmsdmx,ncolcode,maxcolcode,
     -    ixdelsh,iydelsh,inc,scalefac,indexr,ixshuffle,ixshuffleref,
     -    title,0,' ',0,itrajname,xlab,lxlab,xtraj,yclab,nyclab,lyclab,
     -    1,MAX2D,1,MAX2D,MAX2D,ipspage,0)
        if (navg .gt. 1) then
          rewind iouttemp
          read (iouttemp,end=999)
     -      ((rmsd2d(i,j),i=1,nframe),j=1,nframeref)
          close (iouttemp,status='delete')
        end if
        if (nframe .ge. 50) iydel=iydel-50
        if (nframe .lt. 50) iydel=iydel-40
        ixcent=amax1(0.0,(scalefac*(nframe*inc)-80*ncolcode)/2)
        call colcodeminmax(iw1,20+ixcent,-iydel,nrep,ncolcode,
     -    maxcolcode,rmsdmn,rmsdmx)
        iydel=iydel-40
        write (iw1,1005) 50,iydel
        write (iw1,1002) rmsdmn,rmsdmx
        write (iw1,1005) 50+250,iydel
        if (noopt2d .eq. 0) then
          if (iedit .eq. 0) write (iw1,1009) 'all atoms'
          if (iedit .eq. 1) write (iw1,1009) 'a limited set of atoms'
        else
          write (iw1,1016)
        end if
c       Draw lines in the matrix delineating the clusters
        print *
        if (max0(nclx,ncly) .le. 10 .and. min0(nclx,ncly) .gt. 0) then
          write (iw1,1018)
c          print *,'NCLX,NCLY=',nclx,ncly,' NFRAME=',nframe
c          write (6,9876) (ilastclx(i),i=1,nclx)
c9876      format('ILASTCLX:',/,(20i4))
          call rgbcolor(iw1,9)
          lw=2
          if (max0(nclx,ncly) .gt. 10) lw=1
          lw=lw/scalefac
          write (iw1,1012) lw
          if (max0(nclx,ncly) .gt. 1 .and. scalefac .ne. 1.0)
     -      write (iw1,1013) scalefac,scalefac
          write (iw1,*) 'np'
c         ixdell=ixdelsh
          ixdell=25
          if (nframeplot*scalefac .lt. 550) ixdell=ixdell+40/scalefac
c    -      ixdell=ixdelsh+40/scalefac
          do ic=1,nclx
            if (ilastclx(ic) .ne. nframe) then
              write (iw1,1005) ixdell+ilastclx(ic)*inc,iydelsh
              write (iw1,1006) ixdell+ilastclx(ic)*inc,
     -          iydelsh+nframeref*inc
            end if
          end do
          do ic=1,ncly
            if (ilastcly(ic) .ne. nframe) then
              write (iw1,1005) ixdell,iydelsh+ilastcly(ic)*inc
              write (iw1,1006) ixdell+nframe*inc,
     -          iydelsh+ilastcly(ic)*inc
            end if
          end do
          write (iw1,*) 'sk'
          if (max0(nclx,ncly) .gt. 1 .and. scalefac .ne. 1.0)
     -        write (iw1,1013) 1.0/scalefac,1.0/scalefac
        else
          if (max0(nclx,ncly) .gt. 0) write (6,1004)
        end if
        write (iw1,*) 'showpage'
        if (matplotonly .eq. 1) then
          if (noclose .eq. 0) close (iw1)
          return
        end if
      end if
      if (normsavg .eq. 0 .or. nframe .lt. 25) then
c       2D RMSD map, dont calculate normalization
        write (iw0,1000)
        do il=1,nframe-1
          rmssum=0.0
          rmssum2=0.0
          rmssum4=0.0
          do i=1,nframe-il
            rmssum=rmssum+rmsd2d(i,i+il)
            rmssum2=rmssum2+rmsd2d(i,i+il)**2
            rmssum4=rmssum4+rmsd2d(i,i+il)**4
          end do
          res(1,il,8)=rmssum/(nframe-il)
          res(2,il,8)=rmssum2/(nframe-il)
          res(1,il,9)=sqrt(abs(rmssum2/(nframe-il)-res(1,il,8)**2))
          res(2,il,9)=sqrt(abs(rmssum4/(nframe-il)-res(2,il,8)**2))
          write (iw0,1001) il,((res(k,il,ii),ii=8,9),k=1,2)
        end do
        if (noplotdist .eq. 0) then
          call arminmax2(res(1,1,8),1,nframe-1,2,armin1,armax1,armin2,
     -      armax2,0,2)
          call roundlim(armax1,y1div,ny1div)
          call roundlim(armax2,y2div,ny2div)
          call plot2fun(iw1,2,xtraj,res(1,1,8),res(1,1,9),nframe-1,
     -      0.0,0.0,0,0.0,y1div,ny1div,0.0,y2div,ny2div,title,80,' ',1,
     -      xlab,lxlab,'2D-map averaged RMSD',20,
     -      '2D-map averaged MSD',19,1,0,6,2,0,0,1,0,0,
     -      ipspage,1,1,0)
        end if
      end if
      call zeroiti(nrmsd,0,100)
      if (normsavg .eq. 0 .and. nframe .ge. 25) then
        rcmaxx=rcmax**2
        ircmax=rcmax**2
        rcmax=ircmax+1
        if (nframe .ge. 100) then
          rcdiv=rcmaxx/99.9999
        else if (nframe .ge. 50) then
          rcdiv=rcmaxx/49.9999
        else
          rcdiv=rcmaxx/24.9999
        end if
        do il=1,nframe-1
          do i=1,nframe-il
            ix=rmsd2d(i,i+il)**2/rcdiv+1
            if (ix .gt.100) ix=100
            nrmsd(ix)=nrmsd(ix)+1
          end do
        end do
        pmax1=0.0
        pmax2=0.0
c       dcoef=rcmax/sqrt(float(nframe))
        dcoef=rcmaxx/float(nframe)
        do ix=1,100
          refrmsd(ix)=(nframe-(ix*rcdiv)/dcoef)/
     -      float((nframe*(nframe-1))/2)
        end do
        write (iw0,1007) ' '
        do ix=1,100
          res(1,ix,9)=float(nrmsd(ix))/float((nframe*(nframe-1))/2)
          if (pmax1 .lt. res(1,ix,9)) pmax1=res(1,ix,9)
          xrmsd(ix)=ix*rcdiv
          if (ix .gt. 98) then
            res(2,ix,9)=res(2,98,9)
          else
            res(2,ix,9)=res(1,ix,9)/refrmsd(ix)
          end if
          if (pmax2 .lt. res(2,ix,9)) pmax2=res(2,ix,9)
          write (iw0,1008) ix,xrmsd(ix),res(1,ix,9),' ',res(2,ix,9)
        end do
        if (noplotdist .eq. 0) then
          call roundlim(pmax1,y1div,ny1div)
          call roundlim(pmax2,y2div,ny2div)
          call plot2fun(iw1,2,xrmsd,res(1,1,9),res(1,1,9),100,0.0,
     -      rcmax/10.,10,0.0,y1div,ny1div,0.0,y2div,ny2div,title,80,
     -      ' ',1,'MSD',3,' MSD distribution',17,
     -      ' MSD/Reference ideal distribution',33,1,
     -      0,6,2,1,0,1,0,0,ipspage,1,1,0)
        end if
      else if (normsavg .eq. 1) then
c       Cross RMSD map
        ircmax=rcmax+1.0
        rcmax=ircmax
        rcdiv=rcmax/99.9999
        do i1=1,nframe
          do i2=1,nframeref
            ix=rmsd2d(i2,i1)/rcdiv+1
            if (ix .gt.100) ix=100
            nrmsd(ix)=nrmsd(ix)+1
          end do
        end do
        pmax1=0.0
        write (iw0,1007) ' R'
        do ix=1,100
          res(1,ix,9)=float(nrmsd(ix))/float(nframe*nframeref)
          if (pmax1 .lt. res(1,ix,9)) pmax1=res(1,ix,9)
          xrmsd(ix)=ix*rcdiv
          write (iw0,1008) ix,xrmsd(ix),res(1,ix,9)
        end do
        call roundlim(pmax1,y1div,ny1div)
        call plot2fun(iw1,1,xrmsd,res(1,1,9),res(1,1,9),100,0.0,
     -    rcmax/10.,10,0.0,y1div,ny1div,0.0,y2div,ny2div,title,80,' ',1,
     -    'RMSD',4,'RMSD distribution',17,' ',1,1,
     -    0,6,2,1,0,1,0,0,ipspage,1,1,0)
      end if
      if (noclose .eq. 0) close (iw1)
      return
999   write (6,1014)
      stop
      return
1000  format(/,' Average RMSD & MSD as a function of time (frame ',
     -  '# difference) over the 2D map:',/)
1001  format(' RMSD(',i5,')=',f10.5,' SD=',f10.5,
     -  ' MSD=',f10.5,' SD=',f10.5)
1002  format('( Range of the RMSD scale:',f8.3,'  -',f8.3,') show')

1004  format(' NOTE: There are too many clusters for drawing ',
     -  'cluster-delineating lines')
1005  format(2i5,' m')
1006  format(2i5,' l')
1007  format(/,'Distribution of the',a,'MSD values:',/)
1008  format(i4,' p(RMSD=',f9.2,')=',f8.6,a,'p(refRMSD)=',f10.5)
1009  format('( Overlay is based on ',a,') show')
1010  format(' RMSD ',a,' used for plotting=',f10.2,' A')
1011  format(' Range of ',a,' values: [',f12.5,',',f12.5,'] A')
1012  format(i2,' lw')
1013  format(f10.6,f10.6,' scale')
1014  format(' PROGRAM ERROR: could not restore original matrix from ',
     -  'the temporary file rij.tmp')
1015  format(' Details of the structure ensemble as a single cluster:')
1016  format('( RMSD was calculated without overlay) show')
1017  format(1x,a,' RMSD ',a,' is zero - ',a,' items are identical',a,/,
     -  ' Number of zero RMSDs=',i10)
1018  format('% Draw cluster separator lines')
1019  format('Number of frames averaged=',i2)
      end
