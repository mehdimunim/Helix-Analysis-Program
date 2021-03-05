      subroutine plotbond(ips,nhbtot,nhbplot,nbresplot,nclust,ianc_anc,
     -  ifhbclust,ilhbclust,index,ihb,ihbon,ihbstart,rdclust,iclust,
     -  xtrajmax,nydd,itit,ntit,xlab,nxlab,ylab,nylab,nohead,percmin,
     -  percmax,numres,minresdist,maxresdist,distmax,bondname,lbondname,
     -  nobullet,nhneigmin,nbondmax,naabondmax,icorrtyp,icorrtrans,
     -  nframeav,correxp,hblimfac,angmin,lablim,ibondtype,ibondlab,
     -  itemp1,ixresno,irrix,ixres,nrrbond,atnames,
     -  resnames,npspages,ipspage,mxbonds,maxrsd,maxrec)
c*****Plot the time-course of hydrogen, hydrophobic bonds & salt bridges
      dimension ianc_anc(mxbonds),ifhbclust(mxbonds),ilhbclust(mxbonds),
     -  index(mxbonds),ihb(mxbonds),ihbon(mxbonds),ihbstart(mxbonds),
     -  itemp1(mxbonds),ixresno(maxrsd),irrix(maxrsd),
     -  ixres(maxrec),nrrbond(mxbonds)
      character*(*) itit,xlab,ylab,bondname
      character*8 atnames(maxrec),resnames(maxrsd)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      character*41 clstyp
      common /cluster_typ/ nclstyp,inumclst(9),ireadcutoff(9),
     -  lclstyp(9),clstyp(9)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      dimension icolrb(2),lcorrtyp(3)
      character*1 corrtrans(3)
      character*36 corrtyp(3)
      character*37 bondlab(1)
      data icolrb /1,6/,lcorrtyp /36,28,28/
      data corrtrans /'-','/',' '/,corrtyp /
     -  'Correlate both the on and off states',
     -  'Correlate only the on states        ',
     -  'Correlate averaged on states        '/
c     nxd, nyd: number of tics in the x, y axes, resp
c     itit: array containing the title; ntit: number of characters in itit
c     xlab, ylab: labels of the x, y axes
c     nxlab, ylab: number of characters in the labels of the x,y axes
c     ibondlab=1: bond label; ibondlab=2: res pair label
      nhbtot=min0(nhbtot,30*MAXITEMS)
c     write (06,*) 'PLOTBOND nydd,nhbplot,nhbtot,nbresplot=',
c    -  nydd,nhbplot,nhbtot,nbresplot
c     print *,'PLOTBOND nclust,rdclust,numres,maxresdist=',
c    -                  nclust,rdclust,numres,maxresdist
c     print *,'PLOTBOND nframetot=',nframetot,' xtrajmax=',xtrajmax
c      write (40,8711) (ihbtores(i),i=1,nhbtot)
c8711  format(' IHBTORES:',/,(15i5))
      if (iclust .gt. 0 .and. nclust .eq. 0) return
      ntotplot=nhbplot
      if (ibondlab .eq. 2) ntotplot=nbresplot
      if (ntotplot .le. 20) then
        xm=650
        ym=amax1(200.0,350.0*float(ntotplot)/20.0)
        xm0=75
        xm00=-40
        xmm=0.9*xm
      else
        xm=675
        ym=495
        xmm=0.9*xm
        xm0=0.075*xm
        xm00=xm0-0.15*xmm
      end if
      ymm09=0.9*ym
      ymm=amin1(ymm09,float(nhbplot*10))
      if (ibondlab .eq. 1) then
        xm0=xm0+0.05*xmm
        xm00=xm0-0.20*xmm
      end if
      yshiftinc=12
      ym0=2.0*yshiftinc+5.0
      if (xlab(1:1) .eq. ' ') then
        nfrtot=nframetot
        call roundlimint(nfrtot,ixd,nxd)
      else
c       Round up nframetot based on xtrajmax
        frtot=xtrajmax
        call roundlim(frtot,xdiv,nxd)
        nfrtot=float(nframetot)*frtot/xtrajmax+1
        ixd=nfrtot/nxd+1
      end if
c     xfac=xmm/float(nxd*ixd)
      xfac=xmm/(float(nxd)*xdiv*float(nframetot)/xtrajmax)
      nyd=nydd
      if (ibondlab .eq. 2) then
        if (iclust .eq. 0) then
          call trnsfi(itemp1,index,nhbtot)
          call indexit(index,1,nhbtot,0)
        end if
      end if
      call roundlimint(ntotplot,iyd,nyd)
      ymax=nyd*iyd
      yfac=ymm/float(iyd*nyd)
      lw=1
      if (ntotplot .le. 100) lw=2
      if (ntotplot .le. 50) lw=3
      call zeroiti(ihbon,0,nhbtot)
      call zeroiti(ihbstart,0,nhbtot)
      if (iclust .eq. 0) call zeroiti(nrrbond,0,ntotplot)
      nbondmax=0
      naabondmax=0
      nhbsumsum=0
      do ifr=1,nframetot
        call readbitc(ires(1,ifr),ihb,nhbtot,30,MAXITEMS)
        if (ibondlab .lt. 2) then
          nhbsum=0
          nhbaasum=0
          do ibx=1,nhbtot
c           ibx=index(ib)
            if (ihb(ibx) .eq. 1) then
              nhbsum=nhbsum+1
              if (ianc_anc(ibx) .eq. 1) nhbaasum=nhbaasum+1
            end if
          end do
          scres(1,ifr)=nhbsum
          scres(2,ifr)=nhbaasum
          if (nbondmax .lt. nhbsum) nbondmax=nhbsum
          if (naabondmax .lt. nhbaasum) naabondmax=nhbaasum
          nhbsumsum=nhbsumsum+nhbsum
        else
          if (iclust .eq. 0) then
            do ib=1,ntotplot
              nrrbond(ib)=nrrbond(ib)+ihb(ib)
            end do
          end if
        end if
        if (ifr .eq. 1) then
          call trnsfi(ihbon,ihb,ntotplot)
          if (nohead .eq. 0) call psheader(ips,itit,ntit,
     -      -30,-130,830,830,npspages,ipspage)
          write (ips,3000) 10
          ipspage=ipspage+1
          write (ips,1019) ipspage
          write (ips,*) '-90 rotate'
c         write (ips,1001) -xm*1.1,0.1*ym,' translate'
          write (ips,1001) -xm*1.1,0.03*ym,' translate'
c         Drawing bounding box
          write (ips,*) 'np'
          write (ips,1001) xm0,ym0,' m'
          write (ips,1001) xm0,ym0+ymm+4,' l'
          write (ips,1001) xm0+xmm,ym0+ymm+4,' l'
          write (ips,1001) xm0+xmm,ym0,'   l'
          write (ips,1001) xm0,ym0,' l'
          write (ips,*) 'sk'
          write (ips,*) 'np'
          yshift=yshiftinc
          if (distmax .gt. 0.0 .or. nhneigmin .gt. 0) then
            write (ips,1001) xm0,ym0+ymm+yshift,' m'
            write (ips,1005) bondname(1:lbondname),distmax
            yshift=yshift+yshiftinc
          end if
          if (ibondtype .eq. 1 .or. ibondtype .eq. 2) then
            write (ips,1001) xm0,ym0+ymm+yshift,' m'
            if (ibondtype .eq. 2) write (ips,1004) nhneigmin
            if (ibondtype .eq. 1) write (ips,1015) hblimfac,angmin
            yshift=yshift+yshiftinc
          end if
          if (percmin .gt. 0.0 .or. percmax .lt. 100.0) then
            write (ips,1001) xm0,ym0+ymm+yshift,' m'
            write (ips,1007) percmin,percmax
            yshift=yshift+yshiftinc
          end if
          if (minresdist .gt. 0 .or. maxresdist .lt. numres-1) then
            write (ips,1001) xm0,ym0+ymm+yshift,' m'
            write (ips,1008) minresdist,maxresdist
            yshift=yshift+yshiftinc
          end if
          if (ibondlab .eq. 2) then
            write (ips,1001) xm0,ym0+ymm+yshift,' m'
            write (ips,1018) bondname(1:lbondname)
            yshift=yshift+yshiftinc
          end if
          write (ips,1001) xm0,ym0+ymm+yshift,' m'
          if (iclust .eq. 0) then
            write (ips,1009) bondname(1:lbondname)
          else
            if (rdclust .gt. 0.0) then
              write (ips,1010) clstyp(iclust)(1:lclstyp(iclust)),rdclust
            else
              write (ips,1012) clstyp(iclust)(1:lclstyp(iclust)),nclust
            end if
            yshift=yshift+yshiftinc
            write (ips,1001) xm0,ym0+ymm+yshift,' m'
            write (ips,1013) corrtyp(icorrtyp)(1:lcorrtyp(icorrtyp))
            if (icorrtyp .eq. 3) write (ips,1014) nframeav
            write (ips,1023) corrtrans(icorrtrans),correxp
          end if
          yshift=yshift+yshiftinc
          write (ips,1001) xm0,ym0+ymm+yshift,' m'
c         write (ips,1017) trajnam(1:ltrajnam),ifirst,ilast,incr_traj
          call write_traj_lim(ips,' ',1,1,incr_tr,1)
          yshift=yshift+yshiftinc
          write (ips,1001) xm0,ym0+ymm+yshift,' m'
          call psshow(ips,'Title: ',7)
          call psshow(ips,itit,ntit)
          write (ips,*) 'sk'
          write (ips,3000) 12
          call rgbcolor(ips,9)
          write (ips,*) 'np'
c         write (ips,1001) xm0+0.45*xmm,ym0-0.09*ymm09,' m'
          write (ips,1001) xm0+0.45*xmm,ym0-2.0*yshiftinc,' m'
          call psshow(ips,xlab,nxlab)
          write (ips,*) 'sk'
          if (ibondlab .eq. 0 .or. ntotplot .gt. lablim) then
c           Don't skip y label to leave room for res-res info
            write (ips,1001) xm0-0.06*xmm,ym0+0.4*ymm,' translate'
            write (ips,*) '90 rotate'
            write (ips,*) 'np'
            write (ips,1001) 0.0,0.0,' m'
            call psshow(ips,ylab,nylab)
            write (ips,*) 'sk'
            write (ips,*) '-90 rotate'
            write (ips,1001) -(xm0-0.06*xmm),-(ym0+0.4*ymm),' translate'
          else if (ntotplot .gt. lablim) then
           if (ibondlab .eq. 2) write (6,1022) lablim,'residue pairs'
           if (ibondlab .eq. 1) write (6,1022) lablim,
     -       bondname(1:lbondname)
          end if
          write (ips,*) 3,' setlinewidth'
          write (ips,*) 'np'
          write (ips,*) 1,' setlinewidth'
          call rgbcolor(ips,9)
          ifloatx=0
          do ix=1,nxd
            write (ips,*) 'np'
            write (ips,1001) xm0+xmm*float(ix)/float(nxd),ym0,' m'
            write (ips,1001) 0.0,+0.01*ymm,' rlineto'
            write (ips,*) 'sk'
            write (ips,*) 'np'
            write (ips,1001) xm0-0.02*xmm+
     -        xmm*float(4*ix-1)/float(4*nxd),ym0-yshiftinc,' m'
            if (xlab(1:1) .eq. ' ') then
              write (ips,1002) (ix*ixd)
            else
c             if (xdiv .ge. 100.0) then
              if (xdiv .ge. 1.0) then
                ixxdiv=ix*xdiv+0.01
                write (ips,1002) ixxdiv
              else
                write (ips,1003) (ix*xdiv)
              end if
            end if
          end do
          if (ibondlab .eq. 0 .or. ntotplot .gt. lablim) then
            do iy=1,nyd
              write (ips,*) 'np'
              write (ips,1001) xm0,ym0+float(iy)/float(nyd)*ymm,' m'
              write (ips,1001) +0.01*xmm,0.0,' rlineto'
              write (ips,*) 'sk'
              write (ips,*) 'np'
              write (ips,1001) xm0-0.07*xmm,ym0-0.01*ymm+
     -          float(iy)/float(nyd)*ymm,' m'
              write (ips,1002) iy*iyd
            end do
          else
c           Print info
            write (ips,1021) 8
            do iy=1,ntotplot
              if (mod(iy,iyd) .eq. 0) then
                write (ips,*) 'np'
                write (ips,1001) xm0,ym0+float(iy)/float(nyd*iyd)*ymm,
     -            ' m'
                write (ips,1001) +0.01*xmm,0.0,' rlineto'
                write (ips,*) 'sk'
              end if
              write (ips,*) 'np'
              write (ips,1001) xm00,
     -          ym0+float(iy)/float(nyd*iyd)*ymm,' m'
              call makebondlab(iy,iy,iy-1,ibondlab,bondlab,lbondlab,
     -          irrix,ixresno,ixres,index,resnames,atnames,1,mxbonds,
     -          maxrsd,maxrec)
              write (ips,1013) bondlab(1)(1:lbondlab)
            end do
            write (ips,1021) 12
          end if
          if (iclust .gt. 0) then
c           Plot cluster marker bars
            write (ips,*) 12,' setlinewidth'
            xclust0=xm0+1.02*xmm
            ncp=0
            do ic=1,nclust
              if (ilhbclust(ic)-ifhbclust(ic) .gt. 0) then
                ncp=ncp+1
                xinc=1.00
c               if (mod(ncp,2) .eq. 0) xinc=1.01
                write (ips,*) 'np'
                write (ips,1001) xinc*xclust0,ym0+yfac*ifhbclust(ic),'m'
                write (ips,1001) xinc*xclust0,ym0+yfac*ilhbclust(ic),'l'
                write (ips,*) 'sk'
              end if
            end do
          end if
c         write (ips,*) 2,' setlinewidth'
          write (ips,*) lw,' setlinewidth'
        end if
c       write (ips,*) lw,' setlinewidth'
        nlwp=0
        do ib=1,ntotplot
          ibx=index(ib)
          if (nobullet .eq. 0) then
            if (ihbstart(ibx) .eq. 0 .and. ihb(ibx) .eq. 1) then
              ihbstart(ibx)=1
              write (ips,*) 'np'
              call rgbcolor(ips,9)
              write (ips,1016) xm0+xfac*float(ifr),ym0+yfac*float(ib),
     -          float(lw),0.0,360.0
              write (ips,*) 'fill'
              write (ips,*) 'sk'
            end if
          end if
          if (ihb(ibx) .eq. 1 .and. ihbon(ibx) .eq. 0) then
            ihbon(ibx)=ifr
          else if ((ihb(ibx) .eq. 0  .or. ifr .eq. nframetot) .and.
     -              ihbon(ibx) .gt. 0) then
c           Draw line between ihbon(ibx) and ifr
            ym=ym0+yfac*float(ib)
            write (ips,*) 'np'
            call rgbcolor(ips,-icolrb(ianc_anc(ibx)+1))
            write (ips,1001) xm0+xfac*float(ihbon(ibx)),ym,' m'
            write (ips,1001) xm0+xfac*float(ifr),ym,' l'
            write (ips,*) 'sk'
            ihbon(ibx)=0
          end if
        end do
      end do
c     Print tick on the right y axis
      call rgbcolor(ips,9)
      write (ips,*) 1,' setlinewidth'
      do iy=1,nyd
        write (ips,*) 'np'
        write (ips,1001) xm0+xmm,ym0+float(iy)/float(nyd)*ymm,' m'
        write (ips,1001) -0.01*xmm,0.0,' rlineto'
        write (ips,*) 'sk'
      end do
c     Print color code
      write (ips,*) lw,' setlinewidth'
      axshift=0.0
      if (nbondmax .gt. naabondmax) then
        write (ips,*) 'np'
        call rgbcolor(ips,-icolrb(1))
c       write (ips,1001) xm0+10.0,ym0-0.10*ymm09,' m'
        write (ips,1001) xm0+10.0,ym0-0.08*ymm09,' m'
        write (ips,1001) 20.0,0.0,' rlineto'
        write (ips,1001) xm0+30.0,ym0-0.08*ymm09,' m'
        write (ips,1011) ' : anchor -- non-anchor'
        write (ips,*) 'sk'
        axshift=1.0
      end if
      if (naabondmax .gt. 0) then
        write (ips,*) lw,' setlinewidth'
        write (ips,*) 'np'
        call rgbcolor(ips,-icolrb(2))
        write (ips,1001) xm0+10.0+axshift*0.6*xmm,ym0-0.08*ymm09,' m'
        write (ips,1001) 20.0,0.0,' rlineto'
        write (ips,1001) xm0+30.0+axshift*0.6*xmm,ym0-0.08*ymm09,' m'
        write (ips,1011) ' : anchor -- anchor'
        write (ips,*) 'sk'
      end if
      call rgbcolor(ips,9)
      write (ips,*) 1,' setlinewidth'
      write (ips,*) 'showpage'
      if (ibondlab .eq. 2) call trnsfi(index,itemp1,nhbtot)
      return
1001  format(2f8.1,1x,a)
1002  format('(',i8,') show',/,'sk')
1003  format('(',f8.2,') show',/,'sk')
1004  format('( Minimum number of hydrogens on hydrophobic carbons=',
     -  i1,') show')
1005  format('( Length limit of a ',a,' bond=',f5.2,' A) show')
1007  format('(Minimum and maximum % occurrences to show=',
     -  f4.1,f6.1,') show')
1008  format('(Minimum and maximum interresidue distances to show=',
     -  i4,' and ',i5,' residues ) show')
1009  format('(Unclustered plot, original ',a,'-bond order ) show')
1010  format('(Clustering method: ',a,'; distance measure cutoff=',f5.2,
     -  ' ) show')
1011  format('(',a,') show')
1012  format('(Clustering method: ',a,' Number of clusters=',i2,')show')
1013  format('(',a,' ) show')
1014  format('(# of frames averaged=',i4,' ) show')
1015  format('(Hydrogen-bond parameters: hblimfac=',f6.3,
     -  ' angmin=',f8.2,' deg ) show')
1016  format(5f12.4,' arc')
c1017 format('(Trajectory analyzed:',a,' Frames',i7,' - ',i6,
c    -  ' Increment:',i5,') show')
1018  format('(Residue aggregated time course of ',a,' bonds) show')
1019  format('%%Page: 1 ',i4)

1021  format('/Helvetica findfont',/,i2,' scalefont',/,'setfont')
1022  format(' NOTE: more than',i3,1x,a,' - no room for the names')
1023  format('(; distance measure=(1',a,'corr)^',f4.2,' ) show')
3000  format('/m { moveto } def',/,'/l { lineto } def',/,
     -  '/np { newpath } def',/, '/sk { stroke } def',/,
     -  '/f { fill } def',/,'/lw { setlinewidth } def',/,
     -  '/Helvetica findfont',/,i2,' scalefont',/,'setfont')
      end
