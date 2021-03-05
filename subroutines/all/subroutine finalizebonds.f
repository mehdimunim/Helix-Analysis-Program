      subroutine finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
     -  numres,nresslt,npspages,ipspage,nres2d,ibondcorr,iresbondcorr,
     -  nhneigmin,hblimfac,angmin,rhph_sltbmax,inamcol1,inamcol2,
     -  irescol1,irescol2,ifhb2d,ilhb2d,nhbdist,rhbdist,iframeunit,
     -  framefac,title,ltitle,xtrajlab,lxtrajlab,xtraj,value,ifa_s,
     -  ila_s,ih,cv,temp,index,iresno,ifres,isegno,atnames,resnames,
     -  ixres,ixresno,ixsegno, indexa,indexs,index2d,ianc_anc,isc,
     -  ixclst,irepav,irepmx,irepeng,irepkm,rmsdlim,engcl,it1,it2,it3,
     -  it4,it5,irrix,itemp1,itemp2,itemp3,nrrbond,line,bondtype,
     -  lbondtype,ibondtype,label2d,iselfanc,ianchor2,iresshift,
     -  ifailbond,nbondavg,inpfile,namleni,maxbondf,nmcmaxbond,ncolcode,
     -  maxcolcode,maxbondcount,mxbonds,maxrsd,maxfrm,maxrec,mx2d)
      dimension index(maxrec),iresno(maxrec),ifres(maxrec),
     -  isegno(maxrec),ixres(maxrec),ixresno(maxrsd),ixsegno(maxrsd),
     -  indexs(maxrec),index2d(mxbonds),ianc_anc(mxbonds),engcl(mx2d),
     -  it1(mxbonds),it2(mxbonds),it3(mxbonds),it4(mxbonds),
     -  it5(mxbonds),irrix(maxrec),itemp1(maxrec),itemp2(maxrec),
     -  itemp3(maxrec),nrrbond(mxbonds),ifhb2d(mxbonds),ilhb2d(mxbonds),
     -  irepav(mx2d),irepmx(mx2d),irepeng(mx2d),irepkm(mx2d),
     -  nhbdist(mxbonds),rhbdist(mxbonds),value(mxbonds),ifa_s(mxbonds),
     -  ila_s(mxbonds),ih(maxrec),cv(maxrec),temp(maxrec),isc(maxrec),
     -  ixclst(mxbonds),rmsdlim(mxbonds),indexa(maxrec),xtraj(maxfrm)
      character*8 atnames(maxrec),resnames(maxrsd)
      character*80 title,label2d(mx2d),label,remark
      character*(*) xtrajlab,inpfile
      character* 132 line(maxrec)
      character*(*) bondtype
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      parameter (MAXCOPY1=MAXCOPY-1)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY1),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      parameter (MAXBONDS=10000)
      dimension nhbpers(MAXBONDS),itrackf(MAXBONDS),itrackl(MAXBONDS),
     -  maxlenon(MAXBONDS),maxlenoff(MAXBONDS),itrack(MAXFRAMES),
     -  xtrajmod(MAXFRAMES),data(MAXFRAMES)
      dimension c(3,1),ct1(3,1),ct2(3,1)
      real*8 dsum,dsum2,dsumm,dsumm2,bsum,bsum2,aabsum,aabsum2
      character*37 yclab(MAXBONDS)
      character*80 plotlab
      character*200 aucfile
c     print *,'FINALIZEBOND nbondavg=',nbondavg
      if (ifailbond .ne. 0) return
      if (nbfound .eq. 0) then
        write (6,2001)
        write (iw0,2001)
        return
      end if
c     Turn on sorting in clique clustering
      nosortindex2d=0
      lablim=50
      icorrtyp=0
      icorrtrans=0
      call zeroiti(nhbpers,0,MAXBONDS)
      call zeroiti(itrackf,0,MAXBONDS)
      call zeroiti(itrackl,0,MAXBONDS)
      plotlab(1:lbondtype)=bondtype(1:lbondtype)
      plotlab(lbondtype+1:lbondtype+7)=' bond #'
      lplotlab=lbondtype+7
c     print *,'FINALIZEBOND nbfound,nframe,ianchor2,nbresfound=',
c    -  nbfound,nframe,ianchor2,nbresfound
c     print *,'Number of ',bondtype(1:lbondtype),' bonds found=',nbfound
      if (nbfound .eq. 0) return
      icorrbothstart=0
      call askyn('Do you want to end the track statistics at '//
     -  'the last on frame',60,1,-1,iuselaston,135,0)
      if (ibondcorr+iresbondcorr .gt. 0)
     -  call askyn('Do you want to correlate only after both '//
     -    'tracks started',55,1,-1,icorrbothstart,136,0)
      call trackstat(nbfound,nhbdist,nhbpers,maxlenon,maxlenoff,
     -  itrackf,itrackl,it1,it2)
8005  call indexit(index2d,1,nbfound,0)
      write (iw0,2007) bondtype(1:lbondtype),maxbondf,nmcmaxbond
      write (6,2007) bondtype(1:lbondtype),maxbondf,nmcmaxbond
      call printbonddist(nframe,nbfound,nhbdist,index2d,
     -  ianc_anc,naabond,iw0,bondtype,lbondtype,iselfanc,maxbondcount,
     -  mxbonds)
      nbfoundorig=nbf ound
      call filterbonds(n,nbfound,nhbdist,rhbdist,nhbpers,itrackf,
     -  itrackl,itrack,maxlenon,maxlenoff,iframeunit,framefac*increment,
     -  nframe,index2d,bondtype,lbondtype,line,index,inamcol1,inamcol2,
     -  irescol1,irescol2,iresno,isegno,ixres,ixresno,ixsegno,ifres,
     -  numres,nresslt,ianc_anc,isc,bondtype,lbondtype,percminb,
     -  percmaxb,minresdistb,maxresdistb,nochange,cv,iuselaston,iauc,
     -  iaucw,inpfile,namleni,temp,it1,it2,it3,it4,irrix,itemp2,iw0,
     -  maxrec,maxrsd,mxbonds,MAXFRAMES,MAXCOPY)
c     cv contains the % of frames an atom is forming a bond
      call bondsum(nbfoundorig,nbfound,index2d,ixres,atnames,resnames,
     -  isegno,iresno,ifres,it1,iw0,MAXBONDS,maxrsd,maxrec)
      call makebondlab(1,nbfound,0,1,yclab,lyclab,irrix,ixresno,
     -  ixres,index2d,resnames,atnames,nbfound,mxbonds,maxrsd,maxrec)
      if (iaucw .gt. 0) then
        call getname(aucfile,laucfile,
     -    'Name of the file to write the bond autocorrelations',51,80,
     -    '',0,0,0,0)
        call print_auc(aucfile,laucfile,nbfound,nframe,yclab,lyclab)
        write (6,2002) aucfile(1:laucfile)
        write (iw0,2002) aucfile(1:laucfile)
      end if
      rhbcorrclust=0.0
      nhbcorrclust=1
      print *,'NFRAME,XTRAJ(1),XTRAJ(NFRAME)=',
     -  nframe,xtraj(1),xtraj(nframe)
      call plotbond(iw1,nbfoundorig,nbfound,nbresfound,nhbcorrclust,
     -  ianc_anc,ifhb2d,ilhb2d,index2d,indexs,it1,it2,rhbcorrclust,0,
     -  xtraj(nframe),10,title,ltitle,xtrajlab,lxtrajlab,plotlab,
     -  lplotlab,0,percminb,percmaxb,numres,minresdistb,maxresdistb,
     -  rhph_sltbmax,bondtype,lbondtype,0,nhneigmin,nbondmax,naabondmax,
     -  icorrtyp,icorrtrans,nframeav,correxp,hblimfac,angmin,lablim,
     -  ibondtype,1,it3,ixresno,itemp3,ixres,nrrbond,atnames,
     -  resnames,npspages,ipspage,mxbonds,maxrsd,maxrec)
c     Plot the total # of bonds for successive frames
      call roundlimint(naabondmax+1,iyd2,nyd2)
      call roundlimint(nbondmax+1,iyd,nyd)
      nplot=2
      iyinc=0
      if (iselfanc .eq. 0 .or. ianchor2 .eq. 1
     -   .or. naabond .eq. nbfoundorig) then
        nplot=1
        if (iselfanc .eq. 1) iyinc=1
      end if
      dsum=0.d0
      dsumm=0.d0
      dsum2=0.d0
      dsumm2=0.d0
      nbmin=10000
      nbmax=0
      naabmin=10000
      naabmax=0
      do i=1,nframe
        dsum=dsum+scres(1,i)
        dsumm=dsumm+scres(2,i)
        dsum2=dsum2+scres(1,i)**2
        dsumm2=dsumm2+scres(2,i)**2
        if (nbmin .gt. scres(1,i)) nbmin=scres(1,i)
        if (nbmax .lt. scres(1,i)) nbmax=scres(1,i)
        if (naabmin .gt. scres(2,i)) naabmin=scres(2,i)
        if (naabmax .lt. scres(2,i)) naabmax=scres(2,i)
c       write (iw0,*) i,' scr1=',scres(1,i),' dsum=',dsum
      end do
      write (iw0,2003) ' ',bondtype(1:lbondtype),dsum/dfloat(nframe),
     -  dsqrt(dsum2/dfloat(nframe)-(dsum/dfloat(nframe))**2),nbmin,nbmax
      if (iselfanc .eq. 1) write (iw0,2003) ' anchor-anchor ',
     -  bondtype(1:lbondtype),dsumm/dfloat(nframe),
     -  dsqrt(dsumm2/dfloat(nframe)-(dsumm/dfloat(nframe))**2),
     -  naabmin,naabmax
      bsum=0.d0
      aabsum=0.d0
      bsum2=0.d0
      aabsum2=0.d0
      do ifr=1,nframe
        bsum=bsum+scres(1,ifr)
        bsum2=bsum2+scres(1,ifr)**2
        aabsum=aabsum+scres(2,ifr)
        aabsum2=aabsum2+scres(1,ifr)**2
        data(ifr)=scres(1,ifr)
      end do
      avgbond=bsum/dfloat(nframe)
      avgaabond=aabsum/dfloat(nframe)
      bsd=dsqrt(bsum2/dfloat(nframe)-(bsum/dfloat(nframe))**2)
      aabsd=dsqrt(aabsum2/dfloat(nframe)-(aabsum/dfloat(nframe))**2)
      call batchmean(nframe,0,data,'Average # of bonds',18,iw0,0,av,sd,
     -  ci)
      if (nbondavg .gt. 1) then
        nfr=nframe/nbondavg
        incr=0
        do iaa=1,nfr
          nbsum=0
          nbaasum=0
          do ia=incr+1,incr+nbondavg
            nbsum=nbsum+scres(1,ia)
            nbaasum=nbaasum+scres(2,ia)
          end do
          scres(1,iaa)=float(nbsum)/float(nbondavg)
          scres(2,iaa)=float(nbaasum)/float(nbondavg)
          xtrajmod(iaa)=xtraj(incr+nbondavg)
          incr=incr+nbondavg
        end do
        write (remark,2004) nbondavg,avgbond,bsd
        lremark=65
      else
        nfr=nframe
        call trnsfr(xtrajmod,xtraj,nframe)
        write (remark,2005) avgbond,bsd
        lremark=22
      end if
      if (ci .ne. 999.0) then
        write (remark(lremark+1:lremark+15),2008) ci
        lremark=lremark+15
      end if
      if (nplot .eq. 2) then
        write (remark(lremark+1:lremark+14),2006) avgaabond,aabsd
        lremark=lremark+25
      end if
      llabel=16+lbondtype
      label(1:llabel)='Number of '//bondtype(1:lbondtype)//' bonds'
      call plot2fun(iw1,nplot,xtrajmod,scres,scres,nfr,0.0,0.0,0,0.0,
     -  float(iyd),nyd,0.0,float(iyd2),nyd2,title,80,remark,lremark,
     -  xtrajlab,lxtrajlab,label,llabel,'Number of anchor-anchor bonds',
     -  29,1,0,6,2,1,0,0,iyinc,0,ipspage,1,1,0)
      if (ibondcorr .gt. 0 .and. nbfoundorig .gt. mx2d) then
        write (6,2000) nbfoundorig,mx2d
        ibondcorr=0
      end if
      if (ibondcorr .gt. 0) then
c       Calculate, plot and print the dot product of bonds
        ireadwrite=1
        ifhb2d(1)=1
        ilhb2d(1)=nbfound
        call bondcorrsum(nbfoundorig,0,scpmin,scpmax,it1,it2,it3,
     -    correxp,icorrtyp,icorrtrans,0,icorrbothstart,nframeav,0,
     -    iw0,mx2d)
        call bondcorrprint(nbfound,iw0,line,index,iresno,index2d,
     -    correxp,icorrtyp,icorrtrans,nframeav,nhbdist,ifhb2d,ilhb2d,
     -    ixclst,rhbcorrclust,nhbcorrclust,0,inamcol1,inamcol2,irescol1,
     -    irescol2,nhneigmin,hblimfac,angmin,rhph_sltbmax,percminb,
     -    percmaxb,minresdistb,maxresdistb,numres,1,nframe,bondtype,
     -    lbondtype,0,ixresno,resnames,maxrec,mxbonds)
        call plotbondcorr(iw1,nbfound,yclab,lyclab,0,index2d,
     -    title,temp,itemp2,lablim,ncolcode,maxcolcode,ipspage,0)
        rhbcorrclustdef=(scpmax-scpmin)/2.0
        call clusterdistr(nbfound,iw0,rmsdlim,scpmin,scpmax,nhbdist,it1,
     -    it2,it3,it4,ifhb2d,ilhb2d,nhbcorrclust,indexa,index2d,
     -    ixclst,it5,value,ifa_s,ila_s,ih,temp,rhbcorrclustdef,
     -    rhbcorrclust,res(1,1,11),0,'Correlation-based distance',26,0,
     -    0,irepav,irepmx,irepeng,irepkm,engcl,c,ct1,ct2,1,27,
     -    iclstyp,nrrbond,1,label2d,80,0,1,nosortindex2d,1,mxbonds,
     -    maxframe)
c       print *,'Plotting the clustered ',bondtype(1:lbondtype),' bonds'
        call plotbond(iw1,nbfoundorig,nbfound,nbresfound,nhbcorrclust,
     -    ianc_anc,ifhb2d,ilhb2d,index2d,indexs,it1,it2,rhbcorrclust,
     -    iclstyp,xtraj(nframe),10,title,ltitle,xtrajlab,lxtrajlab,
     -    'Bond #',6,1,percminb,percmaxb,numres,minresdistb,maxresdistb,
     -    rhph_sltbmax,bondtype,lbondtype,1,nhneigmin,nbondmax,
     -    naabondmax,icorrtyp,icorrtrans,nframeav,correxp,hblimfac,
     -    angmin,lablim,ibondtype,1,it3,ixresno,itemp3,ixres,
     -    nrrbond,atnames,resnames,npspages,ipspage,mxbonds,maxrsd,
     -    maxrec)
c       print *,'Printing the clustered ',bondtype(1:lbondtype),' bonds'
        call bondcorrprint(nbfound,iw0,line,index,iresno,index2d,
     -    correxp,icorrtyp,icorrtrans,nframeav,nhbdist,ifhb2d,ilhb2d,
     -    ixclst,rhbcorrclust,nhbcorrclust,1,inamcol1,inamcol2,irescol1,
     -    irescol2,nhneigmin,hblimfac,angmin,rhph_sltbmax,percminb,
     -    percmaxb,minresdistb,maxresdistb,numres,0,nframe,bondtype,
     -    lbondtype,0,ixresno,resnames,maxrec,mxbonds)
        call askyn(
     -    'Do yo want to repeat clustering with different filters',
     -     54,1,-1,irepscan,0,0)
        if (irepscan .eq. 1) then
          nbfound=nbfoundorig
          go to 8005
        end if
      end if
      call mapbondstorespairs(nbfound,nbfoundorig,nbresfilt,nbresfound,
     -  nframe,ixres,it1,bondtype,lbondtype,line,index,irescol1,
     -  irescol2,nresslt,percminb,percmaxb,minresdistb,maxresdistb,
     -  percminr,percmaxr,minresdistr,maxresdistr,ifres,
     -  iresno,isegno,nochange,index2d,ifa_s,ila_s,itemp1,itemp2,itemp3,
     -  iw0,mxbonds,maxrsd,maxrec)
c     if (iauc+iresbondcorr .gt. 0)
c    -  call condensetracks(nbfoundorig,nbresfilt,it1,it2,iw0,mxbonds)
      call condensetracks(nbfoundorig,nbresfilt,it1,it2,iw0,mxbonds)
      call res_res_bond(nres2d,nbfound,index2d,nhbdist,iresno,ixres,
     -  ifres,isegno,irrix,itemp2,itemp3,title,ltitle,bondtype,
     -  lbondtype,line,index,irescol1,irescol2,iresshift,iw0,iw1,
     -  ipspage,ncolcode,maxcolcode,mxbonds,maxrsd,maxrec)
      call plotbond(iw1,nbfoundorig,nbfound,nbresfilt,nhbcorrclust,
     -  ianc_anc,ifhb2d,ilhb2d,index2d,indexs,it1,it2,rhbcorrclust,0,
     -  xtraj(nframe),10,title,ltitle,xtrajlab,lxtrajlab,
     -  'Residue pair #',14,1,percminr,percmaxr,numres,minresdistr,
     -  maxresdistr,rhph_sltbmax,bondtype,lbondtype,0,nhneigmin,
     -  nbondmax,naabondmax,icorrtyp,icorrtrans,nframeav,correxp,
     -  hblimfac,angmin,lablim,ibondtype,2,it3,ixresno,irrix,ixres,
     -  nrrbond,atnames,resnames,npspages,ipspage,mxbonds,maxrsd,maxrec)
      if (iresbondcorr .gt. 0) then
        call bondcorrsum(nbfoundorig,nbresfilt,scpmin,scpmax,it1,it2,
     -    it3,correxp,icorrtyp,icorrtrans,1,icorrbothstart,
     -    nframeav,1,iw0,mx2d)
        call indexit(itemp1,1,nbresfilt,0)
        call makebondlab(1,nbresfilt,0,2,yclab,lyclab,irrix,ixresno,
     -    ixres,itemp1,resnames,atnames,nbresfilt,mxbonds,maxrsd,maxrec)
        call plotbondcorr(iw1,nbresfilt,yclab,lyclab,1,itemp1,title,
     -    temp,itemp2,lablim,ncolcode,maxcolcode,ipspage,0)
        rhbcorrclustdef=(scpmax-scpmin)/2.0
        ifhb2d(1)=1
        ilhb2d(1)=nbresfound
        call clusterdistr(nbresfilt,iw0,rmsdlim,scpmin,scpmax,nhbdist,
     -    it1,it2,it3,it4,ifhb2d,ilhb2d,nhbcorrclust,indexa,irrix,
     -    ixclst,it5,value,ifa_s,ila_s,ih,temp,rhbcorrclustdef,
     -    rhbcorrclust,res(1,1,11),0,'Correlation-based distance',26,0,
     -    0,irepav,irepmx,irepeng,irepkm,engcl,c,ct1,ct2,1,27,
     -    iclstyp,nrrbond,1,label2d,80,0,1,nosortindex2d,0,mxbonds,
     -    maxframe)
        call bondcorrprint(nbresfilt,iw0,line,index,iresno,irrix,
     -    correxp,icorrtyp,icorrtrans,nframeav,nhbdist,ifhb2d,ilhb2d,
     -    ixclst,rhbcorrclust,nhbcorrclust,1,inamcol1,inamcol2,irescol1,
     -    irescol2,nhneigmin,hblimfac,angmin,rhph_sltbmax,percminr,
     -    percmaxr,minresdistr,maxresdistr,numres,1,nframe,bondtype,
     -    lbondtype,1,ixresno,resnames,maxrec,mxbonds)
        call plotbond(iw1,nbfoundorig,nbfound,nbresfilt,nhbcorrclust,
     -    ianc_anc,ifhb2d,ilhb2d,ixclst,indexs,it1,it2,rhbcorrclust,
     -    iclstyp,xtraj(nframe),10,title,ltitle,xtrajlab,lxtrajlab,
     -    'Residue pair #',14,1,percminr,percmaxr,numres,minresdistr,
     -    maxresdistr,rhph_sltbmax,bondtype,lbondtype,0,nhneigmin,
     -    nbondmax,naabondmax,icorrtyp,icorrtrans,nframeav,correxp,
     -    hblimfac,angmin,lablim,ibondtype,2,it3,ixresno,irrix,
     -    ixres,nrrbond,atnames,resnames,npspages,ipspage,mxbonds,
     -    maxrsd,maxrec)
      end if
      call askyn(
     -  'Do you want to calculate residue-residue autocorrelation',56,
     -   1,-1,iauc,0,0)
      if (iauc .eq. 1)
     -  call autocorr_res_res(nbresfilt,nframe,iframeunit,
     -    framefac*increment,itrack,nhbdist,nhbpers,maxlenon,
     -    maxlenoff,itrackf,itrackl,ifres,iresno,isegno,resnames,
     -    irescol2-irescol1+1,it2,it3,iuselaston,iw0,mxbonds,MAXFRAMES,
     -    maxrsd,maxrec)
      return
2000  format(' ERROR: the number of bonds found (',i5,') exceeds the',
     -  ' limit (',i5,')',/,' either increase the parameter MAX2D and ',
     -  ' recompile Simulaid',/,' or select fewer atoms to analyze')
2001  format(' NOTE: no bonds were found')
2002  format(' Full autocorrelation functions were written to file',/,
     -  1x,a)
2003  format(' Average number of',a,a,' bonds/frame=',f6.1,' S.D.=',
     -  f5.1,/,' Range: [',i4,',',i4,']')
2004  format(' Number of bonds averaged over',i5,' frames; <NB>=',f5.1,
     -  ' S.D.=',f5.1)
2005  format(' <NB>=',f5.1,' S.D.=',f5.1)
2006  format('; <NBaa>=',f5.1,' S.D.=',f5.1)
2007  format(' Maximum number of ',a,' bonds=',i5,' at frame# ',i8)
2008  format(' CI(1sig)=',f5.1)
      end
