      subroutine filterbonds(n,nbfound,nhbdist,rhbdist,nhbpers,itrackf,
     -  itrackl,itrack,maxlenon,maxlenoff,iframeunit,framefac,nframe,
     -  indexbond,label,llabel,line,index,inamcol1,inamcol2,irescol1,
     -  irescol2,iresno,isegno,ixres,ixresno,ixsegno,ifres,numres,
     -  nresslt,ianc_anc,isc,bondtype,lbondtype,percmin,percmax,
     -  minresdist,maxresdist,nochange,percbond,iuselaston,iauc,iaucw,
     -  inpfile,namleni,temp,it1,it2,it3,it4,irrix,itemp1,iout,maxrec,
     -  maxrsd,mxbonds,mxframes,mxcopy)
      dimension index(maxrec),nhbdist(mxbonds),rhbdist(mxbonds),
     -  nhbpers(mxbonds),itrackf(mxbonds),itrackl(mxbonds),
     -  itrack(mxframes),maxlenon(mxbonds),maxlenoff(mxbonds),
     -  indexbond(mxbonds),iresno(maxrec),isegno(maxrec),ixres(maxrec),
     -  ixresno(maxrsd),ixsegno(maxrsd),ifres(maxrec),ianc_anc(mxbonds),
     -  isc(maxrec),percbond(maxrec),temp(maxrec),it1(mxbonds),
     -  it2(mxbonds),it3(mxbonds),it4(mxbonds),irrix(maxrec),
     -  itemp1(maxrec)
      character*(*) label,bondtype,inpfile
      character*80 bond,listfile
      character*132 line(maxrec)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      dimension limperc(6),nplim(6),npbsum(6)
      data limperc /1,2,5,10,20,50/
      character*6 bondlab
c     print *,'FILTERBOND NBFOUND=',nbfound
      ioutl=0
      nochange=1
      call getfiltlims(0.0,100.0,1,nresslt-1,percmin,percmax,
     -  minresdist,maxresdist,nresslt,nochange,label,llabel,' ',1,
     -  iout)
      do ip=1,6
        nplim(ip)=max0(1,nframe/(100/limperc(ip)))
        npbsum(ip)=0
      end do
      nhbsum=0
      nbb=0
      naa=0
      call zeroit(percbond,n)
      do ii=1,nbfound
        i=indexbond(ii)
        temp(i)=100.0*float(nhbdist(i))/float(nframe)
        nhbsum=nhbsum+nhbdist(i)
        do ip=1,6
          if (nhbdist(i) .gt. nplim(ip)) npbsum(ip)=npbsum(ip)+1
        end do
        if (ianc_anc(i) .eq. 1) naa=naa+1
        if (isc(ihbpair(1,i))*isc(ihbpair(2,i)) .eq. 0) nbb=nbb+1
      end do
      do ip=1,6
        write (iout,1007) limperc(ip),npbsum(ip)
      end do
      if (naa .gt. 0) write (iout,1008)
      if (nbb .gt. 0) write (iout,1009)
      nhbsumtot=nhbsum
      iallkept=1
      if (nochange .eq. 0) then
        ndel=0
        nmindel=0
        nmaxdel=0
        nmindistdel=0
        nmaxdistdel=0
        do ii=1,nbfound
          i=indexbond(ii)
          p=temp(i)
          ir1=ihbpair(1,i)
          ir2=ihbpair(2,i)
          iresdist=iabs(ixres(ihbpair(1,i))-ixres(ihbpair(2,i)))
          if (iresdist .lt. minresdist .or. iresdist .gt. maxresdist
     -      .or. p .lt. percmin .or. p .gt. percmax) then
            ndel=ndel+1
            if (p .lt. percmin) then
              nmindel=nmindel+1
              it1(nmindel)=i
            else if (p .gt. percmax) then
              nmaxdel=nmaxdel+1
              it2(nmaxdel)=i
            end if
            if (iresdist .lt. minresdist) then
              nmindistdel=nmindistdel+1
              it3(nmindistdel)=i
            else if (iresdist .gt. maxresdist) then
              nmaxdistdel=nmaxdistdel+1
              it4(nmaxdistdel)=i
            end if
          else
            indexbond(ii-ndel)=indexbond(ii)
          end if
        end do
        iallkept=0
        if (ndel .gt. 0) then
          nbfound=nbfound-ndel
          if (nmindel .gt. 0) then
            write (iout,1002) nmindel,'infrequent'
            do ib=1,nmindel
              call bonddescr(it1(ib),ihbpair,line,index,iresno,isegno,
     -          inamcol1,inamcol2,irescol1,irescol2,bond,lbond,bondlab,
     -          lbondlab,nbb,ianc_anc,isc,ia1,ia2,ir1,ir2,maxrec,
     -          MAXBONDS)
              write (iout,1005) 'BD',ib,it1(ib),bond(1:lbond),
     -          temp(it1(ib)),bondlab(1:lbondlab),
     -          rhbdist(it1(ib))/nhbdist(it1(ib))
              nhbsum=nhbsum-nhbdist(it1(ib))
            end do
          end if
          if (nmaxdel .gt. 0) then
            write (iout,1002) nmaxdel,'persistent'
            do ib=1,nmaxdel
              call bonddescr(it2(ib),ihbpair,line,index,iresno,isegno,
     -          inamcol1,inamcol2,irescol1,irescol2,bond,lbond,bondlab,
     -          lbondlab,nbb,ianc_anc,isc,ia1,ia2,ir1,ir2,maxrec,
     -          MAXBONDS)
              write (iout,1005) 'BD',ib,it2(ib),bond(1:lbond),
     -          temp(it2(ib)),bondlab(1:lbondlab),
     -          rhbdist(it2(ib))/nhbdist(it2(ib))
              nhbsum=nhbsum-nhbdist(it2(ib))
            end do
          end if
          if (nmindistdel .gt. 0) then
            write (iout,1002) nmindistdel,'close'
            do ib=1,nmindistdel
              call bonddescr(it3(ib),ihbpair,line,index,iresno,isegno,
     -          inamcol1,inamcol2,irescol1,irescol2,bond,lbond,bondlab,
     -          lbondlab,nbb,ianc_anc,isc,ia1,ia2,ir1,ir2,maxrec,
     -          MAXBONDS)
              write (iout,1005) 'BD',ib,it3(ib),bond(1:lbond),
     -          temp(it3(ib)),bondlab(1:lbondlab),
     -          rhbdist(it3(ib))/nhbdist(it3(ib))
              nhbsum=nhbsum-nhbdist(it3(ib))
            end do
          end if
          if (nmaxdistdel .gt. 0) then
            write (iout,1002) nmaxdistdel,'distant'
            do ib=1,nmaxdistdel
              call bonddescr(it4(ib),ihbpair,line,index,iresno,isegno,
     -          inamcol1,inamcol2,irescol1,irescol2,bond,lbond,bondlab,
     -          lbondlab,nbb,ianc_anc,isc,ia1,ia2,ir1,ir2,maxrec,
     -          MAXBONDS)
              write (iout,1005) 'BD',ib,it4(ib),bond(1:lbond),
     -          temp(it4(ib)),bondlab(1:lbondlab),
     -          rhbdist(it4(ib))/nhbdist(it4(ib))
              nhbsum=nhbsum-nhbdist(it4(ib))
            end do
          end if
          write (6,1003) label(1:llabel),nbfound
        else
          write (iout,*) 'All bonds are kept'
          iallkept=1
        end if
      end if
      write (iout,1004) nbfound,label(1:llabel)
      call zeroiti(irrix,0,n)
      call zeroiti(itemp1,0,numres)
      iaucw=0
      loffmin=0
      iauctype=0
      nseg_scr=0
      nreusemax=nframe
      call askyn('Do you want to calculate bond autocorrelation',45,
     -  1,-1,iauc,0,0)
      if (iauc .eq. 1)
     -  call auc_params(iauctype,lastframeinp,loffmin,nreusemax,
     -    nseg_scr,iaucw,nframe,iframeunit,framefac,'Bond',4,
     -    iuselaston,iout,mxframes)
      nauc_extra=0
      call askyn(
     -  'Do you want a list of bonded atom pairs',35,1,-1,ilistfile,0,7)
      if (ilistfile .eq. 1) then
        ioutl=50
        call changeext(inpfile,listfile,namleni,llistfile,'pls',3,0,0)

        call openfile(50,0,' ',1,'new',listfile,llistfile,notfnd,0,1,1,
     -    1,0)
        write (6,1012) listfile(1:llistfile)
        write (iout,1012) listfile(1:llistfile)
      end if
      do ib=1,nbfound
        i=indexbond(ib)
        call bonddescr(i,ihbpair,line,index,iresno,isegno,inamcol1,
     -    inamcol2,irescol1,irescol2,bond,lbond,bondlab,lbondlab,nbb,
     -    ianc_anc,isc,ia1,ia2,ir1,ir2,maxrec,MAXBONDS)
        if (ilistfile .eq. 1) write (ioutl,1001) ia1,ia2
        irrix(ia1)=irrix(ia1)+nhbdist(i)
        irrix(ia2)=irrix(ia2)+nhbdist(i)
        itemp1(ir1)=itemp1(ir1)+nhbdist(i)
        itemp1(ir2)=itemp1(ir2)+nhbdist(i)
c        write (06,9781) ib,i,nhbpers(i),nhbdist(i)
c9781    format(' i=',i5,' ib=',i5,'  nhbpers=',i5,' nhbdist=',i5)
        write (iout,1005) 'BF',ib,i,bond(1:lbond),temp(i),
     -    bondlab(1:lbondlab),rhbdist(i)/nhbdist(i)
        if (iauc .eq. 1)
     -    call persistence(nhbdist(i),nhbpers(i),itrackf(i),itrackl(i),
     -      maxlenon(i),maxlenoff(i),iframeunit,framefac,nframe,
     -      iuselaston,iout)
        if (iuselaston .eq. 1) then
          lentrack=itrackl(i)-itrackf(i)+1
        else
          lentrack=nframe-itrackf(i)+1
        end if
        percon=float(nhbdist(i))/float(lentrack)
        if (iauc .gt. 0) then
          call getbondtrack(i,itrack,ifirstframe,lastframe,30,nframe)
          call autocorr(ib,i,itrack,ifirstframe,lastframe,iframeunit,
     -      framefac,iauctype,lastframeinp,nreusemax,percon,loffmin,
     -      nseg_scr,nauc_extra,nframe,iout,mxframes)
        end if
      end do
      if (ilistfile .eq. 1) close (ioutl)
      call checkdim(mxcopy+nauc_extra,mxcopy,'MAXCOPY',7,
     -  'number of AUCs',14,iout)
      do ia=1,n
        percbond(ia)=float(irrix(ia))/float(nframe)
      end do
      write (iout,1006) bondtype(1:lbondtype),'=',
     -  float(nhbsumtot)/float(nframe)
      if (iallkept .eq. 0) write (iout,1006) bondtype(1:lbondtype),
     -  ' without the deleted ones=',float(nhbsum)/float(nframe)
      write (iout,1019) 'atom'
      do ia=1,n
        ir=ixres(ia)
        if (irrix(ia) .gt. 0) write (iout,1011)
     -    line(index(ia))(inamcol1:inamcol2),ia,
     -    line(index(ifres(ir)))(irescol1:irescol2),ixresno(ir),
     -    ixsegno(ir),float(irrix(ia))/float(nframe)
      end do
      write (iout,1019) 'residue'
      do ir=1,numres
        if (itemp1(ir) .gt. 0) write (iout,1010)
     -    line(index(ifres(ir)))(irescol1:irescol2),ixresno(ir),
     -    ixsegno(ir),float(itemp1(ir))/float(nframe)
      end do
      return
1001  format(2i8)
1002  format(' === List of ',i4,1x,a,' bonds filtered out:')
1003  format(' Number of ',a,' bonds left after filtering=',i5)
1004  format(/,' === List of ',i4,1x,a,' bonds left after filtering:')
1005  format(1x,a,'#',i5,' (',i5,')',a,' % formed=',f6.1,' %',a,' <d>=',
     =  f5.2)
1006  format(' Total number of ',a,' bonds per frame',a,f6.1)
1007  format(' # of bonds occuring more than ',i3,'% of time:',i4)
1008  format(' Lines ending with AA describe anchor-anchor bonds')
1009  format(' Labels BB,SS,BS,SB refer to backbone-backbone, ',
     -  'sidechain-sidechain,',/,
     -  ' backbone-sidechain, sidechain-backbone bonds, resp.')
1010  format(' Residue ',a,i5,' C/S',i2,': <Nbond>/frame=',f6.3)
1011  format(' Atom ',a,i7,' Residue ',a,i5,' C/S',i2,
     -  ': <Nbond>/frame=',f6.3)
1012  format(' List of contact atom pair indices are written to file',/,
     -  1x,a)

1019  format(/,' === Number of bonds formed per frame for each ',a)
      end
