      subroutine makeoniom(inpcrdtyp,inpfile,outfiletmp,namleni,
     -  namlentmp,n,nslt,naslv,islvw,c,index,iatnum,ifchrg,indexn,
     -  indexo,indexs,mmtype,isegno,ixres,ih,ibnd,ifgtyp,nneig,nneiga,
     -  nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,molresflag,
     -  hblimfac,angmin,ipotcol1,ipotcol2,iqcol1,iqcol2,irescol1,
     -  irescol2,nrescol,inamcol1,inamcol2,nnamcol,maxrepconf,innlist,
     -  nconfig,line,radtodeg,maxng,maxrsd,maxrec)
      character*200 inpfile,outfiletmp
      character* 132 line(maxrec),ansline
      dimension c(3,n),index(n),iatnum(n),ifchrg(n),indexn(n),indexo(n),
     -  indexs(n),mmtype(n),isegno(n),ixres(n),ih(n),ibnd(maxng,maxrec),
     -  ifgtyp(maxrec),nneig(n),ineig(maxng,n),nhbneig(n),nneiga(n),
     -  nhneig(n),nnneig(n),ncneig(n),nsneig(n),npneig(n),
     -  molresflag(maxrsd)
      character*1 aanames1
      character*2 mmodtoamb
      character*3 aanames3
      common /atnamcon/ mmodtoamb(100),aanames1(58),aanames3(58),
     -  naanames,nnanames,nnammnames,nnames,ixwatnam
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      dimension ihml(3),qhmlsum(3),nonredlist(6)
      character*1 HMtyp,gcent,hml(3),s1
      character*2 ambtyp
      character*4 potnam
      character*8 resnam,atomnam
      character*80 prtline
      data hml /'H','M','L'/
c     Generate output filenames
      call changeext(inpfile,outfiletmp,namleni,namlentmp,'onm',3,0,0)
      nosolv=1
      if (n .gt. nslt)
     -  call askyn('Do you want to include the solvents too',39,0,-1,
     -    nosolv,0,0)
      nats=nslt
      if (nosolv .eq. 0) nats=n
      print *,'ONIOM input file: ',outfiletmp(1:namlentmp)
      call openfile(21,0,' ',1,'new',outfiletmp,namlentmp,notfnd,
     -  0,1,1,0,0)
      call quiz(HMtyp,iHMtyp,' ','ONIOM',5,'region definition mode',22,
     -  0,5,6,0)
      if (HMtyp .eq. 'l') then
        print *,'Specify list of H atoms'
        call getlist(indexn,nhigh,1,nslt,1,nslt)
        print *,'Specify list of M atoms'
        call getlist(indexo,nmid,1,nslt,1,nslt)
      else if (HMtyp .eq. 's') then
        call getint('Atom number at the center of the quantum part',45,
     -    1,1,nslt,iqcent,56)
        call getreal('Radius of the H region',22,5.0,rh,1,56)
        call getreal('Radius of the M region',22,10.0,rm,1,56)
      end if
      do i=1,3
        ansline(1:37)='Do you want to keep the '//hml(i)//' layer fixed'
        call askyn(ansline(1:37),37,1,-1,ihml(i),0,0)
      end do
      if (ipotcol1*iqcol1 .gt. 0) then
        mmcformat=1
      else
        mmcformat=2
      end if
      potnam='    '
      resnam='        '
      nunknown=0
      nunknownr=0
      qsum=0.0
      qslt=0.0
      natprt=0
      iconvfile=0
      ideoxymmc=-1
      if (inpcrdtyp .gt. ioins) then
        print *,'ERROR: input format does not have atom names'
        stop
      end if
      if (iqcent .gt. 0) then
        nhigh=0
        nmid=0
        nlow=0
        do ia=1,nats
          if (dist2(c(1,iqcent),c(1,ia)) .le. rh**2) then
            ih(ia)=1
            nhigh=nhigh+1
            indexn(nhigh)=ia
          else if (dist2(c(1,iqcent),c(1,ia)) .le. rm**2) then
            ih(ia)=3
            nmid=nmid+1
          else
            ih(ia)=5
            nlow=nlow+1
          end if
        end do
        call zeroiti(mmtype,0,nats)
        call askyn('Do you want to prevent non C-C bond break',41,1,-1,
     -    icconly,0,0)
        if (icconly .eq. 1) then
          call nnlist(nslt,nslt,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -      nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -      irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -      molresflag,hblimfac,angmin,0,ibnd,ifgtyp,isegno,ixres,
     -      maxrepconf,0,0,radtodeg,0,maxng,maxng,maxrsd,maxrec)
          nhighprev=nhigh
          nhighprev0=nhigh
          iscan=1
          do while (iscan .eq. 1)
            do iaa=1,nhighprev
              ia=indexn(iaa)
              do jaa=1,nneig(ia)
                ja=ineig(jaa,ia)
                if (ih(ja) .ne. 1 .and.
     -              (iatnum(ia) .ne. 6 .or. iatnum(ja) .ne. 6)) then
c                 Add ja to the inner list
                  if (ih(ja) .eq. 3) nmid=nmid-1
                  if (ih(ja) .eq. 5) nlow=nlow-1
                  nhigh=nhigh+1
                  indexn(nhigh)=ja
                  ih(ja)=1
                end if
              end do
            end do
            if (nhigh .eq. nhighprev) iscan=0
            nhighprev=nhigh
          end do
          print *,'Added ',nhigh-nhighprev0,' atoms to H list to avoid',
     -      ' breaking non C-C bonds'
          do iaa=1,nhigh
            ia=indexn(iaa)
            do jaa=1,nneig(ia)
              ja=ineig(jaa,ia)
              if (ih(ja) .ne. 1) mmtype(ia)=ja
            end do
          end do
        end if
        nm=0
        do ia=1,nats
          if (ih(ia) .eq. 3) then
            nm=nm+1
            indexo(nm)=ia
          end if
        end do
        if (nm .ne. nmid) print *,'PROGRAM ERROR nm NE nmid:',nm,nmid
        if (icconly .eq. 1) then
          nmidprev=nmid
          nmidprev0=nmid
          iscan=1
          do while (iscan .eq. 1)
            do iaa=1,nmidprev
              ia=indexo(iaa)
              do jaa=1,nneig(ia)
                ja=ineig(jaa,ia)
                if (ih(ja) .gt. 3 .and.
     -              (iatnum(ia) .ne. 6 .or. iatnum(ja) .ne. 6)) then
c                 Add ja to the inner list
                  nlow=nlow-1
                  nmid=nmid+1
                  indexo(nmid)=ja
                  ih(ja)=3
                end if
              end do
            end do
            if (nmid .eq. nmidprev) iscan=0
            nmidprev=nmid
          end do
          print *,'Added ',nmid-nmidprev0,' atoms to M list to avoid',
     -      ' breaking non C-C bonds'
          do iaa=1,nmid
            ia=indexo(iaa)
            do jaa=1,nneig(ia)
              ja=ineig(jaa,ia)
              if (ih(ja) .eq. 5) mmtype(ia)=ja
            end do
          end do
        end if
      else
        do ia=1,nats
          ih(ia)=5
        end do
        do ia=1,nhigh
          ih(indexn(ia))=1
        end do
        do ia=1,nmid
          ih(indexo(ia))=3
        end do
      end if
      call trnsfi(indexs,indexn,nhigh)
      call trnsfi(indexs(nhigh+1),indexo(nhigh+1),nmid)
      ia=nhigh+nmid
c     indexn,indexo,indexs contain the list of H, and M, atoms, resp.
      if (iord .eq. 1) then
        do iaa=1,nats
          if (ih(iaa) .eq. 5) then
            ia=ia+1
            indexs(ia)=iaa
          end if
        end do
        if (ia .ne. nats) print *,'PROGRAM ERROR: ia NE nats:',ia,nats
c       indexs contains the atomindices sorted by H, M, L
      else
        call indexit(indexs,1,nats,0)
      end if
      potnam='    '
      resnam='        '
      nunknown=0
      nunknownr=0
      qsum=0.0
      qslt=0.0
      natprt=0
      nfix=0
      iconvfile=0
      ideoxymmc=-1
      call askyn('Do you want to reorder by ONIOM type',36,1,-1,iord,0,
     -  0)
      if (nconfig .eq. 1)
     -  call getconvdat(2,5,5,nconvdat,iconvtyp,n,line,index,
     -    inamcol1,inamcol2,nnamcol,irescol1,irescol2,igrpinfo,maxrec)
      call zeroit(qhmlsum,3)
      do iaa=1,nats
        ia=indexs(iaa)
        atomnam(1:nnamcol)=line(index(ia))(inamcol1:inamcol2)
        if (iqcol1 .gt. 0)
     -    call readreal(line(index(ia)),iqcol1,iqcol2,qslt)
        resnam(1:nrescol)=line(index(ia))(irescol1:irescol2)
        if (nrescol .gt. 3) call blankout(resnam,4,nrescol)
        if (nrescol .gt. 4) call leftadjustn(resnam,resnam,nrescol)
        if (mmcformat .eq. 1) then
          if (inpcrdtyp .eq. ioins) then
            ambtyp=line(index(ia))(ipotcol1:ipotcol2)
          else if (inpcrdtyp .eq. iommod) then
            call readint(line(index(ia)),ipotcol1,ipotcol2,mmodtyp,4,1,
     -        irerr)
            ambtyp=mmodtoamb(mmodtyp)
          end if
          potnam(1:2)=ambtyp
        else if (mmcformat .eq. 2) then
          call PDBtommc(resnam,atomnam,potnam,qslt,ia,gcent,igr,
     -      ideoxymmc,iconvfile,nunknown,nunknownr)
        end if
c_C-CT-q   3f13.7,2x,a1
c          ^ col 21
        layer=(ih(ia)+1)/2
        s1=hml(layer)
        qhmlsum(layer)=qhmlsum(layer)+qslt
        ncol=62
        write (prtline,1266) (c(k,ia),k=1,3),s1
        if (mmtype(ia) .ne. 0) then
          prtline(63:64)=' H'
          write (prtline(65:70),1011) indexs(mmtype(ia))
          ncol=ncol+8
        end if
        ian=max0(1,ianum(atomnam,0,nnamcol))
        prtline(2:3)=iatnm2(ian)
        ic=4
        if (prtline(3:3) .eq. ' ') ic=ic-1
        prtline(ic:ic)='-'
        prtline(ic+1:ic+4)=potnam
        icn=ic+4
        do while (ic .le. icn .and. prtline(ic:ic) .ne. ' ')
          ic=ic+1
        end do
        prtline(ic:ic)='-'
        if (qslt .lt. 0.0) then
          write (prtline(ic+1:ic+7),1267) qslt
          ic=ic+7
        else
          write (prtline(ic+1:ic+6),1268) qslt
          ic=ic+6
        end if
        qsum=qsum+qslt
        layer=(ih(ia)+1)/2
        if (ihml(layer) .eq. 1) then
          nfix=nfix+1
          prtline(ic+2:ic+3)='-1'
        end if
        natprt=natprt+1
        write (21,1010) prtline(1:ncol)
      end do
      call askyn('Do you want to generate topology input too',42,
     -  1,-1,itop,0,0)
      if (itop .eq. 1) then
        call nnlist(nslt,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -    nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -    molresflag,hblimfac,angmin,0,ibnd,ifgtyp,isegno,ixres,
     -    maxrepconf,0,0,radtodeg,0,maxng,maxng,maxrsd,maxrec)
        call bondord(iatnum,mmtype,n,nneig,ineig,nhneig,ibnd,maxng,c,
     -    index,ncneig,nsneig,inamcol1,inamcol2,irescol1,irescol2,line,
     -    nconfig,maxrepconf,maxrec)
        if (nosolv .eq. 0) then
          call askyn('Do you want to add H-H bonds to the solvent',43,
     -      1,-1,ihh,0,0)
          if (ihh .eq. 1) then
            nw=(n-nslt)/naslv
            do iw=1,nw
              incr=nslt+(iw-1)*naslv
              do ia1=incr+1,incr+naslv-1
                atomnam(1:nnamcol)=line(index(ia1))(inamcol1:inamcol2)
                ian1=max0(1,ianum(atomnam,0,nnamcol))
                do ia2=ia1+1,incr+naslv
                  atomnam(1:nnamcol)=line(index(ia2))(inamcol1:inamcol2)
                  ian2=max0(1,ianum(atomnam,0,nnamcol))
                  if (ian1 .eq. 1 .and. ian2 .eq. 1) then
                    nneig(ia1)=nneig(ia1)+1
                    nneig(ia2)=nneig(ia2)+1
                    ineig(nneig(ia1),ia1)=ia2
                    ineig(nneig(ia2),ia2)=ia1
                    ibnd(nneig(ia1),ia1)=1
                    ibnd(nneig(ia2),ia2)=1
                  end if
                end do
              end do
            end do
          end if
        end if
        write (21,1010)
        do iaa=1,n
          ia=indexs(iaa)
          nn=nneig(ia)
          if (nn .gt. 6) then
            write (6,1270) ia,nn,(ineig(indexs(ja),ia),ja=1,nn)
            nn=6
          else
            nng=0
            do in=1,nn
              if (ineig(indexs(in),ia) .gt. ia) then
                nng=nng+1
                nonredlist(nng)=indexs(in)
              end if
            end do
            write (21,1260) ia,
     -      (ineig(nonredlist(ja),ia),ibnd(nonredlist(ja),ia),ja=1,nng)
          end if
        end do
      end if
      close(21)
      write (6,1269) outfiletmp(1:namlentmp),
     -  nhigh,nmid,nlow,qhmlsum,qsum,nunknown,natprt,nfix
      call trnsfi(isegno,ih,n)
      return
1010  format(a)
1011  format(i6)
1260  format(i5,6(i6,i3,'.0'))
1266  format(20x,3f13.7,2x,a1)
1267  format(f7.4)
1268  format(f6.4)
1269  format(' Oniom input file is written to file ',/,7x,a,/
     -  ' Number of atoms in the high, middle and low regions:',3i7,/,
     -  ' Charge sums: H=',f8.4,' M=',f8.4,' L=',f8.4,' Total=',f8.4,
     -  /,' Number of unassigned atoms=',i5,/,
     -  ' Total number of atoms printed=',i6,/,
     -  ' Number of atoms kept fixed=',i6)
1270  format(' ERROR: Atom ',i5,' has more than 6 neighbours - ',
     -  '(format has to be changed):',/,(10i8))
      end
