      subroutine convtraj(ioutrajtyp,inptrajtyp,mmctrajtyp,inpcrdtyp,n,
     -  nslt,naslv,islvw,iatnum,ifchrg,innlist,ih,c,co,c1,c2,rprox,cv,
     -  iresno,ixres,ifres,ilres,ires_ref,ifres_ref,ilres_ref,atnames,
     -  resnames,segnames,charge,isegno,isegno_ref,segid4,ixseg1,ixseg2,
     -  marker,asterisk,ntitlin,title,ireseq,iresnrestart,iresidrestart,
     -  inamcol1,inamcol2,irescol1,irescol2,iseqncol1,iseqncol2,is1,is2,
     -  iqcol1,iqcol2,idcol,nneig,nneiga,nhbneig,ineig,nhneig,nnneig,
     -  ncneig,nsneig,npneig,numres,numslv,resnamslv,line,blankline,
     -  mmtype,ibnd,index,indexn,indexo,indexs,indexdel,index1,index2,
     -  indexsup,atw,edge,ifree,namin,namout,resin,resout,icntrl,
     -  icntrlw,iha,molresflag,altcol,inscol,ninsres,atwtemp,temp,idupl,
     -  itemp1,hblimfac,angmin,iconfsel,numsel,nextconfsel,molsltlim,
     -  nmolslt,nmolsltnoion,nooptquiz,minresflag,iqspaceask,
     -  isuperimpask,iwriteatsym,radtodeg,tofac,irepcall,ipbcinp,
     -  iwriteonly,lentest,maxconfsel,maxrepconf,pi,maxrsd,maxng,maxrec)
      character*1 asterisk,altcol(maxrec),inscol(maxrec)
      character*4 segnames(maxrsd),segid4(maxrsd)
      character* 132 line(maxrec),blankline
      character*80 title,lineinp,trtitle(32)
      dimension iconfsel(maxconfsel)
      dimension icntrl(20),icntrlw(20),molsltlim(3,maxrsd)
      dimension nneig(maxrec),ineig(maxng,maxrec),iatnum(maxrec),
     -  ifchrg(maxrec),c(3,maxrec),co(maxrec),c1(3,maxrec),c2(maxrec),
     -  rprox(maxrec),cv(maxrec),charge(maxrec),nhbneig(maxrec),
     -  nneiga(maxrec),nhneig(maxrec),nnneig(maxrec),ncneig(maxrec),
     -  nsneig(maxrec),npneig(maxrec),isegno(maxrec),isegno_ref(maxrec),
     -  ixseg1(maxrsd),ixseg2(maxrsd),iresno(maxrec),ixres(maxrec),
     -  ifres(maxrec),ilres(maxrec),ires_ref(maxrec),ifres_ref(maxrec),
     -  ilres_ref(maxrec),mmtype(maxrec),ibnd(maxng,maxrec),
     -  index(maxrec),indexn(maxrec),indexo(maxrec),indexs(maxrec),
     -  indexdel(maxrec),index1(maxrec),index2(maxrec),indexsup(maxrec),
     -  atw(maxrec),temp(maxrec),edge(3),ifree(maxrec),ih(maxrec),
     -  molresflag(maxrsd),atwtemp(maxrec),idupl(maxrsd),itemp1(maxrec)
      character*4 namin(maxrec),namout(maxrec),
     -  resin(maxrec),resout(maxrec)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      common /pbcrotmat/ torot_ac(3,3),torot_ca(3,3),tofac_ac,tofac_ca
      character*200 trajnam,outfile
      common /trajname/ trajnam,outfile,ltrajnam,namleno,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      character*8 resnamslv,atnames(maxrec),resnames(maxrsd)
      character*6 marker(16)
      character*1 xyz
      common /axislab/ xyz(3)
      character*11 trajformatname
      common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
      dimension edgexyz(3),trajrot(3,3),shiftpbc(3),c0(3),crm(3),
     -  rot(3,3)
      character*200 reffile
      character*1 ans,seqtyp,reftyp,chainid,chainid_prev
      character*4 resnam4,atnam4,resnam4n,atnam4n,segid,segid_prev
      character*8 resold,namold
      data inpt /70/,iout /71/,namleno /0/,isort /0/,ishiftpbc /0/,
     -  ishiftocrn /0/,ishiftocnt /0/,ifirst /1/, ilast /999999/,
     -  increment /1/,edgexyz /3*0.0/,iacent /0/,itrajsuperimp /0/
c     print *,'CONVTRAJ ioutrajtyp,inptrajtyp,inpcrdtyp,n=',
c    -  ioutrajtyp,inptrajtyp,inpcrdtyp,n
      iusepbc=0
      icellalt=0
      etot=0.0
      call indexit(indexs,1,n,0)
      nnamcol=inamcol2-inamcol1+1
      nrescol=irescol2-irescol1+1
      natsdel=0
      do ia=1,nslt
        atwtemp(ia)=atw(ia)
        if (molresflag(ixres(ia)) .gt. minresflag) atwtemp(ia)=0.0
      end do
c     Conversion - set up sorting
      ans='u'
103   if (nooptquiz .eq. 0) call quiz(ans,isort,'u',' ',0,
     -  'output atom order',17,0,5,6,82)
      nomatch=0
      if (ans .eq. 'u') then
c       Leave atom order unchanged
        ireorder=0
        if (inptrajtyp .ne. ioutrajtyp .and. iwriteonly .eq. 1)
     -    write (6,1014)
        call indexit(indexs,1,n,0)
        nwrmax=n
        iedit=0
        if (nooptquiz .eq. 0)
     -    call askyn('Do you want to select atoms to write',36,1,-1,
     -      iedit,00,0)
        if (iedit .gt. 0) then
          call select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
     -      is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,irescol1,
     -      irescol2,iqcol1,iqcol2,charge,iatnum,nneig,ineig,indexdel,
     -      6,maxng,maxrec)
          natsdel=0
          nwrmax=0
          do ia=1,n
            if (indexdel(ia) .eq. 1) then
              indexs(ia)=0
              natsdel=natsdel+1
            else
              nwrmax=ia
            end if
          end do
          if (natsdel .gt. 0) then
            ireorder=1
            write (6,1001) natsdel,n-natsdel,nwrmax
          else
            write (6,*) 'NOTE: all atoms are kept'
          end if
        end if
      else if (ans .eq. 'r') then
        ireorder=1
        nwrmax=n
c       Read index file
        call quiz(seqtyp,ieqtyp,'s',' ',0,'sequence file type',18,0,5,6,
     -    119)
        call zeroiti(indexn,0,n)
        lreffile=0
        if (seqtyp .eq. 's') then
c         Open the .sort file
          call openfile(69,0,'index',5,'old',reffile,lreffile,
     -      notfnd,1,1,1,0,0)
          if (notfnd .gt. 0) go to 103
          do ia=1,n
            read(69,1007,end=300) indexs(ia)
            if (indexs(ia) .eq. 0) then
              nomatch=nomatch+1
              write (6,1002) ia
            else
              indexn(indexs(ia))=1
            end if
          end do
        else if (seqtyp .eq. 'p') then
c         Read new order from a PDB file' ATOM records
c         Open the .pdb file
          call openfile(69,0,'index PDB',9,'old',reffile,lreffile,
     -      notfnd,1,1,1,0,0)
          if (notfnd .gt. 0) go to 103
          ia=0
          do while (.true.)
            read(69,1000,end=410) lineinp
            if (lineinp(1:4) .eq. 'ATOM' .or.
     -          lineinp(1:6) .eq. 'HETATM') then
              ia=ia+1
              read (lineinp(7:11),*) indexs(ia)
              if (indexs(ia) .eq. 0) then
                nomatch=nomatch+1
                write (6,1002) ia
              else
                indexn(indexs(ia))=ia
              end if
            end if
          end do
410       if (ia .eq. 0) then
            write (6,1024) reffile(1:lreffile),'PDB'
            go to 103
          end if
          print *,'Read ',ia,' atom records'
        else if (seqtyp .eq. 'c' .or. seqtyp .eq. 'l') then
c         Open the Charmm .CRD file
          call openfile(69,0,'index CRD',9,'old',reffile,lreffile,
     -      notfnd,1,1,1,0,0)
          if (notfnd .gt. 0) go to 103
          nl=0
          do while (lineinp(1:1) .eq. '*' .or. nl .eq. 0)
            read(69,1000,end=500) lineinp
            nl=nl+1
          end do
          if (nl .eq. 1) then
            write (6,1024) reffile(1:lreffile),'Charmm'
            go to 103
          end if
500       if (seqtyp .eq. 'c') read (lineinp(1:5),*) natread
          if (seqtyp .eq. 'l') read (lineinp(1:10),*) natread
          print *,'Reading ',natread,' atom reacords'
          do ia=1,natread
            read(69,1000,end=400) lineinp
            if (seqtyp .eq. 'c') read (lineinp(1:5),*) indexs(ia)
            if (seqtyp .eq. 'l') read (lineinp(1:10),*) indexs(ia)
            if (indexs(ia) .eq. 0) then
              nomatch=nomatch+1
              write (6,1002) ia
            else
              indexn(indexs(ia))=ia
            end if
          end do
        end if
400     nz=0
        do ia=1,n
          if (indexn(ia) .eq. 0) then
            nz=nz+1
            write (6,1005) ia
          end if
        end do
        close (69)
        if (nomatch .gt. 0) write (6,1006) nomatch,reffile(1:lreffile)
        if (nz .gt .0) then
          write (6,1004) nz,reffile(1:lreffile)
          stop
        end if
      else if (ioutrajtyp .gt. 3) then
        print *,'Sorry, this option currently only works with PDB or ',
     -    'Charmm or binary MMC files'
        stop
      else
c       Make a match
        ireorder=1
        nwrmax=n
c       Read output structure for atom and residue name info
        print *,'Specify a structure with the output ordering'
        lreffile=0
        call openfile(69,0,'reference order',15,'old',reffile,lreffile,
     -    notfnd,1,1,1,0,0)
        if (reffile(lreffile-3:lreffile) .eq. '.pdb') then
          ireftyp=iobpdb
          print *,'PDB format assumed'
        else if (reffile(lreffile-3:lreffile) .eq. '.CRD') then
          ireftyp=iocha
          print *,'Charmm CRD format assumed'
        else
          call quiz(reftyp,ireftyp,' ','template',8,
     -      'file form (PDB or CRD)',22,0,5,6,0)
        end if
        if (ireftyp .eq. iocha .or. ireftyp .eq. iochaex) then
c         Charmm .CRD
          nccol=4
          lineinp(1:1)='*'
          segid_prev='    '
          nsegm=0
          do while (lineinp(1:1) .eq. '*')
            read (69,1000,end=1069) lineinp
          end do
          if (ireftyp .eq. iocha) call readint(lineinp,1,5,nref,4,1,
     -      irerr)
          if (ireftyp .eq. iochaex) call readint(lineinp,1,10,nref,4,1,
     -      irerr)
          do ia=1,nref
            read (69,1000,end=1069) lineinp
            resout(ia)='    '
            resout(ia)(1:nrescol)=
     -        lineinp(iirescol(1,ireftyp):iirescol(2,ireftyp))
            namout(ia)=lineinp(iinamcol(1,ireftyp):iinamcol(2,ireftyp))
            read (lineinp(iiresncol(1,ireftyp):iiresncol(2,ireftyp)),*)
     -        ires_ref(ia)
            segid='    '
            segid=lineinp(iisegcol(1,ireftyp):iisegcol(2,ireftyp))
            if (segid .ne. segid_prev) then
              nsegm=nsegm+1
              segid_prev=segid
            end if
            isegno_ref(ia)=nsegm
          end do
        else
c         PDB
          nccol=1
          nref=0
          chainid_prev=' '
          nsegm=0
          do while (.true.)
            read (69,1000,end=101) lineinp
            if (lineinp(1:4) .eq. 'ATOM' .or.
     -          (iha .eq. 1 .and. lineinp(1:6) .eq. 'HETATM')) then
              nref=nref+1
              resout(nref)='    '
              resout(nref)(1:nrescol)=
     -          lineinp(iirescol(1,ireftyp):iirescol(2,ireftyp))
              namout(nref)=
     -          lineinp(iinamcol(1,ireftyp):iinamcol(2,ireftyp))
              read (lineinp(iiresncol(1,ireftyp):iiresncol(2,ireftyp)),
     -          *) ires_ref(nref)
              chainid=lineinp(iisegcol(1,ireftyp):iisegcol(2,ireftyp))
              if (chainid .ne. chainid_prev) then
                nsegm=nsegm+1
                chainid_prev=chainid
              end if
              isegno_ref(nref)=nsegm
            end if
          end do
        end if
101     close (69)
        if (nref .gt. n) then
          write (6,1018) nref,n
          stop
        end if
        call getdupindex(nsegm,segid4,ixseg2)
        print *,'Establishing the input and output name conventions'
        call initnamconv(noconv)
c       Convert input names to output's convention
        nrecdel=0
        nrch=0
        nach=0
        isegno_prev=0
        nsegm=0
        do ia=1,n
          if (isegno(ia) .ne. isegno_prev) then
            nsegm=nsegm+1
            segid4(nsegm)='    '
            segid4(nsegm)=line(index(ia))
     -        (iisegcol(1,inpcrdtyp):iisegcol(2,inpcrdtyp))
            isegno_prev=isegno(ia)
          end if
          resold='     '
          resold(1:nrescol)=line(index(ia))(irescol1:irescol2)
          namold='     '
          namold(1:nnamcol)=line(index(ia))(inamcol1:inamcol2)
          atnam4=namold(1:4)
          if (namold(5:5) .ne. ' ') atnam4=namold(2:5)
          resnam4=resold(1:4)
          if (resold(5:5) .ne. ' ') resnam4=resold(2:5)
          if (noconv .eq. 0) then
            call namconv(nrescol,resnam4,atnam4,resnam4n,atnam4n,nrch,
     -        nach,line(index(ia)),idcol,nrecdel)
            resin(ia)=resnam4n
            namin(ia)=atnam4n
          else
           resin(ia)=resnam4
           namin(ia)=atnam4
          end if
        end do
        call getdupindex(nsegm,segid4,ixseg1)
c       cv and rprox will not be used in maketraj/writeconf calls
        call residue_contig(n,iresno,isegno,index1,ixseg1,indexo,
     -    indexn,indexs,indexdel,ih,maxrsd,maxrec)
        call set_res_lim(iresno,n,ifres,ilres,nres,index1,1,maxrsd,
     -    maxrec)
c        do ir=1,nres
c          do ia=ifres(ir),ilres(ir)
c            write (77,8751) ir,ia,index1(ia),resin(index1(ia))
c8751        format(' ir=',i4,' ia=',i6,' index1=',i6,' resin=',a)
c          end do
c        end do
        call residue_contig(nref,ires_ref,isegno_ref,index2,ixseg2,
     -    indexo,indexn,indexs,indexdel,ih,maxrsd,maxrec)
        call set_res_lim(ires_ref,nref,ifres_ref,ilres_ref,nres_ref,
     -    index2,1,maxrsd,maxrec)
c        do ir=1,nres_ref
c          do ia=ifres_ref(ir),ilres_ref(ir)
c            write (78,8752) ir,ia,index2(ia),resout(index2(ia))
c8752        format(' ir=',i4,' ia=',i6,' index2=',i6,' resout=',a)
c          end do
c        end do
c       Set up matching array
        reffile=reffile(1:lreffile)//'.sort'
        lreffile=lreffile+5
        call openfile(69,0,'order',5,'new',reffile,lreffile,notfnd,
     -    0,1,1,0,0)
        call zeroiti(indexn,0,n)
        call zeroiti(indexs,0,n)
        ir_ref_prev=0
        ir_ref=1
        ifres_ref(nres_ref+1)=ifres_ref(nres_ref)
        do ir=1,nres
c         Find the next matching residue in the template
          ir_ref=ir_ref_prev+1
c         print *,'Start checking ir,ir_ref=',ir,ir_ref
          do while (resin(index1(ifres(ir))) .ne.
     -              resout(index2(ifres_ref(ir_ref)))
     -      .and. ir_ref .le. nres_ref)
            ir_ref=ir_ref+1
          end do
          if (resin(index1(ifres(ir))) .ne.
     -        resout(index2(ifres_ref(ir_ref)))) go to 210
          ir_ref_prev=ir_ref
          do iaa=ifres(ir),ilres(ir)
            ia=index1(iaa)
            do jaa=ifres_ref(ir_ref),ilres_ref(ir_ref)
              ja=index2(jaa)
              if (indexn(ja) .eq. 0) then
                if (namin(ia) .eq. namout(ja)) then
                  indexn(ja)=ia
                  go to 200
                end if
              end if
            end do
            nomatch=nomatch+1
c           if (ia .le. nref) write (6,1020) ia,resin(ia),namin(ia)
200         continue
          end do
        end do
210     if (ir_ref .gt. nres_ref) write (6,1010) ir
        do ia=1,n
          indexs(ia)=indexn(ia)
        end do
        noorder=0
        do ia=2,nref
          if (indexs(ia) .ne. indexs(ia)+1) noorder=noorder+1
        end do
        if (noorder .eq. 0) then
          print *,'NOTE: matching kept the original atom order'
          if (n .eq. nref) ireorder=0
        end if
        do ia=1,n
          if (indexs(ia) .gt. 0) then
            lineinp=line(index(indexs(ia)))(1:80)
            resold='     '
            resold(1:nrescol)=lineinp(irescol1:irescol2)
            namold='     '
            namold(1:nnamcol)=lineinp(inamcol1:inamcol2)
            write (69,1021) ia,resout(ia),namout(ia),
     -        resold, namold,indexs(ia)
          else
            write (69,1021) ia,resout(ia),namout(ia),
     -        '****','****',indexs(ia)
          end if
        end do
        write (6,1023) reffile(1:lreffile)
        close (69)
        if (nomatch .gt. 0) then
          write (6,1008) nomatch,reffile(1:lreffile)
          print *,'You may try to modify ',reffile(1:lreffile),
     -      ' and use it to specify the new order'
          call askyn('Do you want to try matching again',33,1,1,newm,0,
     -      0)
          if (newm .gt. 0) go to 103
          print *,'Output trajectory will contain only the matched ',
     -      'atoms'
        else
          print *,'All atoms matched successfully'
        end if
      end if
      call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
     -  inpcrdtyp,ifirst,ilast,increment,maxconf,
     -  ninconf,noutconf,natom,nfreeat,ifree,icntrl,1,mmctrajtyp,
     -  trajnam,ltrajnam,'input trajectory',16,iconfsel,numsel,
     -  1-irepcall,0,0,0,icellfound,notfnd,0,0,lentest_ok,0,ianaltyp,
     -  mx2d,maxconfsel,maxrec)
102   if (irepcall .eq. 0) then
        if (icellfound .eq. 1) then
          write (6,1012)
          call askyn('Do you want to override the cell size/shape',43,
     -      1,-1,ipbcinp,0,0)
        end if
        if (ipbcinp .eq. 1 .or. icellfound .eq. 0) then
          call setpbccell('',0,edge,edge_gen,cell,ncell,cellalt,
     -      ixyzhex,npbc,ioppbc,iusepbc,vol,nw,rinscr,rcirc,0)
          call trnsfr(cell0,cell,3*ncell)
        end if
      end if
      irot_to=0
      if (icellfound .eq. 1 .and. ipbcinp .eq. 0) then
c       Generate cell info from trajectory cell sizes
        if (irepcall .eq. 0) call pbctype(ioppbc,npbc,ixyzhex,0)
        if (inptrajtyp .eq. 1) then
          if (ioppbc .eq. 5) then
c           Skewed hexagon
            edge(1)=xtlabc0(ixcrd(3))
            edge(2)=xtlabc0(ixcrd(2))
            edge(3)=xtlabc0(ixcrd(1))
            if (ioutrajtyp .ne. inptrajtyp) print *,'WARNING: Skewed ',
     -        'hexagon cell parameters will not be properly converted'
          else
            do k=1,3
              edge(k)=xtlabc0(ixcrd(k))
            end do
          end if
        else
          call trnsfr(edge,box0,3)
        end if
        tofac=1.0
        ichangeconv=0
        if (ioppbc .eq. 6) then
          call askyn(
     -      'Is the input PBC TO cell in the Amber convention',48,1,-1,
     -      ichangeconv,0,0)
          if (ichangeconv .eq. 1) then
c           TO,  Amber to Charmm
            irot_to=-1
            tofac=tofac_ac
            write (6,1013) 'Charmm'
          end if
        else if (ioppbc .eq. 7) then
          call askyn(
     -      'Is the input PBC TO cell in the Charmm convention',49,1,-1,
     -      ichangeconv,0,0)
          if (ichangeconv .eq. 1) then
c           TO, Charmm to Amber
            irot_to=1
            tofac=tofac_ca
            write (6,1013) 'Amber'
          end if
        end if
        if (ichangeconv .eq. 0 .and. (ioppbc .eq. 6 .or. ioppbc .eq. 7))
     -    write (6,*) 'NOTE: TO cell orientation read will be kept'
        if (ioppbc .eq. 6) then
          do k=1,3
            edgexyz(k)=edge(k)/2.0
          end do
        else if (ioppbc .eq. 7) then
          do k=1,3
            edgexyz(k)=tofac_ac*edge(k)/2.0
          end do
        else
          call trnsfr(edgexyz,edge,3)
        end if
c       Recreate the cell from the first frame's cell size
        call crorgn(edgexyz,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -    ixyzhex,rinscr,rcirc)
        call trnsfr(cell0,cell,3*ncell)
        write (6,1016) edgexyz
        if (ioppbc .gt. 0) iusepbc=1
      end if
      if (ioppbc .eq. 6 .and. inptrajtyp .eq. 2) then
        print *,'Truncated Octahedron type conflicts with input ',
     -    'trajectory type'
        call askyn('Do you want to continue',23,1,-1,icont,20,0)
        if (icont .eq. 0) go to 102
      end if
      itrajrot=0
      call unitmat(trajrot)
      if (nooptquiz .eq. 0) then
        if (irot_to .ne. 0) write (6,1015)
        call askyn(
     -    'Do you want to rotate each snapshot of the trajectory',
     -    53,1,-1,itrajrot,97,0)
        if (itrajrot .gt. 0) then
          itrajsuperimp=0
          call genrot(trajrot,pi,iax,angle)
        else if (isuperimpask .eq. 1) then
          call askyn(
     -   'Do you want to superimpose each frame to the input structure',
     -      60,1,-1,itrajsuperimp,0,6)
          if (itrajsuperimp .gt. 0) then
            write (6,1025)
            call askyn('Do you want to select atoms for overlay',39,
     -          1,-1,ieditoverlay,69,0)
            if (ieditoverlay .gt. 0) then
              call select(line,nrecdel,idcol,asterisk,n,nslt,index,
     -          ixres,is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,
     -          irescol1,irescol2,iqcol1,iqcol2,charge,iatnum,
     -          nneig,ineig,indexdel,6,maxng,maxrec)
c             indexdel: 0 for atoms to use for the overlay, 1 for the rest
c             print *,'nfinalov=',nfinalov
c             write (6,9173) (indexdel(i),i=1,nfinalov)
c9173          format(80i1)
              call masktolist(indexsup,indexdel,n,nfinalov,0)
            else
c             Get a condensed list of atoms selected to write
              ndel=0
              do ia=1,nwrmax
                if (indexs(ia) .eq. 0) then
                  ndel=ndel+1
                else
                  indexsup(ia-ndel)=indexs(ia)
                end if
              end do
              nfinalov=nwrmax-ndel
c             print *,'nfinalov=',nfinalov,' nwrmax=',nwrmax,' nd=',ndel
            end if
          endif
        end if
        ipbcreset=0
        if (itrajsuperimp .eq. 1) then
          call trnsfr(co,c,3*n)
        else if (iusepbc .gt. 0) then
          call askyn(
     -      'Do you want to reset all molecules into the PBC cell',52,
     -      -1,1,ipbcreset,0,0)
        else
          call setpbccell(
     -      'Do you want to reset all molecules into the PBC cell',52,
     -      edge,edge_gen,cell,ncell,cellalt,ixyzhex,npbc,
     -      ioppbc,ipbcreset,vol,nw,rinscr,rcirc,0)
          call trnsfr(cell0,cell,3*ncell)
        end if
        if (ipbcreset .gt. 0 .and.
     -      (ioppbc .eq. 6 .or. ioppbc .eq. 7)) then
          call askyn('Do you want to try both TO orientations',39, -1,1,
     -      icellalt,137,0)
        end if
      end if
      if (irot_to .eq. -1) then
        call matprod(trajrot,torot_ac,trajrot)
        itrajrot=1
      else if (irot_to .eq. 1) then
        itrajrot=1
        call matprod(trajrot,torot_ca,trajrot)
        itrajrot=1
      end if
      imcenter=0
      if (ipbcreset .gt. 0) then
        call setpbcdim(ioppbc,ixyzhex,ixyzexcld,ixyzincld,xyz)
        call setrepats(nmolslt,molsltlim,nslt,nneig,ineig,
     -    molresflag,minresflag,indexn,indexo,maxng,maxrsd)
        call quiz(ans,iansrun,'o',' ',0,
     -    'origin of the modified cell',27,0,5,6,0)
        if (ans .eq. 'o') then
          ishiftpbc=0
        else if (ans .eq. 'c') then
          call getxyz('',0,' component',10,999999.0,shiftpbc,0,0)
          ishiftpbc=1
        else if (ans .eq. 'a') then
          call getint('Atom index of the atom to be at the center',42,
     -      1,1,n,iacent,0)
          ishiftpbc=2
        end if
        if (ans .ne. 'a' .and. nmolslt .gt. 1) then
          call askyn(
     -      'Do you want to specify the solute molecule in the center',
     -        56,1,-1,imcenterset,0,0)
          if (imcenterset .eq. 1) call getint('Solute molecule number',
     -      22,1,1,nmolslt,imcenter,0)
        end if
        call prtcell(ioppbc,edge,edge_gen,r,vol,nw,1)
      else if (ioppbc .eq. 1 .or. ioppbc .eq. 2) then
        idef=+1
        if (inptrajtyp .eq. ioutrajtyp) idef=-1
        if (ioutrajtyp .eq. 2 .and. nooptquiz .eq. 0) then
c         Option to shift center to box corner
          call askyn(
     -      'Do you want to shift the origin to the box corner',
     -      49,1,idef,ishiftocrn,0,0)
          if (ishiftocrn .gt. 0) shfac=+1.0
        end if
        if (ishiftocrn .eq. 0 .and. inptrajtyp .eq. 2) then
c         Option to shift center from box corner
          call askyn(
     -      'Do you want to shift the origin to the box center',
     -      49,1,idef,ishiftocnt,0,0)
          if (ishiftocnt .gt. 0) shfac=-1.0
        end if
      end if
      do k=1,3
        box(k)=edge(k)
        xtlabc(ixcrd(k))=edge(k)
        xtlabc(ixang(k))=90.0
        edgexyz(k)=edge(k)/2.0
      end do
      if (noboxinfoar .eq. 0 .and. ipbcreset .eq. 1)
     -  write (6,1009)
      if (numsel .gt. 0) nextconfsel=1
      ifttyp=iftrajtyp(ioutrajtyp)
      if (ioutrajtyp .eq. 3 .and. mmctrajtyp .eq. 1) ifttyp=2
      call openfile(iout,0,'output trajectory',17,'new',outfile,namleno,
     -  notfnd,0,ifttyp,1,0,0)
      if (ioutrajtyp .eq. 1) then
        ifound=0
        i=ntitltr-1
        do while (i .gt. 1 .and. ifound .eq. 0)
          if (trtitle(i)(1:41) .eq.
     -       'Simulaid converted the trajectory in file') then
            trtitle(i+1)=trajnam(1:ltrajnam)
            ifound=1
          end if
          i=i-2
        end do
        if (ifound .eq. 0 .and. ntitltr .lt. 30) then
          trtitle(ntitltr+1)=
     -      'Simulaid converted the trajectory in file'
          trtitle(ntitltr+2)=trajnam(1:ltrajnam)
          ntitltr=ntitltr+2
        end if
      end if
      if (ioutrajtyp .le. 2) then
        if (ioutrajtyp .eq. 2) ntitltr=1
        write (6,1026) title,(trtitle(i),i=1,ntitltr)
        if (ioutrajtyp .eq. 2) then
          call askyn('Do you want to use the structure title for both',
     -      47,1,1,irep,0,0)
          if (irep .eq. 1) trtitle(1)=title
        else
          call askyn('Do you want to add the structure title',
     -      38,1,1,irep,0,0)
          if (irep .eq. 1) then
            ntitltr=min0(32,ntitltr+1)
            trtitle(ntitltr)=title
          end if
        end if
      end if
      nntest=0
      nr50=0
      nr50t=0
      ielim=1
      neskip=0
      do while (ninconf .lt. ilast)
c       Read a conformation
        call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,nslt,
     -    trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,ltrajnam,
     -    natom,nfreeat,ifree,icntrl,c,ninconf,noutconf,
     -    increment,inpcrdtyp,ietot,etot,ifail,ifirst,ilast,iconfsel,
     -    numsel,maxrepconf,nmc,lentest,tofac,maxconf,maxconfsel,maxrec)
c       r12=sqrt(dist2(c(1,molsltlim(3,1)),c(1,molsltlim(3,2))))
c       if (r12 .gt. 40.0) then
c         write (77,*) 'NMC=',nmc,' R12=',r12
c         nr50=nr50+1
c       end if
        if (ifail .gt. 0) stop
        icsel=0
        if (numsel .eq. 0) then
          if (ninconf .ge. ifirst .and.
     -      mod(ninconf-ifirst,increment) .eq. 0) icsel=1
          ifpr=ifirst
        else
c         See if this configuration is on the list
          ifpr=iconfsel(nextconfsel)
          if (ninconf .eq. ifpr) then
            icsel=1
            nextconfsel=nextconfsel+1
          else if (ifpr .eq. 0) then
c           List exhausted
            ninconf=ilast+increment
          end if
        end if
        if (ietot .eq. 1) then
          if (ielim .eq. 1) then
            call askyn('Do you want to give an energy threshold (max)',
     -        45,1,-1,iemax,0,0)
            if (iemax .eq. 1)
     -        call getreal('Energy maximum to use',21,0.0,emax,0,0)
            ielim=0
          end if
          if (iemax .eq. 1 .and. etot .gt. emax) then
            if (icsel .eq. 1) neskip=neskip+1
            icsel=0
          end if
        end if
        if (ninconf .eq. 1 .and. ifpr .gt. 25) print *,
     -    'Searching for first configuration to write - wait ...'
        if (icsel .gt. 0) then
c         Write a (converted) conformation
          noutconf=noutconf+1
          if (nntest .eq. 0) then
            call comparetop(c,n,nneig,ineig,iatnum,innlist,nslt,
     -        naslv,cell,ncell,ioppbc,maxng,maxrec)
            nntest=1
          end if
          if (ipbcreset .eq. 1) then
c           Reset into PBC
            if (noboxinfoar .eq. 0) call updatecell(inptrajtyp,edge)
            call systemcenter(n,nmolslt,nmolsltnoion,molsltlim,c,c2,
     -        indexn,atwtemp,cell,ncell,cellalt,icellalt,ixyzexcld,
     -        ixyzincld,nslt,naslv,iacent,imcenter,imcenter_set,
     -        noutconf,maxrec,maxrsd)
c           r12=sqrt(dist2(c(1,molsltlim(3,1)),c(1,molsltlim(3,2))))
c           if (r12 .gt. 40.0) then
c             write (78,*) 'NMC=',nmc,' R12=',r12
c             nr50t=nr50t+1
c           end if
cDB            call debug_c(noutconf,line,index,nslt,cell,ncell,
cDB     -        molsltlim(1,1),molsltlim(2,1),
cDB     -        molsltlim(1,2),molsltlim(2,2),maxrec)
          end if
          if (itrajsuperimp .gt. 0) then
            call bestoverlay(nfinalov,indexsup,indexsup,co,c,atw,0.d0,
     -        c1,c2,temp,rot,c0,crm,0,0.001,0,maxrec)
            call shiftmol(c,n,crm,c1,-1.0)
            call rotate_c(c1,n,rot,c2,'OVERLAY',7)
            call shiftmol(c2,n,c0,c,+1.0)
          else if (itrajrot .gt. 0) then
            call rotate_c(c,n,trajrot,c,'CONVTRAJ',8)
          end if
          if (ishiftocrn+ishiftocnt .ne. 0)
     -      call shiftmol(c,n,edgexyz,c,shfac)
          call maketraj(noutconf,maxconf,c,c1,rprox,cv,ih,nslt,n,naslv,
     -      islvw,ntitltr,trtitle,line,altcol,inscol,ninsres,inpcrdtyp,
     -      0,inptrajtyp,ioutrajtyp,mmctrajtyp,iout,iatnum,ifchrg,
     -      innlist,iresno,ixres,atnames,resnames,segnames,charge,
     -      isegno,marker,ntitlin,title,ireseq,iresnrestart,
     -      iresidrestart,nneig,nneiga,nhbneig,ineig,nhneig,nnneig,
     -      ncneig,nsneig,npneig,numres,numslv,resnamslv,blankline,
     -      mmtype,ibnd,index,indexn,indexo,indexs,icntrlw,icntrl,ifree,
     -      molresflag,idupl,itemp1,ireorder,icellfound,natsdel,nwrmax,
     -      hblimfac,angmin,iqspaceask,keeprem,iwriteatsym,radtodeg,
     -      etot,nooptquiz,maxrepconf,maxng,mxres,maxrec)
          if (nextconfsel .gt. numsel) go to 9999
        end if
        if (ilast .ne. 999999 .and. ilast .ge. 1000) then
          ilast10=ilast/10
          if (mod(ninconf,ilast10) .eq. 0)
     -      print *,(ninconf/ilast10)*10,'% done'
        end if
      end do
      if (neskip .gt. 0) write (6,1027) emax,neskip,ninconf,noutconf
c     write (6,*) 'NR50=',nr50,' NR50t=',nr50t
      return
300   write (6,1003) ia,reffile(1:lreffile),nomatch
      return
9999  write (6,1011) noutconf
      return
1069  write (6,1019)
      stop
1000  format(a)
1001  format(' Editing deleted',i7,' atoms. Number of atoms left=',i8,
     -  /,' Highest index of kept atoms=',i8)
1002  format(' Target atom',i6,' has no original atom assigned to it')
1003  format(' ERROR: file ',a,' contained only',i7,' matches',/,
     -  8x,'Number of missing matches found=',i5)
1004  format(' ERROR:',i5,' original atoms are unassigned in file ',a)
1005  format(' Original atom',i6,' is not assigned a new index')
1006  format(' NOTE:',i5,' target atoms missing matches in file ',a,/,
     -  6x,'- output trajectory will contain only the matcehd atoms')
1007  format(59x,i7)
1008  format(' Matching failed for ',i5,' atoms',/,
     -  ' NOTE: see the file ',a,' for more info.',/,
     -  7x,'You can edit this file to complete the matches')
1009  format(' Note: Make sure that the PBC information read from ',
     -  'the trajectory',/,
     -  ' is in the same convention as the PBC size specified above')
1010  format(' WARNING: run out of template residues at the',i6,'-th',
     -  ' input residue')
1011  format(' Completed writing of the new trajectory with ',i6,
     -  ' configurations')
1012  format(' Periodic cell information was found')
1013  format(' NOTE: TO cell will be rotated to conform to the ',a,
     -  ' convention')
1014  format(' WARNING: trajectory interconversion leaves ',
     -  'the atom order unchanged')
1015  format(' NOTE: periodic cell will be rotated to follow Charmm vs',
     -  ' Amber/NAMD convention',/,7x,'You can add an additional ',
     -  'rotation.')
1016  format(' Cell parameters are set from the trajectory to',/,3f10.5)
c1017  format(' Nconf=',i8,' rc1,rcmin1=',f7.1,f8.1,' rc2,rcmin2=',
c     -  f7.1,f8.1,/,' n11-2,n21-2=',i4,3i5,' vtot1=',f10.1,' vtot2=',
c     -  f10.2,' rmsd1=',f7.1,' rmsd2=',f7.1,' isw=',i1,/)
1018  format(' ERROR: number of atoms in the reference system (',
     -  i5,')',/,8x,'is more than the input system (',i5,')')
1019  format(' ERROR: reference configuration file is incomplete')
c1020  format(' ERROR: no match for atom ',i5,' res=',a5,' name=',
c     -  a5)
1021  format(' Target: #=',i7,' R=',a5,' N=',a5,
     -       ' Orig: R=',a5,' N=',a5,' #=',i7)
1023  format(' Matching found between the two structures is written ',
     -  'to file:',/,5x,a)
1024  format(' ERROR: file ',a,' is not a ',a,' file')
1025  format(' NOTE: superimposing to the input structure contradicts ',
     -  'PBC resets.',/,' If PBC resets are needed, perform it first ',
     -  ' followed by an other run',/,' performing the superimposition',
     -  ' without declaring any',/,' mobile ion or mobile solute ',
     -  'molecule')
1026  format(' Current structure title:',/,1x,a,/,
     -  ' Current trajectory title:',(/,1x,a))
1027  format(' Number of frames filtered out with E > ',e12.5,':',i4,/,
     -  ' Number of frames read=',i6,' Number of frames written=',i6)
      end
